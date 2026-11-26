"""Additive primary-point ensemble using the unmodified accepted live producer."""
from pathlib import Path
import argparse
import ctypes as ct
import csv
import gzip
import hashlib
import io
import json
import math
import subprocess
import sys
import tarfile
import tempfile
import types

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = '6512dbc91202bb413f6bb4d188afeafcb3c69a65'
NRECORDS = 32
BATCH = 8
EXPOSURE = 1e-4
DT = 5e-8
LABEL = 'g13-pass2-primary-portable-v1|'


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def table(path, rows):
    if not rows:
        raise ValueError('Empty table')
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accepted():
    sys.path.insert(0, str(HERE.parent))
    import coarse
    import live
    return coarse, live


def mix(x):
    x = np.asarray(x, dtype=np.uint64)
    with np.errstate(over='ignore'):
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    return x ^ (x >> np.uint64(31))


def keys(record, label, step, size):
    prefix = int.from_bytes(hashlib.sha256((LABEL+label+'|'+str(record)).encode()).digest()[:8], 'little')
    return np.ascontiguousarray(mix(np.uint64(prefix) ^ (np.uint64(step)*np.uint64(262144)+np.arange(size, dtype=np.uint64))))


class Portable:
    def __init__(self):
        source = HERE/'portable.c'
        self.temporary = tempfile.TemporaryDirectory(prefix='g13-portable-')
        library = Path(self.temporary.name)/'portable.so'
        subprocess.run(['cc', '-O2', '-std=c99', '-fPIC', '-shared', '-ffp-contract=off', '-fno-fast-math', str(source), '-lm', '-o', str(library)], check=True)
        self.lib = ct.CDLL(str(library))
        i64 = np.ctypeslib.ndpointer(dtype=np.int64, flags='C_CONTIGUOUS')
        u64 = np.ctypeslib.ndpointer(dtype=np.uint64, flags='C_CONTIGUOUS')
        f64 = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        self.lib.portable_binomial.argtypes = [ct.c_int64, i64, f64, u64, i64]
        self.lib.portable_uniform.argtypes = [ct.c_int64, u64, f64]
        self.lib.portable_multinomial.argtypes = [ct.c_int64, ct.c_int64, f64, ct.c_uint64, i64]
        self.lib.portable_poisson.argtypes = [ct.c_double, ct.c_uint64]
        self.lib.portable_poisson.restype = ct.c_int64

    def uniform(self, record, label, step, size):
        out = np.empty(size)
        self.lib.portable_uniform(size, keys(record, label, step, size), out)
        return out

    def poisson(self, lam, record, interval):
        return self.lib.portable_poisson(lam, int(keys(record, 'total', interval, 1)[0]))

    def multinomial(self, n, p, record, label, step):
        p = np.ascontiguousarray(p, dtype=float)
        out = np.empty(len(p), dtype=np.int64)
        self.lib.portable_multinomial(int(n), len(p), p, int(keys(record, label, step, 1)[0]), out)
        return out

    def binomial(self, n, p, record, label, step):
        n = np.ascontiguousarray(n, dtype=np.int64)
        p = np.ascontiguousarray(np.broadcast_to(p, n.shape), dtype=float)
        out = np.empty_like(n)
        self.lib.portable_binomial(n.size, n.ravel(), p.ravel(), keys(record, label, step, n.size), out.ravel())
        return out


class CountRng:
    def __init__(self, portable, records):
        self.p, self.records = portable, records
        self.multicall = self.bincall = 0

    def multinomial(self, n, p):
        k = self.multicall
        self.multicall += 1
        r = len(self.records)
        if k < 2*r:
            interval, rep = divmod(k, r)
            return self.p.multinomial(n, p, self.records[rep], 'time', interval)
        step, rep = divmod(k-2*r, r)
        return self.p.multinomial(n, p, self.records[rep], 'site', step)

    def binomial(self, n, p):
        k = self.bincall
        self.bincall += 1
        step, rem = divmod(k, 2*len(self.records))
        rep, stage = divmod(rem, 2)
        return self.p.binomial(n, p, self.records[rep], ('species', 'sign')[stage], step)


class TotalRng:
    def __init__(self, portable, records):
        self.p, self.records, self.call = portable, records, 0

    def poisson(self, lam):
        interval, rep = divmod(self.call, len(self.records))
        self.call += 1
        return self.p.poisson(lam, self.records[rep], interval)


class MarkRng:
    def __init__(self, portable, records, burn):
        self.p, self.records, self.burn = portable, records, burn
        self.call = 0
        self.sha = [hashlib.sha256() for _ in records]
        self.psha = [hashlib.sha256() for _ in records]
        self.marks = np.zeros((len(records), 3), dtype=np.int64)
        self.low = None

    def binomial(self, n, p):
        step, stage = divmod(self.call, 2)
        self.call += 1
        out = np.empty_like(n)
        for r, record in enumerate(self.records):
            out[r] = self.p.binomial(n[r], np.broadcast_to(p, n.shape)[r], record, ('low', 'high')[stage], step)
            self.sha[r].update(out[r].astype('<i8', copy=False).tobytes())
            if stage:
                q = np.floor(np.broadcast_to(p, n.shape)[r]*2**40+.5).astype('<u8')
                self.psha[r].update(q.tobytes())
        if stage == 0:
            self.low = out
        elif step >= self.burn:
            self.marks += np.stack([self.low.sum(axis=(1, 2, 3)), (n-out).sum(axis=(1, 2, 3)), out.sum(axis=(1, 2, 3))], axis=-1)
        return out


def prepare():
    c, l = accepted()
    _, stopping = c.profiles_and_stopping()
    model = c.assemble('grains', c.DIMS[0]/8, 5000., 250000.)
    metadata = l.mode_group_velocities(model)
    eps, uz, group = l.probes(model)
    st = stopping['PdD', 2.5]
    sm = c.source_matrix(model, [c.deposition_depth(t) for t in st['tracks']])
    ep = np.array([st['tracks'][0]['components_eV'][1], 2*st['tracks'][1]['components_eV'][1]])
    weights = np.array([.85, .15])*ep
    weights /= weights.sum()
    reference = l.source_reference(model, sm, np.repeat(weights/36, 36), eps, uz, .01*138e-6*st['host_recoil_J_per_true_C'])
    eh = ep*l.QE
    eh *= st['host_recoil_J_per_true_C']/(np.dot([.85, .15], eh)/l.QE)
    receiver = l.off_beam_receiver(model)
    regions = []
    for side in (.0025, .01, float(c.DIMS[0])):
        for z0, z1 in ((0., 50e-6), (25e-6, 75e-6), (175e-6, 225e-6)):
            e, u = l.linear_region(model, side, z0, z1)
            regions.append((side, z0, z1, e, u))
    arrays = dict(omega=model['omega'], group_velocity=model['group_velocity'], eigenvectors=model['v'], eps=eps, uz=uz, group=group, eh=eh, source=sm,
                  covariance=reference[0], ecanon=np.array(reference[1]), canonical_eps=reference[2], canonical_uz=reference[3], reference_residual=np.array(reference[4]),
                  region_bounds=np.array([r[:3] for r in regions]), region_projectors=np.array([r[3:] for r in regions]), site_weights=np.array([l.beam_profile(s) for s in range(3)]))
    for i, a in enumerate(receiver):
        arrays['receiver_'+str(i)] = np.asarray(a)
    burn = round(5*np.max(1/(1/50e-6+model['group_velocity']/.2))/DT)
    steps = burn+round(EXPOSURE/DT)
    portable = Portable()
    class BeamRng:
        def __init__(self):
            self.call = 0
        def random(self, size):
            k = self.call
            self.call += 1
            return np.array([portable.uniform(r, 'beam', k, 1)[0] for r in range(size)])
        def uniform(self, lo, hi, size):
            return lo+(hi-lo)*self.random(size)
    l.REPS = NRECORDS
    arrays['time_weights'] = np.array([l.beam_timing(steps, t, BeamRng()) for t in range(3)])
    arrays['burn'] = np.array(burn)
    np.savez_compressed(HERE/'inputs.npz', **arrays)
    sources = ['g13_rogue/pass2/live.py', 'g13_rogue/pass2/coarse.py', 'g4_beam/merit/profiles.json', 'g11_cmns/collective/part1/stopping.json']
    import scipy
    dump(HERE/'source_ledger.json', dict(base=BASE, accepted_sources={p:digest(ROOT/p) for p in sources}, canonical_inputs_sha256=digest(HERE/'inputs.npz'),
        generation_runtime=dict(numpy=np.__version__, scipy=scipy.__version__), distribution_upstream='https://github.com/numpy/numpy/blob/v2.1.2/numpy/random/src/distributions/distributions.c',
        distribution_upstream_sha256='461755e7126c677b596ba0c388a33c83789e93834687db788c170d85704556b1',
        retained_modes=len(model['omega']), operator_residual=model['residual'], mass_error=model['mass_error'], reference_residual=reference[4],
        grain_side_m=c.DIMS[0]/8, basis_degree=[5,5,2], burn_steps=burn, exposure_s=EXPOSURE, dt_s=DT, independent_records=NRECORDS,
        scope=dict(energy_keV=2.5,current_A=138e-6,intrinsic_energy_lifetime_s=50e-6,energy_length_m=.2,edge_energy_reflectivity=1.,eta=.01,fmin_Hz=5000.,fmax_Hz=250000.)))
    print(json.dumps(dict(prepared=True, modes=len(model['omega']), burn_steps=burn)), flush=True)


def load_input():
    ledger = json.loads((HERE/'source_ledger.json').read_text())
    if ledger['base'] != BASE or digest(HERE/'inputs.npz') != ledger['canonical_inputs_sha256']:
        raise AssertionError('Canonical numeric input identity changed')
    for path, expected in ledger['accepted_sources'].items():
        if digest(ROOT/path) != expected:
            raise AssertionError('Accepted source identity changed: '+path)
    return dict(np.load(HERE/'inputs.npz', allow_pickle=False))


def run_arm(data, portable, records, spatial, timing, beta, gauge=False):
    c, l = accepted()
    l.REPS, l.EXPOSURE, l.DT, l.ETA = len(records), EXPOSURE, DT, .01
    model = {k:data[k] for k in ('omega', 'group_velocity')}
    reference = (data['covariance'], float(data['ecanon']), data['canonical_eps'], data['canonical_uz'], float(data['reference_residual']))
    source = dict(eh=data['eh'], s=(data['source'][:,:36], data['source'][:,36:]))
    eps, uz = data['eps'], data['uz']
    receiver = tuple(data['receiver_'+str(i)] for i in range(7))
    regions = [dict(side=b[0],z0=b[1],z1=b[2],projectors=p) for b,p in zip(data['region_bounds'],data['region_projectors'])]
    if gauge:
        signs = np.where(np.arange(len(model['omega'])) % 2, -1., 1.)
        source = dict(eh=source['eh'], s=tuple(signs[:,None]*s for s in source['s']))
        eps, uz = eps*signs, uz*signs
        both = np.tile(signs,2)
        reference = (reference[0]*both[:,None]*both[None,:], *reference[1:])
        receiver = (receiver[0]*signs,receiver[1]*signs,*receiver[2:4],receiver[4]*signs[:,None]*signs[None,:],*receiver[5:])
        regions = [dict(r,projectors=np.array(r['projectors'])*signs) for r in regions]
    count, total = CountRng(portable, records), TotalRng(portable, records)
    mark = MarkRng(portable, records, int(data['burn']))
    generators = {4104:count,8104:mark,14109:None,24109:total}
    real_np, real_profile, real_timing = l.np, l.beam_profile, l.beam_timing
    class NumpyProxy:
        random = types.SimpleNamespace(default_rng=lambda seed:generators[seed])
        def __getattr__(self,name):
            return getattr(np,name)
    l.np = NumpyProxy()
    l.beam_profile = lambda _:data['site_weights'][spatial]
    l.beam_timing = lambda steps, level, rng:data['time_weights'][timing,:,records].T
    try:
        result = l.simulate(model,source,data['group'],eps,uz,reference,receiver,regions,2.5,138e-6,50e-6,.2,beta,spatial,timing,1.)
    finally:
        l.np, l.beam_profile, l.beam_timing = real_np, real_profile, real_timing
    rows = []
    for r, record in enumerate(records):
        parents = int(result['impact']['all']['parents_by_record'][r])
        row = dict(record=record, spatial=spatial, timing=timing, beta=beta, exposure_s=EXPOSURE, parents=parents,
                   mark_low=int(mark.marks[r,0]),mark_middle=int(mark.marks[r,1]),mark_high=int(mark.marks[r,2]),
                   draw_sha256=mark.sha[r].hexdigest(), probability_sha256=mark.psha[r].hexdigest())
        for name in ('A_disturbed','C_any','C_strain_positive','C_strain_negative','C_displacement_positive','C_displacement_negative'):
            n = int(result['impact'][name]['parents_by_record'][r])
            row[name+'_parents'] = n
            row[name+'_fraction'] = n/parents
            row[name+'_rate_s'] = n/EXPOSURE
        b = l.occurrence(result['B'][:,r:r+1],DT)
        for key in ('time_fraction','complete_events','boundary_censored_events','initiation_rate_s'):
            row['B_'+key] = b[key]
        rows.append(row)
    ledger = dict(spatial=spatial,timing=timing,beta=beta,records=records,**{k:result[k] for k in ('source_work_J','incoming_host_J','marked_self_J','damping_loss_J','energy_balance_residual_J','max_budget_violation_J','capped_record_bins')})
    if abs(ledger['energy_balance_residual_J']) > 1e-9*ledger['incoming_host_J'] or ledger['max_budget_violation_J'] > 1e-18:
        raise AssertionError('Accepted source/energy ledger failure')
    tails=[]
    for r,record in enumerate(records):
        for sign,label in ((1,'positive_peak'),(-1,'negative_hole')):
            for component,observable in ((0,'epsilon_xx'),(1,'u_z')):
                for threshold in (3.,5.):
                    mask=sign*result['z'][:,r:r+1,component] >= threshold
                    tails.append(dict(record=record,spatial=spatial,timing=timing,beta=beta,observable=observable,sign=label,threshold=threshold,
                                      **{k:v for k,v in l.occurrence(mask,DT).items() if not k.startswith('zero_event')}))
    return rows,ledger,tails,result


def run_batch(batch):
    records=list(range(batch*BATCH,min((batch+1)*BATCH,NRECORDS)))
    if not records:
        raise ValueError('Batch out of scope')
    data,portable=load_input(),Portable()
    rows,ledgers,tails=[],[],[]
    for spatial in range(3):
        for timing in range(3):
            for beta in (0.,.35):
                a,b,c,_=run_arm(data,portable,records,spatial,timing,beta)
                rows+=a; ledgers.append(b); tails+=c
                print(json.dumps(dict(batch=batch,spatial=spatial,timing=timing,beta=beta)),flush=True)
    table(HERE/f'records_{batch}.csv',rows)
    table(HERE/f'signed_{batch}.csv',tails)
    dump(HERE/f'energy_{batch}.json',ledgers)


def read_rows(pattern):
    paths = sorted(HERE.glob(pattern))
    if not paths:
        paths = [HERE/('records.csv' if pattern.startswith('records') else 'signed.csv')]
    return [row for path in paths for row in csv.DictReader(path.open())]


def summarize():
    from scipy.stats import t
    rows=read_rows('records_[0-3].csv')
    expected={(r,s,h,b) for r in range(NRECORDS) for s in range(3) for h in range(3) for b in (0.,.35)}
    keyed={(int(r['record']),int(r['spatial']),int(r['timing']),float(r['beta'])):r for r in rows}
    if len(rows)!=len(expected) or set(keyed)!=expected:
        raise AssertionError('Missing or duplicated scoped records')
    metrics=[k for k in rows[0] if k.endswith('_fraction') or k.endswith('_rate_s')]
    def stats(x):
        x=np.array(x,dtype=float); sd=float(x.std(ddof=1)); half=float(t.ppf(.975,len(x)-1)*sd/math.sqrt(len(x)))
        return dict(n_records=len(x),mean=float(x.mean()),sample_sd=sd,mean_ci95_low=float(x.mean()-half) if sd else '',
                    mean_ci95_high=float(x.mean()+half) if sd else '',ci95_halfwidth=half if sd else '',
                    interval_note='pointwise Student t' if sd else 'constant finite sample; interval unresolved, not a zero population rate')
    summaries=[]; contrasts=[]
    for s in range(3):
        for h in range(3):
            for b in (0.,.35):
                for metric in metrics:
                    x=[float(keyed[r,s,h,b][metric]) for r in range(NRECORDS)]
                    summaries.append(dict(spatial=s,timing=h,beta=b,metric=metric,**stats(x)))
                    if (s,h)!=(0,0):
                        d=[float(keyed[r,s,h,b][metric])-float(keyed[r,0,0,b][metric]) for r in range(NRECORDS)]
                        contrasts.append(dict(kind='beam_minus_uniform_steady',spatial=s,timing=h,beta=b,metric=metric,**stats(d)))
            for metric in metrics:
                d=[float(keyed[r,s,h,.35][metric])-float(keyed[r,s,h,0.][metric]) for r in range(NRECORDS)]
                contrasts.append(dict(kind='feedback035_minus_0',spatial=s,timing=h,beta='',metric=metric,**stats(d)))
                if (s,h)!=(0,0):
                    d=[float(keyed[r,s,h,.35][metric])-float(keyed[r,s,h,0.][metric])-float(keyed[r,0,0,.35][metric])+float(keyed[r,0,0,0.][metric]) for r in range(NRECORDS)]
                    contrasts.append(dict(kind='feedback_by_beam_interaction',spatial=s,timing=h,beta='',metric=metric,**stats(d)))
    table(HERE/'records.csv',[keyed[k] for k in sorted(keyed)])
    table(HERE/'summary.csv',summaries)
    table(HERE/'contrasts.csv',contrasts)
    signed = read_rows('signed_[0-3].csv')
    table(HERE/'signed.csv',signed)
    signed_summary=[]
    for s in range(3):
        for h in range(3):
            for b in (0.,.35):
                for observable in ('epsilon_xx','u_z'):
                    for sign in ('positive_peak','negative_hole'):
                        for threshold in (3.,5.):
                            selected=[r for r in signed if int(r['spatial'])==s and int(r['timing'])==h and float(r['beta'])==b and r['observable']==observable and r['sign']==sign and float(r['threshold'])==threshold]
                            for metric in ('time_fraction','initiation_rate_s'):
                                signed_summary.append(dict(spatial=s,timing=h,beta=b,observable=observable,sign=sign,threshold=threshold,metric=metric,**stats([r[metric] for r in selected])))
    table(HERE/'signed_summary.csv',signed_summary)
    paths=sorted(HERE.glob('energy_[0-3].json'))
    ledgers=[x for p in paths for x in json.loads(p.read_text())] if paths else json.loads((HERE/'energy.json').read_text())
    dump(HERE/'energy.json',ledgers)
    for r in range(NRECORDS):
        all_rows=[keyed[k] for k in expected if k[0]==r]
        if len({int(x['parents']) for x in all_rows})!=1:
            raise AssertionError('Matched physical totals diverged')
        for x in all_rows:
            if sum(int(x['mark_'+m]) for m in ('low','middle','high'))!=int(x['parents']):
                raise AssertionError('Parent/mark conservation failure')
    precision=[x for x in contrasts if x['metric']=='C_any_fraction' and x['kind']=='beam_minus_uniform_steady']
    dump(HERE/'ensemble_verification.json',dict(exact_arms=18,records_per_arm=NRECORDS,total_records=len(rows),count_conservation=True,
        max_relative_energy_residual=max(abs(x['energy_balance_residual_J'])/x['incoming_host_J'] for x in ledgers),
        practical_precision=dict(metric='paired C_any fraction',criterion_ci95_halfwidth=.005,units='fraction (0.5 percentage point)',
            achieved=sum(x['ci95_halfwidth']!='' and x['ci95_halfwidth']<=.005 for x in precision),comparisons=len(precision),
            max_halfwidth=max((x['ci95_halfwidth'] for x in precision if x['ci95_halfwidth']!=''),default='')),
        intervals='Pointwise paired/record Student t, 31 df; no multiplicity or model uncertainty coverage; no zero-rate inference from a sampled zero.'))
    print(json.dumps(dict(arms=18,records=len(rows),summary_rows=len(summaries),contrast_rows=len(contrasts))),flush=True)
    for pattern in ('records_[0-3].csv','signed_[0-3].csv','energy_[0-3].json'):
        for path in HERE.glob(pattern):
            path.unlink()


def portability(output):
    data,portable=load_input(),Portable()
    rows,ledger,_,result=run_arm(data,portable,[0],2,2,.35)
    _,_,_,control=run_arm(data,portable,[0],2,2,.35,gauge=True)
    if not np.array_equal(result['B'],control['B']) or not all(np.array_equal(result['impact'][k]['parents_by_record'],control['impact'][k]['parents_by_record']) for k in result['impact']):
        raise AssertionError('Common modal gauge changed class counts')
    import scipy
    dump(Path(output),dict(runtime=dict(numpy=np.__version__,scipy=scipy.__version__),records=rows,energy=ledger,
                          gauge_max_signed_difference=float(np.max(abs(result['signed']-control['signed']))),
                          gauge_class_counts_equal=True,field_signed_sha256=hashlib.sha256(result['signed'].astype('<f8').tobytes()).hexdigest()))


def package():
    filenames = ('replicate.py','portable.c','NUMPY_LICENSE.txt','inputs.npz','source_ledger.json',
                 'records.csv','summary.csv','contrasts.csv','signed.csv','signed_summary.csv',
                 'energy.json','ensemble_verification.json','portability.json','README.md')
    prefix = 'g13_rogue/pass2/replication/'
    payloads={prefix+name:(HERE/name).read_bytes() for name in filenames}
    manifest=dict(base=BASE,files=[dict(path=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()) for name,data in sorted(payloads.items())])
    text=(json.dumps(manifest,indent=2)+'\n').encode()
    (HERE/'manifest.json').write_bytes(text)
    payloads[prefix+'manifest.json']=text
    stream=io.BytesIO()
    with tarfile.open(fileobj=stream,mode='w',format=tarfile.PAX_FORMAT) as archive:
        for name,data in sorted(payloads.items()):
            info=tarfile.TarInfo(name)
            info.size=len(data);info.mode=0o644;info.mtime=0
            archive.addfile(info,io.BytesIO(data))
    data=gzip.compress(stream.getvalue(),mtime=0)
    (HERE/'capsule.tar.gz').write_bytes(data)
    print(json.dumps(dict(path=prefix+'capsule.tar.gz',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),members=len(payloads))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    sub.add_parser('prepare')
    run=sub.add_parser('run'); run.add_argument('--batch',type=int,required=True)
    sub.add_parser('summarize')
    sub.add_parser('package')
    check=sub.add_parser('portability');check.add_argument('--output',required=True)
    a=p.parse_args()
    if a.action=='prepare': prepare()
    elif a.action=='run': run_batch(a.batch)
    elif a.action=='summarize': summarize()
    elif a.action=='package': package()
    else: portability(a.output)


if __name__=='__main__':
    main()
