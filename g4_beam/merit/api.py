"""Shared pass5 physical merit and paper-axis comparison; SI profile input.

All depth integration is delegated to forward.consumer. No target parameter is
inferred from a material name except the explicitly labelled pass4 fit family.
"""
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
import csv
import json
import numpy as np
from scipy.optimize import minimize_scalar
from g4_beam.forward import model as M
from g4_beam.forward.consumer import forward_yields, stopping_eV_per_m, QE
from g4_beam.forward.profile_thermal import transport_stopping

HERE = Path(__file__).resolve().parent
BEAM = HERE.parent
CAUSES = ('screeningrise', 'surface', 'neutrals')
ENERGIES = np.array([.25,.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5])
FOUR = np.array([.25,1.,2.5,6.5])
MEAS = {'M0':(0,0,0), 'M1':(1,0,0), 'M2':(0,1,0),
        'M3':(1,1,0), 'M4':(0,1,1), 'M5':(1,1,1)}

@dataclass
class Target:
    name: str
    profiles: object
    material: str = 'PdD'
    stopping_fn: object = stopping_eV_per_m
    columns: tuple = ()
    measurement: dict = field(default_factory=dict)
    Ue_shift_keV: object = 0.
    _cache: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        if self.material not in M.MATS:
            raise ValueError('material selects one of the four pass4 fitted families')
        if self.columns and (any(w < 0 for w,_ in self.columns) or
                             abs(sum(w for w,_ in self.columns)-1)>1e-12):
            raise ValueError('Parallel-column area fractions must sum to one')

@lru_cache(None)
def bestfits():
    rows=json.loads((HERE/'fits_transported.json').read_text())
    return {c:min((r for r in rows if r['cause']==c),key=lambda r:r['chi2']) for c in CAUSES}

@lru_cache(None)
def uniform_reference(material='PdD'):
    h=M.ELEM[M.MATERIALS[material]['host']]
    nh=M.host_density(material)*1e6
    layer=(0.,2e-6,nh,nh*M.MATERIALS[material]['x'],h['Z'],h['M'],nh)
    return Target('uniform_closure_'+material,lambda e:([layer],0.,{'source':'pass4 uniform target'}),material)


@lru_cache(None)
def profile_records():
    return json.loads((HERE/'profiles.json').read_text())['records']


@lru_cache(None)
def reference(material='PdD'):
    records={float(r['case']['E_cm_keV']):r for r in profile_records() if r['case']['condition']==material}
    def profile(e):
        r=records[float(e)]
        return r['layers'],r['nads_m2'],dict(r['state'],**r['case'])
    return Target('UC_Davis_transported_'+material,profile,material,transport_stopping)


def _settings(target,cause,bare,channel,fit=None):
    fit=bestfits()[cause] if fit is None else fit
    p=fit['params']; k=fit['current_power']; bh,bl=fit['b_hi'],fit['b_lo']
    fraction=lambda e:M.collected_fraction(e,bh,bl)
    measurement=dict(current_A=lambda e:M.I_cup(e)*1e-6/fraction(e),
                     collected_fraction=fraction,neutral_reference_current_A=M.I_INT*1e-6,
                     extra_current_power=k,channel=channel)
    measurement.update(target.measurement)
    measurement['channel']=channel
    params={'additional_Ue_keV':target.Ue_shift_keV}
    cand={'screeningrise':'sigmoid_tt','surface':'layer','neutrals':'fast_neutrals'}[cause]
    if cause=='screeningrise':params.update(Ue_s=p['Ues_'+target.material],E0=p['E0'])
    if cause=='surface':params.update(Ue_s=p['Ues_'+target.material],d_nm=p['d'])
    if cause=='neutrals':params.update(f_n=p['f_n'],E_n_lab_keV=p['E_n'],delta=int(fit['cand'][-1]))
    if bare:
        # Same beam, including neutrals, and identical target/stopping; remove
        # screening AND common nuclear background for the bare-nuclei denominator.
        cand='fast_neutrals' if cause=='neutrals' else 'none'
        params['additional_Ue_keV']=0.
    return 0. if bare else p['Ue_'+target.material],cand,params,measurement


def evaluate(target, energies, cause, *, bare=False, channel='p'):
    """Cached actual consumer result. Target and profile snapshots are immutable in use."""
    ee=np.atleast_1d(np.asarray(energies,float))
    key=(tuple(ee),cause,bare,channel)
    if key in target._cache:return target._cache[key]
    out=evaluate_fit(target,ee,cause,bestfits()[cause],bare=bare,channel=channel)
    target._cache[key]=out
    return out


def evaluate_fit(target,energies,cause,fit,*,bare=False,channel='p'):
    """Raw shared-consumer yields at one explicit transported-fit row.

    This does NOT divide by an implicit reference. Profile-likelihood consumers
    must evaluate BOTH changed target and reference at the same supplied fit.
    """
    ee=np.atleast_1d(np.asarray(energies,float))
    u,cand,p,measurement=_settings(target,cause,bare,channel,fit)
    columns=target.columns or ((1.,target.profiles),)
    outputs=[]
    for weight,profiles in columns:
        if weight==0:continue
        out=forward_yields(ee,profiles,Ue_keV=u,candidate=cand,parameters=p,
                           measurement=measurement,stopping_fn=target.stopping_fn,
                           reaction='eq7' if bare else 'proton')
        outputs.append((weight,out))
    keys=('physical_yield_per_C','counts_per_C','reported_signal','ion_yield_per_C','neutral_yield_per_C')
    out={name:sum(w*r[name] for w,r in outputs) for name in keys}
    out['reference_signal']=sum(w*r['reference_signal'] for w,r in outputs)
    out['R']=out['reported_signal']/out['reference_signal']
    return out


def merit(target,E,*,cause=None):
    """Dimensionless physical proton yield relative to fixed PdD at SAME CM keV."""
    if cause is None:return {c:merit(target,E,cause=c) for c in CAUSES}
    if cause not in CAUSES:raise ValueError('Unknown plateau cause')
    n=evaluate(target,E,cause)['physical_yield_per_C']
    d=evaluate(reference(),E,cause)['physical_yield_per_C']
    if np.any(d<=0):raise ValueError('Positive physical reference required')
    ratio=n/d
    return float(ratio[0]) if np.ndim(E)==0 else ratio


def eq7_fit(energies,R,weights=None):
    """Literal Eq7 ratio fitted over paper's2.5–6.5keV interval; log RMS."""
    ee=np.asarray(energies); rr=np.asarray(R)
    ok=(ee>=2.5)&(ee<6.5)&(rr>0)&np.isfinite(rr)
    if weights is not None:ok &= np.asarray(weights)>0
    if not np.any(ok):return None,None
    wt=np.ones(np.sum(ok)) if weights is None else np.asarray(weights)[ok]
    def residual(u):return np.log(M.sigma(ee[ok],u)/M.sigma(6.5,u))-np.log(rr[ok])
    result=minimize_scalar(lambda u:np.sum(wt*residual(u)**2),bounds=(0,200),method='bounded',options={'xatol':1e-9})
    return float(result.x),float(np.sqrt(np.mean(residual(result.x)**2)))

@lru_cache(None)
def paper_rows():
    with (BEAM/'forward/errbars.csv').open() as f:
        return list(csv.DictReader(line for line in f if not line.startswith('#')))

@lru_cache(None)
def paper_fit(material,power):
    rows=[r for r in paper_rows() if r['material']==material and float(r['E_cm'])<6.5]
    ee=np.array([float(r['E_cm']) for r in rows]); rr=np.array([float(r['R']) for r in rows])
    rr*= (M.I_cup(6.5)/M.I_cup(ee))**power
    weights=1/np.array([float(r['sigma_lnR']) for r in rows])**2
    return eq7_fit(ee,rr,weights)[0]


def paper_weights(material,energies):
    table={float(r['E_cm']):float(r['sigma_lnR']) for r in paper_rows() if r['material']==material and float(r['E_cm'])<6.5}
    return np.array([1/table[float(e)]**2 if float(e) in table else 0. for e in energies])


def _reported(target,energies,cause,fit=None):
    # Same fitted proton/neutron stitch for baseline and explicit-fit consumers.
    def run(ee,channel='p'):
        return evaluate(target,ee,cause,channel=channel) if fit is None else evaluate_fit(target,ee,cause,fit,channel=channel)
    pp=run(np.r_[energies,4.5,6.5]);nn=run([4.5,6.5],channel='n')
    stitch=nn['reported_signal'][0]/pp['reported_signal'][-2]
    s=(bestfits()[cause] if fit is None else fit)['params']['s']
    neutron=run(energies,channel='n')
    return np.where(energies<=4.,pp['reported_signal'][:-2]*stitch*s,
                    neutron['reported_signal'])/nn['reported_signal'][-1]


def compare(target,*,energies=None):
    """Row-per-cause/energy public table; fixed PdD overlay retains density gain."""
    ee=ENERGIES if energies is None else np.asarray(energies,float)
    rows=[]
    for cause in CAUSES:
        grid=np.unique(np.r_[ee,6.5])
        out=evaluate(target,grid,cause); idx=np.searchsorted(grid,ee); iref=np.searchsorted(grid,6.5)
        y=out['physical_yield_per_C'][idx]; yref=out['physical_yield_per_C'][iref]
        base=evaluate(target,grid,cause,bare=True)['physical_yield_per_C'][idx]
        references={m:evaluate(reference(m),grid,cause)['physical_yield_per_C'] for m in M.MATS}
        bare_refs={m:evaluate(reference(m),grid,cause,bare=True)['physical_yield_per_C'] for m in ('Pd','PdD')}
        reported=_reported(target,ee,cause)
        # The cause at its fitted measurement reading reconstructs the plotted
        # curve. Apply the SAME extra inverse-current powers to this curve and
        # the actual digitized paper curves, without undoing its native fit.
        as_plotted=reported
        ue={j:eq7_fit(ee,as_plotted*(M.I_cup(6.5)/M.I_cup(ee))**j,paper_weights('PdD',ee)) for j in range(3)}
        ue_pd={j:eq7_fit(ee,as_plotted*(M.I_cup(6.5)/M.I_cup(ee))**j,paper_weights('Pd',ee)) for j in range(3)}
        for i,e in enumerate(ee):
            nD=(target.measurement.get('D1_fraction',.85)+2*target.measurement.get('D2_fraction',.15))/QE
            neutral_atoms=0.
            if cause=='neutrals':
                _,_,p,measurement=_settings(target,cause,False,'p')
                current=measurement['current_A'];current=float(current(e) if callable(current) else current)
                neutral_atoms=p['f_n']*measurement['neutral_reference_current_A']/current/QE
            all_D=nD+neutral_atoms
            enhancement=y[i]/base[i] if base[i]>0 else None
            row=dict(target=target.name,cause=cause,E_cm_keV=float(e),physical_reactions_per_true_C=float(y[i]),
                     charged_D_per_true_C=nD,reactions_per_incident_charged_D=float(y[i]/nD),
                     incident_D_including_neutrals_per_true_C=all_D,
                     reactions_per_incident_all_D=float(y[i]/all_D),
                     ion_reactions_per_true_C=float(out['ion_yield_per_C'][idx[i]]),
                     neutral_reactions_per_true_C=float(out['neutral_yield_per_C'][idx[i]]),
                     merit_PdD=float(y[i]/references['PdD'][idx[i]]),ratio_to_Pd=float(y[i]/references['Pd'][idx[i]]),
                     enhancement_over_bare_same_target=enhancement,
                     enhancement_ratio_to_PdD=None if enhancement is None else enhancement/(references['PdD'][idx[i]]/bare_refs['PdD'][idx[i]]),
                     enhancement_ratio_to_Pd=None if enhancement is None else enhancement/(references['Pd'][idx[i]]/bare_refs['Pd'][idx[i]]),
                     fixed_PdD_reference_curve=float(y[i]/references['PdD'][iref]),self_R=float(y[i]/yref),reported_R=float(reported[i]),
                     reference_status='actual transported run-end PdD + adsorbate; refitted causes; same-energy physical proton reference')
            row.update(cause_fit_chi2=bestfits()[cause]['chi2'],cause_fit_dof=bestfits()[cause]['dof'],
                       cause_fit_current_power=bestfits()[cause]['current_power'],
                       effective_fit_scope='literal Eq7;2.5<=E<6.5; PdD paper energy/error weights; matched reference-specific ratios')
            for j in range(3):
                row.update({f'self_R_current_power_{j}':float(as_plotted[i]*(M.I_cup(6.5)/M.I_cup(e))**j),
                            f'Ue_effective_{j}_keV':ue[j][0],f'Ue_fit_log_rms_{j}':ue[j][1],
                            f'Ue_ratio_PdD_{j}':ue[j][0]/paper_fit('PdD',j),f'Ue_ratio_Pd_{j}':ue_pd[j][0]/paper_fit('Pd',j)})
            for m in M.MATS:
                obs=next((float(r['R']) for r in paper_rows() if r['material']==m and abs(float(r['E_cm'])-e)<1e-8),None)
                si=np.flatnonzero(np.isclose(M.SI1_E,e))
                counts=float(M.counts_per_C_table(m)[si[0]]) if len(si) else None
                # Actual SI1 SAME-energy counts ratios, mapped onto the fixed
                # modeled PdD6.5 scale with an explicit model bridge.
                same_ratio=counts/float(M.counts_per_C_table('PdD')[si[0]]) if len(si) else None
                anchor=same_ratio*references['PdD'][idx[i]]/references['PdD'][iref] if same_ratio is not None else None
                row.update({f'model_{m}_merit_PdD':float(references[m][idx[i]]/references['PdD'][idx[i]]),
                            f'experiment_{m}_self_R':obs,f'experiment_{m}_counts_per_cup_C':counts,
                            f'experiment_{m}_fixed_PdD_anchor':anchor,
                            f'experiment_{m}_same_E_ratio_PdD':same_ratio,
                            f'experiment_{m}_anchor_method':'SI1 same-E counts ratio times modeled PdD(E)/PdD6.5' if anchor is not None else None,
                            f'paper_{m}_quoted_Ue_keV':{'Pd':.8,'PdD':1.7,'Ti':.3,'TiD':.9}[m]})
                for j in range(3):row[f'matched_{m}_Ue_effective_{j}_keV']=paper_fit(m,j)
            rows.append(row)
    return rows


def write_csv(path,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
