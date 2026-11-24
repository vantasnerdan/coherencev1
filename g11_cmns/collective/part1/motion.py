"""Saved G10 force-graph transients; no DFT, fitted yields or fictitious PIMD velocities.

Units are eV, angstrom, fs. The crossing observable is not a fusion rate.
Run from any directory with the installed scientific interpreter.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import time

import numpy as np
from scipy.constants import atomic_mass, electron_volt, physical_constants, hbar, c, alpha
from scipy.optimize import minimize, brentq
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from g10_quantum_site import model as M
from g4_beam.transport import bca
from g4_beam.forward import model as BH

HERE = Path(__file__).resolve().parent
MASS = M.MASS['D'] * atomic_mass * 1e10 / electron_volt  # eV fs^2/A^2
MU = MASS / 2
HBAR = hbar / electron_volt * 1e15
QE = electron_volt
AMPLITUDE = .037968  # spatial vector phasor RMS, NOT pair amplitude or time RMS
HW = .04162         # saved Hartree model clock; not the bare G10 Gamma frequency
TAU_FS = 4136.036943628514
DT = .25
STEPS = 8000
REALIZATIONS = 8
SHELL_A = 2.0


def emit(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def csv_write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader(); writer.writerows(rows)


def load_csv(name):
    with (ROOT / 'g12_residual/regression' / name).open() as stream:
        return list(csv.DictReader(stream))


def geometry(m):
    d = M.load(m['source'])
    g = d['profile']['geometry']
    return np.asarray(g['cell_A']), g


def pairs(m):
    cell, _ = geometry(m)
    i, j = np.triu_indices(len(m['K']), 1)
    d = m['positions'][i] - m['positions'][j]
    d -= np.rint(d / np.diag(cell)) * np.diag(cell)
    # All unordered occupied pairs; closest-image convention is frozen at t=0.
    return i, j, d


def minimum(m):
    shape = (len(m['K']), 3)
    result = minimize(lambda q: M.potential(q.reshape(shape), m)[0], np.zeros(np.prod(shape)),
                      jac=lambda q: -M.potential(q.reshape(shape), m)[1].ravel(), method='BFGS',
                      options={'gtol': 1e-9, 'maxiter': 1000})
    q = result.x.reshape(shape)
    if np.max(abs(M.potential(q, m)[1])) > 1e-6:
        raise ValueError('Model stationary reference unresolved')
    return q


def thermal(m, temperature, seed):
    """Independent finite classical canonical chains; exact saved potential in Metropolis.

    Momentum is the physical Maxwell marginal, NOT a ring-polymer bead velocity.
    Finite-chain position equilibration is an approximation, retained as such.
    """
    rng = np.random.default_rng(seed)
    center = minimum(m)
    q = np.tile(center, (REALIZATIONS, 1, 1))
    en = M.potential(q, m)[0]
    accepted = 0
    for _ in range(6000):
        trial = q.copy()
        sites = rng.integers(len(center), size=REALIZATIONS)
        trial[np.arange(REALIZATIONS), sites] += rng.normal(0, .10, (REALIZATIONS, 3))
        e = M.potential(trial, m)[0]
        take = np.log(rng.random(REALIZATIONS)) < -(e-en)/(M.KB*temperature)
        q[take] = trial[take]; en[take] = e[take]; accepted += int(take.sum())
    v = rng.normal(0, np.sqrt(M.KB*temperature/MASS), q.shape)
    return q, v, {'position_sampler': '6000 single-site Metropolis proposals per independent chain',
                  'acceptance': accepted/(6000*REALIZATIONS), 'temperature_K': temperature,
                  'momentum': 'physical classical Maxwell; no PIMD real-time assignment',
                  'position_equilibration_uncertainty': 'finite chains; not a quantum ensemble'}


def local_positive_energy(q, v, ref, m):
    """Positive localization diagnostic, not a unique decomposition of interacting H.

    Local excess kinetic plus local quadratic/quartic Bregman energy. Signed
    intersite interaction is excluded rather than relabelled positive energy.
    """
    x = q-ref
    xx = q*q; rr = xx.sum(-1)
    local = .5*np.einsum('...ia,iab,...ib->...i', q, m['K'], q)
    local += (m['a']-m['b'])/4*(xx*xx).sum(-1)+m['b']/4*rr*rr
    local -= (q*m['f']).sum(-1)
    er = .5*np.einsum('ia,iab,ib->i', ref, m['K'], ref)
    er += (m['a']-m['b'])/4*(ref**4).sum(-1)+m['b']/4*(ref*ref).sum(-1)**2
    er -= (ref*m['f']).sum(-1)
    # Gradient of local energy at reference, necessary for a displaced partial cell.
    grad = np.einsum('iab,ib->ia', m['K'], ref)-m['f']
    grad += ((m['a']-m['b'])*ref*ref+m['b']*(ref*ref).sum(-1,keepdims=True))*ref
    return np.maximum(.5*MASS*(v*v).sum(-1)+local-er-(grad*x).sum(-1), 0.)


def propagate(label, m, q, v, prep, *, pump=0., kicks=None, tracer=None, damping=False):
    q=q.copy(); v=v.copy(); ref=minimum(m)
    nrep, n, _=q.shape
    i,j,d=pairs(m)
    initial = float(np.mean(M.potential(q,m)[0]+.5*MASS*(v*v).sum((1,2))))
    best=np.full(nrep,np.inf); records=[None]*nrep
    peak_energy=np.zeros(nrep); peak_p=np.zeros(nrep)
    min_pr=np.full(nrep,float(n)); retain=np.zeros(nrep); supported=np.zeros(nrep)
    all_supported=np.zeros(nrep); crossing=np.zeros(nrep,int); tagged_crossing=np.zeros(nrep,int)
    first_exit=[None]*nrep; history=[]
    previous=None; longest=np.zeros(nrep); current=np.zeros(nrep)
    impulse_work=np.zeros(nrep); pump_work=np.zeros(nrep); dissipated=np.zeros(nrep)
    qeq=q.copy();veq=v.copy();feq=M.potential(qeq,m)[1]
    tagged_best=np.full(nrep,np.inf);tagged_records=[None]*nrep
    retained_initial=np.zeros(nrep);initial_hot=None;crossing_records=[]
    def external(x,t):
        force=np.zeros_like(x)
        if pump:
            force[:,:,0] -= pump*MASS*(HW/HBAR)**2*np.cos(2*HW/HBAR*t)*x[:,:,0]
        return force
    force=M.potential(q,m)[1]+external(q,0)
    for step in range(STEPS+1):
        t=step*DT
        rvec=d[None]+q[:,i]-q[:,j]
        dv=v[:,i]-v[:,j]
        r=np.linalg.norm(rvec,axis=-1)
        radial=(rvec*dv).sum(-1)/r
        energy=.5*MU*(dv*dv).sum(-1)
        p=MU*np.linalg.norm(dv,axis=-1)
        col=r.argmin(1); rr=np.arange(nrep)
        inside=np.linalg.norm(q,axis=-1)<=m['support_A']
        supported += inside.mean(1)
        all_supported += inside.all(1)
        if tracer is not None:
            tagged_columns=np.flatnonzero((i==tracer)|(j==tracer))
            tagged_col=tagged_columns[r[:,tagged_columns].argmin(1)]
            for a in range(nrep):
                k=int(tagged_col[a])
                if r[a,k]<tagged_best[a]:
                    tagged_best[a]=r[a,k]
                    tagged_records[a]=dict(tagged_closest_A=float(r[a,k]),tagged_time_fs=t,
                        tagged_partner_original=int(m['original_indices'][j[k] if i[k]==tracer else i[k]]),
                        tagged_relative_energy_eV=float(energy[a,k]),tagged_relative_momentum_eV_fs_per_A=float(p[a,k]),
                        tagged_radial_velocity_A_per_fs=float(radial[a,k]),
                        tagged_pair_within_radius_support=bool(inside[a,i[k]] and inside[a,j[k]]))
        for a in range(nrep):
            if not inside[a].all() and first_exit[a] is None:first_exit[a]=t
            k=int(col[a])
            if r[a,k]<best[a]:
                best[a]=r[a,k]
                records[a]=dict(realization=a,closest_A=float(r[a,k]),time_fs=t,
                    pair_local=[int(i[k]),int(j[k])],pair_original=[int(m['original_indices'][i[k]]),int(m['original_indices'][j[k]])],
                    tracer_pair=bool(tracer is not None and (i[k]==tracer or j[k]==tracer)),
                    relative_momentum_eV_fs_per_A=float(p[a,k]),relative_energy_eV=float(energy[a,k]),
                    radial_velocity_A_per_fs=float(radial[a,k]),radial_kinetic_eV=float(.5*MU*radial[a,k]**2),
                    pair_within_radius_support=bool(inside[a,i[k]] and inside[a,j[k]]),
                    maximum_site_displacement_A=float(np.linalg.norm(q[a],axis=-1).max()))
        peak_energy=np.maximum(peak_energy,energy.max(1));peak_p=np.maximum(peak_p,p.max(1))
        epos=local_positive_energy(q,v,ref,m)
        if kicks is not None:
            epos=np.maximum(epos-local_positive_energy(qeq,veq,ref,m),0.)
        if initial_hot is None and np.any(epos.sum(1)>1e-20):
            initial_hot=epos.argmax(1)
        s=epos.sum(1); pr=np.divide(s*s,(epos*epos).sum(1),out=np.full(nrep,float(n)),where=s>1e-25)
        min_pr=np.minimum(min_pr,pr)
        fraction=np.divide(epos.max(1),s,out=np.zeros(nrep),where=s>1e-25)
        concentrated=fraction>=.25
        current=np.where(concentrated,current+DT,0.);longest=np.maximum(longest,current)
        retain += concentrated*DT if step<STEPS else 0.
        if initial_hot is not None and step<STEPS:
            retained_initial += (np.divide(epos[np.arange(nrep),initial_hot],s,
                out=np.zeros(nrep),where=s>1e-25)>=.25)*DT
        if previous is not None:
            crossed=(previous>=SHELL_A)&(r<SHELL_A)&(radial<0)
            crossing += crossed.sum(1)
            if tracer is not None:tagged_crossing += crossed[:,(i==tracer)|(j==tracer)].sum(1)
            for a,k in np.argwhere(crossed):
                crossing_records.append(dict(realization=int(a),sample_time_fs=t,
                    pair_original=[int(m['original_indices'][i[k]]),int(m['original_indices'][j[k]])],
                    separation_A=float(r[a,k]),relative_energy_eV=float(energy[a,k]),
                    relative_momentum_eV_fs_per_A=float(p[a,k]),radial_velocity_A_per_fs=float(radial[a,k]),
                    pair_within_radius_support=bool(inside[a,i[k]] and inside[a,j[k]]),
                    crossing_semantics='sample immediately inside 2A, not a compact-contact crossing or a fusion count'))
        previous=r
        if step%40==0:
            for a in range(nrep):
                k=int(col[a])
                history.append(dict(case=label,realization=a,time_fs=t,closest_A=float(r[a,k]),
                    pair_i=int(i[k]),pair_j=int(j[k]),relative_energy_at_closest_eV=float(energy[a,k]),
                    relative_momentum_at_closest_eV_fs_per_A=float(p[a,k]),radial_velocity_at_closest_A_per_fs=float(radial[a,k]),
                    max_relative_energy_eV=float(energy[a].max()),positive_energy_PR=float(pr[a]),
                    max_local_positive_fraction=float(fraction[a]),fraction_sites_within_radius=float(inside[a].mean()),
                    total_model_energy_eV=float(M.potential(q[a],m)[0]+.5*MASS*(v[a]*v[a]).sum())))
        if step==STEPS:break
        if kicks is not None:
            # Each row (time_fs,site,energy_eV,direction) deposits energy through
            # a prescribed impulse; exact work is measured, not counted twice.
            for tk,site,ek,direction in kicks:
                if step==int(round(tk/DT)):
                    dvel=np.sqrt(2*ek/MASS)*np.asarray(direction)
                    before=.5*MASS*(v[:,site]*v[:,site]).sum(1)
                    v[:,site]+=dvel
                    impulse_work += .5*MASS*(v[:,site]*v[:,site]).sum(1)-before
        old=q.copy();v += .5*DT*force/MASS;q += .5*DT*v
        if damping:
            before=.5*MASS*(v*v).sum((1,2));v*=np.exp(-DT/TAU_FS)
            dissipated += before-.5*MASS*(v*v).sum((1,2))
        q += .5*DT*v
        ef=external(q,t+DT);force=M.potential(q,m)[1]+ef;v+=.5*DT*force/MASS
        pump_work += (external((old+q)/2,t+DT/2)*(q-old)).sum((1,2))
        if kicks is not None:
            veq+=.5*DT*feq/MASS;qeq+=.5*DT*veq
            if damping:veq*=np.exp(-DT/TAU_FS)
            qeq+=.5*DT*veq;feq=M.potential(qeq,m)[1];veq+=.5*DT*feq/MASS
        if not np.isfinite(q).all():raise FloatingPointError(label)
    final=float(np.mean(M.potential(q,m)[0]+.5*MASS*(v*v).sum((1,2))))
    for a,r in enumerate(records):
        r.update(case=label,peak_relative_energy_eV=float(peak_energy[a]),peak_relative_momentum_eV_fs_per_A=float(peak_p[a]),
            min_positive_energy_PR=float(min_pr[a]),duration_fraction_ge_quarter_fs=float(retain[a]),
            longest_localization_interval_fs=float(min(longest[a],STEPS*DT)),
            fraction_site_times_within_radius=float(supported[a]/(STEPS+1)),
            fraction_all_site_times_within_radius=float(all_supported[a]/(STEPS+1)),first_radius_exit_fs=first_exit[a],
            inward_shell_crossings=int(crossing[a]),tagged_inward_shell_crossings=int(tagged_crossing[a]))
        r['initial_hot_site_retention_ge_quarter_fs']=float(retained_initial[a])
        if tagged_records[a] is not None:r.update(tagged_records[a])
    summary=dict(case=label,model=m['name'],preparation=prep,n_sites=n,realizations=nrep,dt_fs=DT,duration_fs=STEPS*DT,
        pump_h=pump,pump_clock_meV=1000*HW,amplitude_decay_fs=TAU_FS if damping else None,
        model_support_radius_A=m['support_A'],support_meaning='radial diagnostic only; full cubic force tensor is itself a closure, not measured on every direction',
        initial_mean_model_energy_eV=initial,final_mean_model_energy_eV=final,
        mean_impulse_work_eV=float(impulse_work.mean()),mean_pump_work_eV=float(pump_work.mean()),mean_dissipated_eV=float(dissipated.mean()),
        energy_balance_error_eV=float(final-initial-impulse_work.mean()-pump_work.mean()+dissipated.mean()),
        shell_A=SHELL_A,sampled_inward_crossings=int(crossing.sum()),tagged_crossings=int(tagged_crossing.sum()),
        pair_exposure_ps=nrep*len(i)*STEPS*DT/1000,
        crossing_resolution_per_pair_ps=1000/(nrep*len(i)*STEPS*DT),
        uncertainty='deterministic cold preparations have one realization; thermal cases eight finite classical chains, not independent trajectory frames or a tail enclosure',
        crossing_records=crossing_records,
        participation_semantics='positive local energy excess relative to simultaneous same-initial-state, same-damping, no-impulse trajectory' if kicks is not None else 'positive local kinetic+onsite Bregman energy above model stationary reference; intersite signed energy excluded',
        standard_fusion_rate_per_s=None,rate_missing='no matched inner radial/contact evolution from these well trajectories')
    return summary, records, history


def pattern(m, kind, target_energy, amplitude=AMPLITUDE):
    n=len(m['K']);q=np.zeros((n,3));ref=minimum(m)
    if kind=='Gamma':q[:,0]=1.
    elif kind=='finite_q':q[:,0]=np.cos(2*np.pi*m['positions'][:,0]/8.168)
    elif kind=='localized':q[m['site'],0]=1.
    elif kind=='seeded_Gamma':
        q[:,0]=1.;q[:,1]=.1*np.cos(2*np.pi*m['positions'][:,0]/8.168)
    else:raise ValueError(kind)
    e0=M.potential(ref,m)[0]
    scale=brentq(lambda s:M.potential(ref+s*q,m)[0]-e0-target_energy,0.,4.)
    return ref+scale*q,dict(pattern=kind,scale_A=float(scale),exact_added_potential_eV=target_energy,
        q_wave_A_inverse=[2*np.pi/8.168,0,0] if kind in ('finite_q','seeded_Gamma') else None,
        phase_radians=0.,amplitude_convention='prepared peak displacement; target energy equals Gamma at declared phasor spatial RMS',
        saved_pump_amplitude_A=amplitude)


def stopping_split(target, ec):
    records=M.load('g4_beam/merit/profiles.json')['records']
    r=next(r for r in records if r['case']['condition']==target and float(r['case']['E_cm_keV'])==ec)
    layers=r['layers'];total=np.zeros(3);tracks=[]
    gx,gw=np.polynomial.hermite_e.hermegauss(5);gw/=gw.sum()
    for multiplicity,fraction in ((1,.85),(2,.15)):
        for g,w in zip(gx,gw):
            e0=2000*ec/multiplicity*(1+.05/2.3548*g)
            energy=e0; dep=np.zeros(3);timeline=[];time_fs=0.;depth=0.
            for layer in layers:
                if energy<=100.00001:break
                zt,zb,nh,nd,z,m,_=layer
                def rhs(zz,y):
                    e=max(y[0],100.)
                    split=np.array([nh*bca.se_ls(e,z)+nd*bca.se_ls(e,1),nh*bca.sn_zbl(e,z,m),nd*bca.sn_zbl(e,1,BH.M_D)])*1e-18
                    speed=np.sqrt(2*e*electron_volt/(M.MASS['D']*atomic_mass))
                    return np.r_[-split.sum(),split,1e15/speed]
                def stop(zz,y):return y[0]-100.
                stop.terminal=True;stop.direction=-1
                sol=solve_ivp(rhs,(zt,zb),np.r_[energy,dep,time_fs],events=stop,rtol=2e-7,atol=1e-8,max_step=2e-9)
                energy=float(sol.y[0,-1]);dep=sol.y[1:4,-1];time_fs=float(sol.y[4,-1]);depth=float(sol.t[-1])
                timeline.append(dict(depth_A=depth*1e10,time_fs=time_fs,energy_eV=energy,
                    electron_eV=float(dep[0]),host_recoil_eV=float(dep[1]),D_recoil_eV=float(dep[2])))
            total+=multiplicity*fraction*w*dep
            if abs(g)<1e-12:tracks.append(dict(species='D+' if multiplicity==1 else 'one D from D2+',
                initial_eV=e0,terminal_unresolved_eV=energy,terminal_time_fs=time_fs,terminal_depth_A=depth,
                components_eV=dep.tolist(),timeline=timeline))
    return dict(target=target,E_cm_keV=ec,Se_J_per_true_C=float(total[0]),host_recoil_J_per_true_C=float(total[1]),
        direct_D_recoil_J_per_true_C=float(total[2]),sum_J_per_true_C=float(total.sum()),
        residual_cutoff_J_per_true_C=115.,initial_J_per_true_C=2000*ec,
        energy_balance_J_per_true_C=float(total.sum()+115.-2000*ec),tracks=tracks,
        interpretation='mean CSDA deposition; actual recoil destinations, order and transfer efficiencies are not measured')


def capture_model():
    """One tracer at one actual deleted lattice site, NOT occupied-site double loading.

    Added well/bonds are an explicit conditional x1 closure. This is neither a
    computed trapping geometry nor a change to the fixed transport profile.
    """
    m=M.partial_cases()[-1];cell,g=geometry(m);d=M.load(m['source'])
    vacancy=int(m['vacancies'][0])
    original=np.asarray(g['original_atom_indices'])
    # original index is in the full pristine source geometry, before deletion.
    pristine=M.load(M.PD)['profile']['geometry']
    pos=np.asarray(pristine['reference_positions_A'])[vacancy]*np.diag(cell)/np.diag(np.asarray(pristine['cell_A']))
    base=M.full();m=dict(m);n=len(m['K'])
    m['K']=np.concatenate([m['K'],base['K'][:1]])
    m['f']=np.concatenate([m['f'],np.zeros((1,3))])
    m['positions']=np.vstack([m['positions'],pos]);m['original_indices']=np.r_[m['original_indices'],vacancy]
    # Derive source-shell radii and means, preserving original occupied J.
    source=M.load(M.PD)['coupling']['shells']
    delta=m['positions'][:-1]-pos;delta-=np.rint(delta/np.diag(cell))*np.diag(cell)
    distances=np.linalg.norm(delta,axis=1);new=np.zeros(n)
    for shell in source:
        radius=shell['distance_A']*np.diag(cell)[0]/8.168
        new[np.isclose(distances,radius,atol=2e-5)]=shell['sum_pulling_J_xx_eV_A2']/shell['distinct_supercell_H_sites']
    m['J']=np.block([[m['J'],new[:,None]],[new[None,:],np.zeros((1,1))]])
    m['name'] += '_one_captured_tag_conditional'
    return m,n,dict(capture='one tagged stopped D at actual deleted original lattice site',vacancy_original_index=vacancy,
        tracer_local_index=n,added_density=0,occupancy='one selected conditional capture; no second D at an occupied site',
        conditional_closure='new-site K and bonds transferred x1, F0=0; no relaxed capture potential supplied',
        cutoff_eV_is_vibrational_energy=False)


def rate_tables():
    # g6_pairlaw matching: allowed-spin normalized s-wave mixture, proton only.
    mu_kg=M.MASS['D']*atomic_mass/2
    s0=BH.BH['p'][0]/1000*1000*electron_volt*1e-28  # keV barn -> J m^2
    a_beam=s0/(np.pi*alpha*mu_kg*c)*1e30
    A=.75*a_beam
    data=load_csv('davis_residuals.csv');si={r['row_id']:r for r in load_csv('SI1_normalizations.csv')}
    add={r['row_id']:r for r in load_csv('additive_rows.csv')}
    rows=[]
    for d in data:
        if d['channel']!='p' or float(d['E_cm_keV'])>2.5:continue
        e=float(d['E_cm_keV']);cf=float(d['cup_focusing_fraction']);current=float(BH.I_cup(e))*1e-6/cf
        raw=float(d['model_corrected_raw_reactions_per_true_C']);r0=float(d['R_model'])
        # Infer the model's own raw neutron anchor, never the measured anchor.
        denom=raw/cf/r0; ref_current=float(BH.I_cup(6.5))
        for power in (0,1):
            readout_factor=(ref_current/float(BH.I_cup(e)))**power
            # Xi_event = integral excess summed contact dt, A^-3 s per stopped D event.
            # C_domain = stationary summed normalized-pair contact, A^-3 per active domain.
            event_per_c=A*1.15/QE
            stationary_per_c=A/current
            a=add[d['row_id']];s=si.get(d['row_id'],{})
            rows.append(dict(target=d['target'],E_cm_keV=e,reading='k'+str(power),
                predicted_emitted_p_per_true_C=None,predicted_emitted_p_per_s=None,
                rate_status='unknown contact; computed shell flux is not a nuclear event count',
                conditional_event_emitted_per_C_per_Xi=event_per_c,
                conditional_event_emitted_per_s_per_Xi=event_per_c*current,
                conditional_stationary_emitted_per_s_per_Cdomain=A,
                conditional_stationary_emitted_per_C_per_Cdomain=stationary_per_c,
                conditional_event_deltaR_per_Xi=event_per_c/cf/denom*readout_factor,
                conditional_stationary_deltaR_per_Cdomain=stationary_per_c/cf/denom*readout_factor,
                modeled_own_n65_readout_per_C=denom,own_reference_is_measured_absolute=False,
                R_observed=float(d['R_observed']),R_model_k=float(d['R_model'])*readout_factor,
                adjacent_G12_Ak=float(d['R_observed'])-float(d['R_model'])*readout_factor,
                SI1_detected_per_s=s.get('detected_per_s_proxy',''),
                SI1_detected_per_C=s.get('detected_per_true_C_model_collection_proxy',''),
                SI1_acceptance='unknown; no subtraction or emitted-rate calibration',
                matched_potential='PdD full or partial conditional local model only' if d['target'].startswith('Pd') else 'no matched Ti(D) anharmonic force graph'))
    csv_write(HERE/'rates.csv',rows)
    return dict(Ap_allowed_spin_A3_per_s=A,Ap_random_nine_spin_A3_per_s=.5*a_beam,
        convention='g6_pairlaw: Gamma_p=A_p rho_contact per normalized s-wave pair; no empirical resonance or attempt frequency',
        event='Xi=integral [sum rho_event(t)-sum rho_matched_eq(t)] dt, A^-3 s; 1.15 incoming D/e, conditional capture probability multiplies once',
        stationary='Cdomain=sum rho_contact over unordered normalized pairs in one driven domain; active-domain count is an explicit multiplier, not inferred from N',
        unknowns=['inner radial Hamiltonian and contact evolution','capture fraction and microscopic stopping assignment','active driven-domain population and lifetime'],
        shell='inward crossings at 2 A, actual joint r,p sampled; no supplied matching transmission from shell to nuclear contact')


def main():
    start=time.perf_counter();HERE.mkdir(parents=True,exist_ok=True)
    full=M.full();partial=M.partial_cases()[-1]
    summaries=[];realizations=[];histories=[]
    refining='--refine' in sys.argv
    retained_producer=None
    if refining:
        prior=json.loads((HERE/'results.json').read_text())
        retained_producer=prior['script_sha256']
        summaries=prior['cases']
        with (HERE/'realizations.csv').open() as f:realizations=list(csv.DictReader(f))
        with (HERE/'history.csv').open() as f:histories=list(csv.DictReader(f))
    def run(*args,**kw):
        label=args[0]
        if refining and not (label.startswith('track_eta') or label.startswith('captured_') or label in ('thermal_300K','partial_300K')):return
        s,r,h=propagate(*args,**kw)
        summaries[:]=[row for row in summaries if row['case']!=label]
        realizations[:]=[row for row in realizations if row['case']!=label]
        histories[:]=[row for row in histories if row['case']!=label]
        summaries.append(s);realizations.extend(r);histories.extend(h)
    zeros=np.zeros((1,len(full['K']),3))
    for temp in (77.,300.):
        q,v,prep=thermal(full,temp,110+int(temp));run('thermal_'+str(int(temp))+'K',full,q,v,prep)
        if temp==300.:thermal_q,thermal_v=q,v
    partial_q,partial_v,prep=thermal(partial,300.,412)
    run('partial_300K',partial,partial_q,partial_v,prep)
    for factor in (1.,2.):
        reference=np.zeros((len(full['K']),3));reference[:,0]=factor*AMPLITUDE
        target=float(M.potential(reference,full)[0])
        for kind in ('Gamma','finite_q','localized','seeded_Gamma'):
            qq,prep=pattern(full,kind,target,factor*AMPLITUDE)
            run('cold_'+kind+'_'+str(factor),full,qq[None],zeros,prep,
                pump=.08 if kind=='seeded_Gamma' else 0.,damping=kind=='seeded_Gamma')
    run('track_matched_300K',full,thermal_q,thermal_v,
        dict(kind='same initial thermal state and damping as event runs, no impulses'),damping=True)
    target=float(M.potential(np.tile([AMPLITUDE,0.,0.],(len(partial['K']),1))+minimum(partial),partial)[0]
                 -M.potential(minimum(partial),partial)[0])
    qq,prep=pattern(partial,'Gamma',target)
    run('partial_coherent_seed',partial,qq[None],np.zeros_like(qq)[None],
        dict(prep,kind='coherent displacement on actual heterogeneous partial K F0 J'),pump=.08,damping=True)
    # Same thermal positions/momenta and damping in all event/control cases.
    splits=json.loads((HERE/'stopping.json').read_text()) if refining else [stopping_split(t,e) for t in ('Pd','PdD','Ti','TiD') for e in (.25,2.5)]
    budget=next(s for s in splits if s['target']=='PdD' and s['E_cm_keV']==.25)
    direct=budget['direct_D_recoil_J_per_true_C']/1.15
    for efficiency in (0.,.001):
        transferred=(budget['Se_J_per_true_C']+budget['host_recoil_J_per_true_C'])/1.15*efficiency
        # Only this declared share goes into the small graph. The rest remains
        # electron/host/other-D deposition, never a hidden second injection.
        local_energy=.001*direct+transferred
        for shape in ('impulsive','delayed'):
            times=[0.] if shape=='impulsive' else [0.,25.,50.,75.]
            direction=np.array([1.,0.,0.]);kicks=[(t,full['site'],local_energy/len(times),direction) for t in times]
            run('track_eta'+str(efficiency)+'_'+shape,full,thermal_q,thermal_v,
                dict(kind='one-ion conditional small-cell deposition',total_ledger_D_recoil_per_D_eV=direct,
                    modeled_direct_share=.001,electron_host_to_D_efficiency=efficiency,
                    nominal_kick_energy_eV=local_energy,shape=shape,
                    uncertainty='nearest site, order and split not measured; impulse work includes thermal cross term',
                    unassigned_deposition='remains electron/host/other D; terminal100eV unresolved, never injected as vibration'),
                kicks=kicks,damping=True)
    captured,tag,prep=capture_model();cq,cv,therm=thermal(captured,300.,713)
    run('captured_tag_equilibrium_300K',captured,cq,cv,dict(prep,thermal_reference=therm),tracer=tag)
    qcapture=np.concatenate([partial_q,np.zeros((REALIZATIONS,1,3))],axis=1)
    vtag=np.random.default_rng(714).normal(0,np.sqrt(M.KB*300/MASS),(REALIZATIONS,1,3))
    vcapture=np.concatenate([partial_v,vtag],axis=1)
    run('captured_tag_postcapture_300K',captured,qcapture,vcapture,
        dict(prep,stage='others from partial equilibrium; one tracer at vacancy reference with Maxwell momentum; conditional post-stopping initialization, no100eV residual injection'),tracer=tag)
    csv_write(HERE/'realizations.csv',realizations);csv_write(HERE/'history.csv',histories)
    emit(HERE/'stopping.json',splits);rates=rate_tables()
    source_paths=['g10_quantum_site/model.py',M.PD,partial['source'],
        'g3_lattice/rerun/part2/updated/part2.md','g3_lattice/rerun/part2/updated/parameters.json',
        'g3_lattice/rerun/part2/updated/pd_summary.json','g12_residual/inputs_ledger.md',
        'g12_residual/regression/davis_residuals.csv','g12_residual/regression/additive_rows.csv',
        'g12_residual/regression/SI1_normalizations.csv','g4_beam/merit/profiles.json','g6_pairlaw/README.md']
    saved=M.load('g3_lattice/rerun/part2/updated/pd_summary.json')['result']['results']
    carrier=next(r for r in saved if r['label']=='random_x1_h0.08')['half_pump_correlations']['amplitude_weighted']
    wave=np.array(carrier['peak_q']['q_cartesian_A_inverse'])
    _,_,bonds=pairs(full);nearest=bonds[np.isclose(np.linalg.norm(bonds,axis=1),np.linalg.norm(bonds,axis=1).min())]
    carrier_pair_amplitude=2*AMPLITUDE*np.abs(np.sin(nearest@wave/2))
    emit(HERE/'results.json',dict(cases=summaries,rate_contract=rates,
        mass_eV_fs2_per_A2=MASS,mu_eV_fs2_per_A2=MU,
        saved_clock_meV=HW*1000,bare_Gamma_meV=1000*HBAR*np.sqrt((full['K'][0,0,0]-full['J'][0].sum())/MASS),
        saved_amplitude_spatial_vector_RMS_A=AMPLITUDE,saved_spatial_phase='full per-site complex field not retained in delivered pump packet; no Gamma inference from scalar RMS',
        constructed_patterns='Gamma, q=(2pi/8.168,0,0) and localized are declared probes, not reconstructions of the saved pump',
        saved_carrier=carrier['peak_q'],saved_Gamma_coherent_fraction=carrier['Gamma_coherent_fraction'],
        saved_carrier_weight=carrier['peak_coherent_fraction'],
        ideal_single_carrier_NN_relative_peak_A=float(carrier_pair_amplitude.max()),
        ideal_single_carrier_NN_relative_peak_energy_eV=float(.5*MU*(HW/HBAR*carrier_pair_amplitude.max())**2),
        carrier_meaning='analytic single-carrier illustration at actual q and scalar amplitude, unknown global phase; full field is not reconstructed and the 8.168A force graph cannot represent a 7.545nm Bloch carrier',
        retained_unaffected_producer_sha256=retained_producer,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),runtime_seconds=time.perf_counter()-start,
        versions=dict(python=sys.version.split()[0],numpy=np.__version__),
        omissions=['No inner pair potential/contact dynamics; shell flux alone does not identify emitted fusion rates.',
                   'Old quartic extrapolation and scalar-shell closures, especially partial .03 A support; wide/connected cuts not supplied.',
                   'Microscopic stopping-to-D transfer, conditional capture preparation and domain population/loss are not measured; no Ti(D) graph.']))
    # One settled changed-consumer proof; this targets COM cancellation and work accounting.
    gamma=[r for r in realizations if r['case']=='cold_Gamma_1.0']
    eq=[s for s in summaries if s['case']=='cold_Gamma_1.0'][0]
    assert max(float(r['peak_relative_energy_eV']) for r in gamma)<1e-24
    assert abs(float(eq['energy_balance_error_eV']))<1e-6
    assert max(abs(s['energy_balance_J_per_true_C']) for s in splits)<.01
    assert rates['Ap_allowed_spin_A3_per_s']>0
    emit(HERE/'consumer_proof.json',dict(uniform_pair_relative_energy_max_eV=max(float(r['peak_relative_energy_eV']) for r in gamma),
        uniform_free_energy_error_eV=float(eq['energy_balance_error_eV']),
        max_stopping_balance_J_per_true_C=max(abs(s['energy_balance_J_per_true_C']) for s in splits),
        no_fusion_rate_from_kinetic_extrema=True,physical_PIMD_velocities_used=False,
        event_excess_uses_simultaneous_matched_reference=True,tagged_closest_approach_separately_retained=True,
        all_realizations_retained=len(realizations),histories_retained=len(histories)))
    print(json.dumps({'cases':len(summaries),'realizations':len(realizations),'runtime_s':time.perf_counter()-start}))


if __name__=='__main__':
    main()
