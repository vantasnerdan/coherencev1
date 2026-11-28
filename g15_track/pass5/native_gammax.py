"""Measured longitudinal Gamma/X increments; named shell-1 continuation.

The incoming source, pressure, G3 room base and heat remain accepted data. This
is a propagation/meeting allocation, not a regenerated hot-track source or a
unique material dispersion. Endpoint operators and their saved runs are kept.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from common import HERE,ENERGIES,BEAMS,KB,inherited,digest,emit,table
import meetings

NAME='native_Gamma_X'
FROZEN_COUPLING='coupling_frozen_launch'
PROPAGATION_ALLOCATIONS=(NAME,FROZEN_COUPLING)
INPUT_DIR='native_Gamma_X_inputs'
EXTRA_SOURCES=('native_gammax.py',*(INPUT_DIR+'/'+n for n in ('folded-response.json','thermal-changes.json','crossover.json','provenance.json')))


class Law:
    allocations=('onsite','coupling',NAME,FROZEN_COUPLING)
    def __init__(self,path4,mesh=16,low_A_power=1.):
        self.mesh=mesh;self.low_A_power=low_A_power
        receipt=json.loads((HERE/INPUT_DIR/'provenance.json').read_text())
        for row in receipt['files']:
            if digest(HERE/row['local'])!=row['sha256']:raise ValueError('Changed selected native Gamma/X input')
        folded=json.loads((HERE/INPUT_DIR/'thermal-changes.json').read_text())
        crossover=json.loads((HERE/INPUT_DIR/'crossover.json').read_text())
        self.folded=[r[f'mesh{mesh}'] for r in folded['rows']]
        low=[r[f'mesh{mesh}']['contrast'] for r in crossover['rows'] if r['temperature_K']<3000]
        knots=[(0.,0.)]+[(r['electronic_excitation_eV_per_D_site'],r['K_change_eV_A2']) for r in low]
        knots += [(r['electronic_excitation_eV_per_D_site'],r['mode_K_change_eV_A2']['gamma_100']) for r in self.folded[1:]]
        self.u,self.gamma=np.asarray(sorted(knots)).T
        self.au=np.array([r['electronic_excitation_eV_per_D_site'] for r in self.folded])
        self.A=np.array([.75*(r['mode_K_change_eV_A2']['x_100']-r['mode_K_change_eV_A2']['gamma_100']) for r in self.folded])
        s,native,*_=inherited(path4);self.lattice=s.inputs()['lattice']
        self.shell=np.asarray(self.lattice['shells'][0]['periodic_alias_offsets_A'])*1e-10
        old=native.Material(16);self.old_u=old.u;self.old_delta=old.dk

    @staticmethod
    def interpolate(u,knots,values):
        u=np.asarray(u,float)
        if np.any(~np.isfinite(u)) or np.any(u< -1e-13) or np.any(u>knots[-1]*(1+1e-12)):
            raise ValueError('Native projection leaves measured excitation support; no extrapolation')
        index=np.clip(np.searchsorted(knots,u,side='right')-1,0,len(knots)-2)
        return np.interp(u,knots,values),np.diff(values)[index]/np.diff(knots)[index]

    def amplitudes(self,u):
        gamma,slope=self.interpolate(u,self.u,self.gamma)
        A,dA=self.interpolate(u,self.au,self.A)
        if self.low_A_power!=1.:
            u=np.asarray(u);low=u<self.au[1];x=np.maximum(u,0.)/self.au[1]
            A=np.where(low,self.A[1]*x**self.low_A_power,A)
            dA=np.where(low,self.A[1]*self.low_A_power*x**(self.low_A_power-1)/self.au[1],dA)
        return gamma,slope,A,dA

    def project(self,q,u,allocation=NAME):
        q=np.asarray(q);u=np.asarray(u);argument=q@self.shell.T
        F=np.cos(argument).mean(axis=-1)
        gradF=-np.einsum('...b,bi->...i',np.sin(argument),self.shell)/len(self.shell)
        if allocation in ('onsite','coupling',FROZEN_COUPLING):
            delta,slope=self.interpolate(u,self.old_u,self.old_delta)
            if allocation=='onsite':return delta,slope,np.zeros(np.broadcast_shapes(q.shape[:-1],u.shape)+(3,))
            return delta*F,slope*F,delta[...,None]*gradF
        if allocation!=NAME:raise ValueError('Unknown native projection allocation')
        gamma,slope,A,dA=self.amplitudes(u)
        return gamma+A*(1-F),slope+dA*(1-F),-A[...,None]*gradF


class FrozenIncomingMaterial(meetings.Material):
    """Same archived onsite incoming packet, changed common propagation law."""
    def __init__(self,accepted,lattice):
        super().__init__(accepted,lattice);self.law=Law(accepted.root)
    def electronic_project(self,u,mode,allocation):
        u=np.asarray(u);q=self.k[mode].reshape((len(mode),)+(1,)*(u.ndim-1)+(3,))
        return self.law.project(q,u,allocation)
    def source_project(self,delta,mode,allocation):
        return delta,0.


class CachedStrain:
    """Reuse exact named source-chart queries; no surrogate interpolation."""
    def __init__(self,sampler,out,prefix):
        self.sampler=sampler;self.path=Path(out)/('.'+prefix+'source_strain_cache.npz')
        names=('live.py','receiver.py','low.py','common.py','low_modes.npz','selected_inputs.json')
        self.epoch=hashlib.sha256(''.join(digest(Path(out)/n) for n in names).encode()).hexdigest()
        self.cache={}
        if self.path.exists():
            with np.load(self.path) as data:
                if str(data['_epoch'])==self.epoch:self.cache={k:data[k] for k in data.files if k!='_epoch'}
    def sample(self,ec,beam,xyz_m,mount='nominal',time_s=2,allocation='onsite',control='joined',closure='nominal_retained'):
        source_allocation='onsite' if allocation in PROPAGATION_ALLOCATIONS else allocation
        xyz=np.asarray(xyz_m,float);times=np.asarray(time_s,float)
        key=hashlib.sha256(json.dumps((ec,beam,mount,source_allocation,control,closure,xyz.shape,times.shape)).encode()+xyz.tobytes()+times.tobytes()).hexdigest()
        if key not in self.cache:
            self.cache[key]=self.sampler.sample(ec,beam,xyz,mount,time_s,source_allocation,control,closure)
            np.savez_compressed(self.path,_epoch=np.array(self.epoch),**self.cache)
        return self.cache[key]


def controls(out,path4):
    from medium import Medium
    s,*_=inherited(path4);a=s.inputs()['lattice']['a_A']*1e-10
    rows=[];error=0.
    for mesh in (12,16):
        law=Law(path4,mesh);q=np.array([[0.,0.,0.],[2*np.pi/a,0.,0.]])
        F=np.cos(q@law.shell.T).mean(axis=1)
        if abs(F[1]+1/3)>1e-12:raise ValueError('G3 first-shell projector does not match the measured X convention')
        for r in law.folded:
            u=r['electronic_excitation_eV_per_D_site'];g=r['mode_K_change_eV_A2']['gamma_100'];x=r['mode_K_change_eV_A2']['x_100']
            actual=law.project(q,np.array([u,u]))[0];error=max(error,float(np.max(abs(actual-[g,x]))))
            rows.append(dict(mesh=mesh,U_eV_per_formula=u,Delta_K_Gamma=g,Delta_K_X=x,Delta_diagonal=(g+x)/2,Delta_opposing_J_at_Gamma=(x-g)/2,
                shell1_A=.75*(x-g),old_onsite_X=g,old_coupling_X=-g/3,
                scope='measured longitudinal Gamma/X increments; shell-1 A is a named nonunique continuation, not a convex endpoint mixture or material vg'))
    if error>1e-12:raise ValueError('Native projection loses an exact Gamma/X datum')
    law=Law(path4);q=np.array([[.19,.31,.23]])*2*np.pi/a;u=np.array([.03]);dk=1e4
    analytical=law.project(q,u)[2];numerical=[]
    for axis in range(3):
        shift=np.eye(3)[axis]*dk;numerical.append((law.project(q+shift,u)[0]-law.project(q-shift,u)[0])/(2*dk))
    gradient_error=float(np.max(abs(np.array(numerical).T-analytical))/max(float(abs(analytical).max()),1e-300))
    if gradient_error>1e-6:raise ValueError('Native finite-q derivative has wrong sign/units')
    zero=law.project(q,np.zeros(1));zero_error=max(float(abs(zero[i]).max()) for i in (0,2))
    if zero_error:raise ValueError('Native projection changes the zero-U acoustic susceptibility')
    table(Path(out)/'native_Gamma_X_constraints.csv',rows)
    emit(Path(out)/'native_Gamma_X_controls.json',dict(exact_Gamma_X_maximum_absolute_error_eV_A2=error,finite_q_gradient_relative_error=gradient_error,
        zero_U_projection_error=zero_error,source='f3114eb77af3343d6203a00efbd59e588cdbeb7f',
        pressure='original accepted pressure retained',incoming='accepted onsite packet birth states and archived inner histories fixed',
        low_U='Gamma uses exact primitive crossover knots; A is linearly interpolated from 0 to its first measured knot, not measured at low U',
        limitations='No unique shell partition or material vg; no regenerated tagged launch response; missing full local U(t,r) source feedback stays unresolved.'))
    return rows


def build(out,path4,ec,narrow=False):
    from medium import Medium
    from live import Strain
    out=Path(out);law=Law(path4);medium=Medium(out,native_projection=law)
    strain=CachedStrain(Strain(out,path4),out,f'E{ec:g}_')
    allocations=('onsite',*PROPAGATION_ALLOCATIONS) if narrow else PROPAGATION_ALLOCATIONS
    selected={(ec,c,a,b) for c in (('nominal_retained',) if narrow else meetings.ARMS) for a in allocations for b in (('hot_bursts',) if narrow else BEAMS)}
    prefix='first_native_' if narrow else f'nativeGX_E{ec:g}_'
    report=meetings.build(out,path4,medium.sample,strain.sample,case_selection=selected,prefix=prefix,allocations=allocations,
        source_allocation='onsite',material_factory=FrozenIncomingMaterial,extra_scientific_sources=EXTRA_SOURCES,
        event_samples=128 if narrow else meetings.EVENT_SAMPLES)
    medium.close();return report


def assemble(out):
    out=Path(out);report=meetings.merge(out,prefix='nativeGX_')
    with np.load(out/'meetings.npz') as base,np.load(out/'nativeGX_meetings.npz') as changed:
        lookup={(float(e),str(c),str(a),str(b)):i for i,(e,c,a,b) in enumerate(zip(base['E_cm_keV'],base['closure'],base['allocation'],base['beam']))}
        indices=[]
        for j,(e,c,b) in enumerate(zip(changed['E_cm_keV'],changed['closure'],changed['beam'])):
            i=lookup[float(e),str(c),'onsite',str(b)]
            for name in ('xyz_m','probe_time_s','probe_temperature_K','probe_launch_time_s','probe_lab_keV','probe_column_weight_m2','probe_time_proposal_Hz'):
                if not np.array_equal(changed[name][j],base[name][i]):raise ValueError('Frozen-launch probe queries do not match the invariant source chart: '+name)
            indices.append(i)
    np.savez_compressed(out/'nativeGX_acoustic_probe.npz',source_row_index=np.array(indices),source_acoustic_sha256=np.array(digest(out/'acoustic_probe.npz')),
        meeting_sha256=np.array(digest(out/'nativeGX_meetings.npz')),
        semantics=np.array('Exact source-row aliases after bytewise matching actual probe coordinates/times/T; DeltaK(q,0)=0 and original pressure/G3 strain coupling are held fixed.'))
    return report


def consume_all(out,path4):
    from yield_consumer import consume,standard_error
    out=Path(out);samples={};base=consume(out,path4,sample_contributions=samples)
    base_checks=json.loads((out/'consumer_checks.json').read_text())
    base_holes=list(csv.DictReader((out/'probe_projection_holes.csv').open()))
    base_projection=list(csv.DictReader((out/'probe_projection_samples.csv').open()))
    projected=consume(out,path4,prefix='nativeGX_',sample_contributions=samples)
    rows=base+[r for r in projected if r['allocation']==NAME]
    table(out/'yield_outcomes.csv',rows)
    for name,old in (('probe_projection_holes.csv',base_holes),('probe_projection_samples.csv',base_projection)):
        added=[r for r in csv.DictReader((out/('nativeGX_'+name)).open()) if r['allocation']==NAME]
        table(out/name,old+added)
    index={(float(r['E_cm_keV']),r['closure'],r['allocation'],r['beam'],r['branch'],r['comparator']):r for r in base+projected}
    controls=[]
    for r in projected:
        if r['allocation']!=NAME:continue
        key=(float(r['E_cm_keV']),r['closure'],r['beam'],r['branch'],r['comparator'])
        def full(a):return key[0],key[1],a,key[2],key[3],key[4]
        on=index[full('onsite')];coupling=index[full(FROZEN_COUPLING)]
        vector,groups=samples[full(NAME)];v_on,g_on=samples[full('onsite')];v_c,g_c=samples[full(FROZEN_COUPLING)]
        if not np.array_equal(groups,g_on) or not np.array_equal(groups,g_c):raise ValueError('Paired frozen-launch error has mismatched physical clusters')
        controls.append(dict(E_cm_keV=key[0],closure=key[1],beam=key[2],branch=key[3],comparator=key[4],
            onsite_frozen_launch_per_D=on['replacement_total_addition_per_D'],coupling_frozen_launch_per_D=coupling['replacement_total_addition_per_D'],
            native_Gamma_X_frozen_launch_per_D=r['replacement_total_addition_per_D'],
            native_minus_onsite_per_D=float((vector-v_on).sum()),native_minus_onsite_paired_MC_error_per_D=standard_error(vector-v_on,groups),
            native_minus_coupling_per_D=float((vector-v_c).sum()),native_minus_coupling_paired_MC_error_per_D=standard_error(vector-v_c,groups),
            scope='Identical incoming onsite packet states, original pressure, probe proposals and signed source chart; changed common propagation law only. Not a regenerated native launch.'))
    table(out/'frozen_launch_controls.csv',controls)
    base_checks.update(case_rows=len(rows),native_projected_meetings_sha256=digest(out/'nativeGX_meetings.npz'),
        native_projected_source_alias_sha256=digest(out/'nativeGX_acoustic_probe.npz'),
        native_projected_checks=json.loads((out/'nativeGX_consumer_checks.json').read_text()),
        native_source='frozen incoming onsite states and original pressure; source feedback is a separate tagged-map sensitivity',
        new_native_input_commit='f3114eb77af3343d6203a00efbd59e588cdbeb7f')
    emit(out/'consumer_checks.json',base_checks);return rows


def receiver_projection(out,path4):
    """Reuse the identical zero-U sound law; recalculate the actual Gamma pole."""
    from medium import Medium
    from sampling import frequency_readout,FREQUENCY_COMPARATORS
    out=Path(out);medium=Medium(out,native_projection=Law(path4))
    rows=list(csv.DictReader((out/'receiver_statistics.csv').open()));rows=[r for r in rows if r['allocation']!=NAME]
    extra=[];histories=[];source_indices=[];energies=[]
    p4,warm=FREQUENCY_COMPARATORS
    for ec in ENERGIES:
        original=list(csv.DictReader((out/f'E{ec:g}_receiver_statistics.csv').open()))
        with np.load(out/f'E{ec:g}_receiver.npz') as saved:
            for i,row in enumerate(original):
                if row['allocation']!='onsite':continue
                history=saved['history'][i];field=(history[:,3],history[:,4],np.zeros((len(history),3)),np.zeros((len(history),3)))
                fraction,zero,one,fixed,frozen=frequency_readout(medium,field,history[:,1],NAME)
                copied=next(r.copy() for r in rows if float(r['E_cm_keV'])==ec and r['closure']==row['closure'] and r['allocation']=='onsite' and r['beam']==row['beam'] and r['receiver']==row['receiver'])
                copied.update(allocation=NAME,fractional_Gamma_frequency_min=float(fraction.min()),fractional_Gamma_frequency_max=float(fraction.max()),
                    fractional_Gamma_frequency_rms=float(np.sqrt(np.mean(fraction*fraction))),plus_over_pass4=float(fraction.max()/p4),minus_over_pass4=float(-fraction.min()/p4),
                    plus_over_warm=float(fraction.max()/warm),minus_over_warm=float(-fraction.min()/warm),
                    positive_time_fraction_pass4=float(np.mean(fraction>p4)),negative_time_fraction_pass4=float(np.mean(fraction< -p4)),absolute_time_fraction_warm=float(np.mean(abs(fraction)>warm)),
                    Gamma_fixed_T_sensitivity_min=float(fixed.min()),Gamma_fixed_T_sensitivity_max=float(fixed.max()),
                    Gamma_frozen_Q_sensitivity_min=float(frozen.min()),Gamma_frozen_Q_sensitivity_max=float(frozen.max()),
                    tail_differential_absolute_error=float(np.max(abs(fraction-fixed*history[:,1]))),returned_Omega_ratio_roundoff=float(np.max(abs(one/zero-1-fraction))),
                    acoustic_source_allocation='onsite exact zero-U invariant',native_projection='f3114eb longitudinal Gamma/X, named shell-1 continuation')
                margin=float(copied['low_innovation_frequency_margin']);error=float(copied['low_innovation_Chebyshev_probability_bound']);sampling=float(copied['fraction_sampling_Hoeffding95_halfwidth'])
                copied.update(positive_fraction_Pass4_conditional_lower=max(0.,float(np.mean(fraction>p4+margin))-error-sampling),
                    positive_fraction_Pass4_conditional_upper=min(1.,float(np.mean(fraction>p4-margin))+error+sampling))
                extra.append(copied);histories.append(np.stack((history[:,0],fraction,zero,one),axis=-1));source_indices.append(i);energies.append(ec)
    table(out/'receiver_statistics.csv',rows+extra)
    np.savez_compressed(out/'native_Gamma_X_receiver.npz',history=np.array(histories),source_energy_keV=np.array(energies),source_receiver_row=np.array(source_indices),
        columns=np.array(['time_s','Gamma_fractional_frequency','strain_off_Omega_rad_s','strained_Omega_rad_s']),
        source_scope=np.array('Strain/stress/T/U use the named E*_receiver.npz onsite row; no duplicate or additional acoustic carrier.'))
    medium.close();return extra


def tagged_sensitivity(out,path4):
    """Exact retained local U snapshots; no Umax-to-mean or inverse-deltaK map."""
    from medium import Medium
    out=Path(out);law=Law(path4);medium=Medium(out,native_projection=law);rows=[];transport=[]
    accepted=meetings.Accepted(path4);lattice=accepted.json('frozen/pass1_inputs.json')['lattice'];material=meetings.Material(accepted,lattice)
    for mesh,power,label in ((16,1.,'mesh16_linear_A'),(12,1.,'mesh12_linear_A'),(16,2.,'mesh16_quadratic_low_A')):
        scenario=Law(path4,mesh,power);q=material.k
        for U in (0.,.004,.012,.026,.03,.07367988296696293):
            T=np.full(len(q),300.);u=np.full(len(q),U);zero=np.zeros((len(q),3))
            medium.native_projection=scenario
            w,v,_=medium._hamiltonian(q,T,u,zero,zero,np.zeros(len(q)),zero,NAME,False)
            quiet=medium._hamiltonian(q,T,np.zeros(len(q)),zero,zero,np.zeros(len(q)),zero,NAME,False)[0]
            for i in range(len(q)):transport.append(dict(scenario=label,U_eV_per_formula=U,mode=i,qx_m_inverse=q[i,0],qy_m_inverse=q[i,1],qz_m_inverse=q[i,2],
                delta_Omega_over_Omega=float(w[i]/quiet[i]-1),vg_x_m_s=v[i,0],vg_y_m_s=v[i,1],vg_z_m_s=v[i,2],
                scope='Computed on retained G3 k under the named Gamma/X shell-1 continuation; not a unique measured material vg'))
        for closure in meetings.ARMS:
            for sp in (0,1):
                name=f'native_E2.5_{closure}_onsite_sp{sp}_m16.npz'
                with np.load(accepted(name)) as data:
                    indices=data['selected_mode_indices'];k=data['k_m_inverse'][indices];base=data['base_K'][indices]
                    weights=data['focus_weight'];weights=weights/weights.sum();times=data['time_s'];native_q=data['variance_q_selected_A2']
                    for ti,U in zip(data['map_time_s'],data['map_U']):
                        valid=(U>=0)&(U<=min(scenario.u[-1],scenario.au[-1]));mass=float(weights[valid].sum())
                        if mass==0:raise ValueError('No supported tagged-map native excitation')
                        uk=U[valid,None];qq=k[None,:,:];new=scenario.project(qq,uk,NAME)[0]
                        old_on=scenario.project(qq,uk,'onsite')[0];old_c=scenario.project(qq,uk,'coupling')[0]
                        Q=np.array([np.interp(ti,times,native_q[:,j]) for j in range(len(indices))])
                        for j,mode in enumerate(indices):
                            wnew=material.pole(base[j]+new[:,j]+material.C*material.variance(base[j]+new[:,j],300.))
                            won=material.pole(base[j]+old_on[:,0]+material.C*material.variance(base[j]+old_on[:,0],300.))
                            wc=material.pole(base[j]+old_c[:,j]+material.C*material.variance(base[j]+old_c[:,j],300.))
                            wt=weights[valid]/mass
                            rows.append(dict(scenario=label,closure=closure,species=sp+1,time_s=ti,mode=int(mode),supported_focus_weight=mass,
                                source_U_max_eV_per_formula=float(U.max()),native_minus_onsite_frequency_fraction=float(wt@(wnew/won-1)),
                                native_minus_coupling_frequency_fraction=float(wt@(wnew/wc-1)),
                                minimum_native_minus_onsite_frequency_fraction=float(np.min(wnew/won-1)),maximum_native_minus_onsite_frequency_fraction=float(np.max(wnew/won-1)),
                                signed_frozen_covariance_potential_contrast_eV_per_D_site=float(.5*Q[j]*(wt@(new[:,j]-old_on[:,0]))),
                                scope='Accepted 2.5-keV local U snapshots and selected native covariance only; instantaneous signed Hamiltonian contrast, not integrated regenerated source work or inferred contact. Full U(t,r) missing.'))
    table(out/'native_Gamma_X_tagged_sensitivity.csv',rows);table(out/'native_Gamma_X_transport.csv',transport)
    emit(out/'native_Gamma_X_readout.json',dict(third_allocation=NAME,source_boundary='Same onsite incoming birth states and archived inner histories; original pressure and heat held fixed',
        local_U='Exact retained spatial U maps at 2.5 keV only; never reconstructed from projected DeltaK or a maximum',
        maximum_tagged_absolute_fractional_frequency_contrast=max(max(abs(r['minimum_native_minus_onsite_frequency_fraction']),abs(r['maximum_native_minus_onsite_frequency_fraction'])) for r in rows),
        source_feedback='Unresolved without full local U(t,r). Snapshot sensitivity is not a regenerated launch, donor-work time integral or physical tail bound.',
        sign='Measured folded diagonal and opposing sublattice coupling both rise; no convex 3:1 endpoint mixture. Raw 1000–1500K sign change is not a mesh-resolved crossover bracket.',
        low_A='Linear 0-to-first-measured-A continuation; mesh12 and quadratic low-A alternatives are named sensitivities.',
        material_vg='Unknown uniquely; only the named shell-1 representative is propagated.'))
    medium.close();return rows
