"""Installed i-PI closed-path trajectories and Gaussian variational references.

Physical position uses bead coordinates. Fictitious bead velocities are never
exported as nuclear momenta. The occupied periodic graph is stored in models.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from xml.sax.saxutils import escape
import numpy as np
from ase.units import Bohr,Hartree
from scipy.optimize import root_scalar
# i-PI loads the custom PES by file path from its temporary run directory.
sys.path.insert(0,str(Path(__file__).resolve().parent))
from model import HERE,KB,MASS,full,partial_cases,strain_cases,potential,kinetic_coefficient,serializable

__DRIVER_CLASS__='BulkDriver'
__DRIVER_NAME__='bulk-site'


class BulkDriver:
    def __init__(self,model_path,**kwargs):
        self.model=json.loads(Path(model_path).read_text())
        for key in ('K','f','J'):self.model[key]=np.array(self.model[key])
        self.virial=np.zeros((3,3))

    def __call__(self,cell,positions):
        x=np.asarray(positions)*Bohr
        # ffdirect passes flattened atom coordinates, optionally batch-first.
        x=x.reshape(x.shape[:-1]+(-1,3)) if x.shape[-1]!=3 else x
        energies,forces=potential(x,self.model)
        if np.ndim(energies):
            return [(float(e/Hartree),f.reshape(-1)*Bohr/Hartree,self.virial,'{}') for e,f in zip(energies,forces)]
        return float(energies/Hartree),forces.reshape(-1)*Bohr/Hartree,self.virial,'{}'


def read_xyz(path,N):
    lines=Path(path).read_text().splitlines();length=N+2
    if len(lines)%length:raise ValueError('Incomplete i-PI trajectory frame')
    rows=[np.array([[float(v) for v in line.split()[1:4]] for line in lines[i+2:i+length]]) for i in range(0,len(lines),length)]
    return np.array(rows)


def block_stats(values,blocks=12):
    values=np.asarray(values);n=len(values)//blocks
    means=values[-n*blocks:].reshape((blocks,n)+values.shape[1:]).mean(axis=1)
    return means.mean(axis=0),means.std(axis=0,ddof=1)/np.sqrt(blocks)


def observable_stats(values,logweights=None):
    """Contiguous-block expectation/SE, retaining ratio covariance for weights."""
    values=np.asarray(values)
    if logweights is None:return block_stats(values)
    weights=np.exp(logweights-logweights.max());bs=len(values)//12;used=12*bs
    shape=(12,bs)+(1,)*(values.ndim-1);w=weights[:used].reshape(12,bs)
    num=(values[:used].reshape((12,bs)+values.shape[1:])*w.reshape(shape)).mean(axis=1)
    den=w.mean(axis=1);value=num.sum(axis=0)/den.sum()
    linear=value+(num-den.reshape((12,)+(1,)*value.ndim)*value)/den.mean()
    return value,linear.std(axis=0,ddof=1)/np.sqrt(12)


def trajectory(model,isotope,T,P,steps=18000,stride=15):
    N=len(model['K']);mass=MASS[isotope]/(510998.95069/931494102.42)
    with tempfile.TemporaryDirectory(prefix='g10-ipi-') as tmp:
        work=Path(tmp);(work/'model.json').write_text(json.dumps(serializable(model)))
        (work/'init.xyz').write_text(f'{N}\n# CELL(abcABC): 40 40 40 90 90 90 positions{{angstrom}} cell{{angstrom}}\n'+''.join('H 0 0 0\n' for _ in range(N)))
        (work/'input.xml').write_text(f'''<simulation verbosity="low">
<output prefix="bulk"><properties stride="{stride}" filename="out">[step,potential{{electronvolt}}]</properties>
<trajectory stride="{stride}" filename="positions" format="xyz" bead="-1">positions{{angstrom}}</trajectory></output>
<total_steps>{steps}</total_steps><prng><seed>731</seed></prng>
<ffdirect name="bulk"><pes>custom</pes><pes_path>{escape(str(Path(__file__).resolve()))}</pes_path><batch_size>{P}</batch_size>
<parameters>{{model_path:{work/'model.json'}}}</parameters></ffdirect>
<system><initialize nbeads="{P}"><file mode="xyz">init.xyz</file><masses mode="manual">{[mass]*N}</masses>
<velocities mode="thermal" units="kelvin">{T}</velocities></initialize>
<forces><force forcefield="bulk"/></forces><ensemble><temperature units="kelvin">{T}</temperature></ensemble>
<motion mode="dynamics"><fixcom>False</fixcom><dynamics mode="nvt"><timestep units="femtosecond">0.35</timestep>
<thermostat mode="pile_l"><tau units="femtosecond">20</tau></thermostat></dynamics></motion></system></simulation>''')
        proc=subprocess.run([str(Path(sys.executable).with_name('i-pi')),'input.xml'],cwd=work,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        if proc.returncode:raise RuntimeError(proc.stdout[-7000:])
        files=sorted(work.glob('bulk.positions_*.xyz'))
        if len(files)!=P:raise RuntimeError(f'Expected {P} bead files, found {[p.name for p in work.iterdir()]}')
        paths=np.stack([read_xyz(p,N) for p in files],axis=1)
        # Discard initial 1/3 of the chain, retain complete contiguous blocks.
        paths=paths[len(paths)//3:]
        properties=np.loadtxt(work/'bulk.out')
    site=model['site'];selected=paths[:,:,site,:]
    field=np.einsum('j,fpja->fpa',model['J'][site],paths)
    _,force=potential(paths,model)
    # Centroid virial estimator: <K>=3N/(2 beta) - .5 <(x-xc).F>.
    kinetic=3/(2/(KB*T))-np.mean(np.sum((paths-paths.mean(axis=1,keepdims=True))*force,axis=-1),axis=(1,2))/2
    moments=np.mean(selected[:,:, :,None]*selected[:,:,None,:],axis=1)
    mean,se=block_stats(moments);center,cse=block_stats(selected.mean(axis=1));km,kse=block_stats(kinetic)
    weights={}
    if model['name'] in ('Pd','Ni'):
        for variant in ('energy','force'):
            other=full(model['name'],variant)
            delta=potential(paths,other)[0]-potential(paths,model)[0]
            weights[variant]=-delta.mean(axis=1)/(KB*T)
    return selected,field,weights,dict(temperature_K=T,isotope=isotope,beads=P,steps=steps,stride=stride,
        saved_frames=len(selected),mean_A=center.tolist(),mean_standard_error_A=cse.tolist(),
        second_moment_A2=mean.tolist(),second_moment_standard_error_A2=se.tolist(),
        mean_square_radius_A2=float(np.trace(mean)),mean_kinetic_per_site_eV=float(km),kinetic_standard_error_eV=float(kse),
        coordinate_support_A=model['support_A'],position_outside_support_fraction=float(np.mean(np.linalg.norm(selected,axis=-1)>model['support_A'])),
        sampling='closed primitive paths; equal-time bead marginal; 12 contiguous time blocks')


def gaussian_lattice(model,mass,T,host_mass=None):
    """Gaussian variational stationarity on the same unclamped graph."""
    N=len(model['K']);K=model['K'].copy();H0=np.block([[K[i] if i==j else -model['J'][i,j]*np.eye(3) for j in range(N)] for i in range(N)])
    # Positive quartic stabilizes soft bare graphs; initialization is a trial
    # variational matrix, not clipping or alteration of the physical Hessian.
    trial=H0+np.eye(3*N)*max(1.,1.-np.linalg.eigvalsh(H0)[0])
    masses=np.full(N,mass if host_mass is None else host_mass);masses[model['site']]=mass
    rootc=np.repeat(np.sqrt([2*kinetic_coefficient(v) for v in masses]),3)
    center=np.zeros((N,3))
    for iteration in range(250):
        e,U=np.linalg.eigh(rootc[:,None]*trial*rootc[None,:])
        if e.min()<=0:raise ValueError('Variational branch lost positive trial covariance')
        q=np.sqrt(e);v=1/(2*q*np.tanh(q/(2*KB*T)))
        weighted=rootc[:,None]*U
        C=(weighted*v)@weighted.T;local=np.array([C[3*i:3*i+3,3*i:3*i+3] for i in range(N)])
        diagonal=np.diagonal(local,axis1=1,axis2=2)+center**2
        hess=K.copy()
        for axis in range(3):hess[:,axis,axis]+=3*model['a']*diagonal[:,axis]+model['b']*(diagonal.sum(axis=1)-diagonal[:,axis])
        for a in range(3):
            for b in range(3):
                if a!=b:hess[:,a,b]+=2*model['b']*(local[:,a,b]+center[:,a]*center[:,b])
        new=np.block([[hess[i] if i==j else -model['J'][i,j]*np.eye(3) for j in range(N)] for i in range(N)])
        # Analytic Gaussian averaged force, including nonzero registered F0.
        meanforce=model['f']-np.einsum('iab,ib->ia',K,center)+model['J']@center
        meanforce-=model['a']*(center**3+3*center*np.diagonal(local,axis1=1,axis2=2))
        for a in range(3):
            for b in range(3):
                if a!=b:meanforce[:,a]-=model['b']*(center[:,a]*(center[:,b]**2+local[:,b,b])+2*center[:,b]*local[:,a,b])
        shift=np.linalg.solve(new,meanforce.ravel()).reshape(N,3)
        error=max(abs(new-trial).max(),abs(meanforce).max())
        trial=.5*(trial+new);center+=.5*shift
        if error<1e-9:break
    if error>=1e-8:raise RuntimeError('Gaussian variational stationarity unresolved')
    site=model['site']
    return dict(case=model['name'],temperature_K=T,mass_amu=mass,host_mass_amu=mass if host_mass is None else host_mass,iterations=iteration+1,
                local_centroid_A=center[site].tolist(),local_covariance_A2=local[site].tolist(),
                local_trial_stiffness_eV_A2=hess[site].tolist(),
                min_bare_graph_eigenvalue_eV_A2=float(np.linalg.eigvalsh(H0)[0]),
                min_trial_eigenvalue_eV_A2=float(np.linalg.eigvalsh(trial)[0]),
                stationarity_residual=float(error),scope='analytic Gaussian variational reference, not a non-Gaussian tail')


def sscha_local(model,mass,T):
    """Installed SSCHA API, translational dimer with site relative mass m."""
    from ase import Atoms
    from ase.calculators.calculator import Calculator,all_changes
    from cellconstructor import Structure,Phonons,Units
    from sscha.Ensemble import Ensemble
    k=float(model['K'][0,0,0]);a=model['a'];b=model['b'];c=2*kinetic_coefficient(mass)
    def residual(trial):
        quantum=np.sqrt(c*trial);v=c/(2*quantum*np.tanh(quantum/(2*KB*T)))
        return trial-k-(3*a+2*b)*v
    trial=root_scalar(residual,bracket=(max(k,.001),10.)).root
    sep=np.array([2.,0,0]);atoms=Atoms('HH',positions=[[0,0,0],sep],cell=[20]*3,pbc=True,masses=[2*mass]*2)
    structure=Structure.Structure();structure.generate_from_ase_atoms(atoms)
    projector=np.kron([[1.,-1.],[-1.,1.]],np.eye(3));dyn=Phonons.Phonons(structure,nqirr=1,force_real=True)
    unit=Units.RY_TO_EV*Units.A_TO_BOHR**2;dyn.dynmats[0]=trial*projector/unit
    local=dict(model,K=np.array([np.eye(3)*k]),f=np.zeros((1,3)),J=np.zeros((1,1)))
    class Calc(Calculator):
        implemented_properties=['energy','forces']
        def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
            super().calculate(atoms,properties,system_changes)
            x=self.atoms.positions[1]-self.atoms.positions[0]-sep
            energy,force=potential(x[None,:],local)
            self.results=dict(energy=float(energy),forces=np.array([-force[0],force[0]]))
    np.random.seed(587);ens=Ensemble(dyn,T,supercell=(1,1,1));ens.generate(2048,evenodd=True)
    ens.compute_ensemble(Calc(),compute_stress=False,verbose=False)
    # Pinned API order: force evaluation -> gradient -> perturb forces -> restore
    # -> Hessian (which internally converts ensemble units). Never reverse it.
    gradient,error=ens.get_preconditioned_gradient(return_error=True,verbose=False)
    saved=ens.forces.copy()  # Installed default storage is Ry/angstrom.
    ens.forces=np.round(saved*Units.RY_TO_EV,4)/Units.RY_TO_EV
    rounded=ens.get_preconditioned_gradient(verbose=False)
    ens.forces=saved
    hessian=ens.get_free_energy_hessian(include_v4=True,get_full_hessian=True,use_symmetries=False,verbose=False)
    matrix=hessian.dynmats[0].real*unit
    return dict(temperature_K=T,mass_amu=mass,analytic_SCHA_stiffness_eV_A2=trial,
                SSCHA_gradient_stiffness_residual_eV_A2=float(np.sum(projector*gradient.real)*unit/np.sum(projector**2)),
                SSCHA_gradient_error_norm_eV_A2=float(np.linalg.norm(error)*unit),
                force_rounding_1e_minus4_eV_A_gradient_change_eV_A2=float(np.linalg.norm(rounded-gradient)*unit),
                free_energy_hessian_relative_eV_A2=((matrix[:3,:3]+matrix[3:,3:]-matrix[:3,3:]-matrix[3:,:3])/4).tolist(),
                hessian_convention='relative coordinate stiffness = dimer nonzero eigenvalue/2 = four-block contraction/4',
                ensemble_configurations=2048,scope='installed SSCHA isolated-site reference; coupled lattice Gaussian stationarity separately')


def summarize_saved(x,field,model,isotope,T,logweights=None):
    """Recover the selected-site consumer from already produced paths."""
    site=model['site'];squared=x*x
    force=-np.einsum('ab,fpb->fpa',model['K'][site],x)+model['f'][site]+field
    force-=((model['a']-model['b'])*squared+model['b']*squared.sum(axis=-1,keepdims=True))*x
    kinetic=1.5*KB*T-.5*np.mean(np.sum((x-x.mean(axis=1,keepdims=True))*force,axis=-1),axis=1)
    stats=lambda values:observable_stats(values,logweights)
    mean,se=stats(np.mean(x[:,:,:,None]*x[:,:,None,:],axis=1))
    center,cse=stats(x.mean(axis=1));km,kse=stats(kinetic)
    outside,_=stats(np.mean(np.linalg.norm(x,axis=-1)>model['support_A'],axis=1))
    _,r2se=stats(np.mean(np.sum(x*x,axis=-1),axis=1))
    return dict(temperature_K=T,isotope=isotope,beads=x.shape[1],steps=18000,stride=15,saved_frames=len(x),
                mean_A=center.tolist(),mean_standard_error_A=cse.tolist(),second_moment_A2=mean.tolist(),
                second_moment_standard_error_A2=se.tolist(),mean_square_radius_A2=float(np.trace(mean)),
                r2_standard_error_A2=float(r2se),
                mean_kinetic_per_site_eV=float(km),kinetic_standard_error_eV=float(kse),
                kinetic_scope='selected-site centroid virial from saved local path and actual neighbor field',
                coordinate_support_A=model['support_A'],
                position_outside_support_fraction=float(outside),
                sampling='closed primitive paths; equal-time bead marginal; 12 contiguous time blocks')


def tracer_paths(source,model,isotope,T):
    """Exact finite-P centroid-scaled mass importance sampling, fixed D host.

    q'=q_centroid+sqrt(m_D/m_tracer)*(q-q_centroid) preserves the selected
    spring action. Its constant Jacobian cancels in normalized expectations;
    the remaining weight is exp[-beta/P sum(V(q')-V(q))].
    """
    from conditioned import quartic
    x=source['beads_A'];field=source['neighbor_field_eV_A'];site=model['site']
    centroid=x.mean(axis=1,keepdims=True)
    changed=centroid+np.sqrt(MASS['D']/MASS[isotope])*(x-centroid)
    def local(y):
        return .5*np.einsum('fpa,ab,fpb->fp',y,model['K'][site],y)-np.sum(y*(model['f'][site]+field),axis=-1)+quartic(y,model['a'],model['b'])
    logmass=-(local(changed)-local(x)).mean(axis=1)/(KB*T)
    weights={'joint':logmass}
    for variant in ('energy','force'):
        v=full('Pd',variant)
        old_delta=quartic(x,v['a'],v['b'])-quartic(x,model['a'],model['b'])
        new_delta=quartic(changed,v['a'],v['b'])-quartic(changed,model['a'],model['b'])
        weights[variant]=logmass+source['log_weights_'+variant]-(new_delta-old_delta).mean(axis=1)/(KB*T)
    return changed,field,weights


def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    read_rows=lambda name:json.loads((output/name).read_text()) if (output/name).exists() else []
    rows=read_rows('trajectory_summary.json');gaussian=read_rows('gaussian_reference.json');sscha=read_rows('sscha_reference.json')
    def save():
        for name,data in [('trajectory_summary.json',rows),('gaussian_reference.json',gaussian),('sscha_reference.json',sscha)]:
            (output/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    jobs=[(full('Pd'),'D',T) for T in (300,77)]+[(full('Ni'),'H',T) for T in (300,77)]
    for m,i,T in jobs:
        P=32 if T==300 else 128;name=f'{m["name"]}{i}_{T}K';path=output/(name+'_paths.npz')
        if not any(r['case']==name for r in rows):
            if path.exists():
                with np.load(path) as z:row=summarize_saved(z['beads_A'],z['neighbor_field_eV_A'],m,i,T)
            else:
                x,field,weights,row=trajectory(m,i,T,P)
                np.savez_compressed(path,beads_A=x,neighbor_field_eV_A=field,
                                    **{'log_weights_'+k:v for k,v in weights.items()})
            row.update(case=name,model=m['name'],exact_case=name,isotope_semantics='homogeneous isotope replacement',host_isotope=i);rows.append(row);save()
        if not any(r['case']==m['name'] and r['temperature_K']==T and r['mass_amu']==MASS[i] for r in gaussian):
            gaussian.append(gaussian_lattice(m,MASS[i],T));save()
        if m['name'] in ('Pd','Ni') and not any(r['case']==name for r in sscha):
            sscha.append(dict(case=name,**sscha_local(m,MASS[i],T)));save()
        print(name,next(r['mean_square_radius_A2'] for r in rows if r['case']==name),flush=True)
    # H/T are explicit tracers in the sampled D host, not homogeneous PdH/PdT.
    for T in (300,77):
        m=full('Pd')
        with np.load(output/f'PdD_{T}K_paths.npz') as source:
            for isotope in ('H','T'):
                name=f'Pd{isotope}_Dhost_{T}K';path=output/(name+'_paths.npz')
                if not any(r['case']==name for r in rows):
                    x,field,weights=tracer_paths(source,m,isotope,T)
                    np.savez_compressed(path,beads_A=x,neighbor_field_eV_A=field,**{'log_weights_'+k:v for k,v in weights.items()})
                    row=summarize_saved(x,field,m,isotope,T,weights['joint'])
                    w=np.exp(weights['joint']-weights['joint'].max())
                    row.update(case=name,model='Pd',exact_case=f'Pd{isotope}_{T}K',host_isotope='D',
                               isotope_semantics='selected-site mass tracer in fixed D host; exact finite-P centroid-scaled importance weights',
                               frame_weight_ess=float(w.sum()**2/(w@w)),source_case=f'PdD_{T}K')
                    rows.append(row);save()
                if not any(r.get('trajectory_case')==name for r in gaussian):
                    gaussian.append(dict(gaussian_lattice(m,MASS[isotope],T,host_mass=MASS['D']),trajectory_case=name));save()
                if not any(r['case']==name for r in sscha):
                    sscha.append(dict(case=name,**sscha_local(m,MASS[isotope],T)));save()
    reference=output/'PdH_300K_hostH_reference_paths.npz'
    if reference.exists() and not any(r['case']=='PdH_300K_hostH_reference' for r in rows):
        with np.load(reference) as source:row=summarize_saved(source['beads_A'],source['neighbor_field_eV_A'],full('Pd'),'H',300)
        row.update(case='PdH_300K_hostH_reference',model='Pd',exact_case='PdH_300K',host_isotope='H',isotope_semantics='genuinely sampled homogeneous H host reference')
        rows.append(row);save()
    # Partial/strain contrast is variational and separate, never a fabricated
    # non-Gaussian distribution or a measured local strained well.
    contrasts=[gaussian_lattice(m,MASS['D'],300) for m in partial_cases()+strain_cases()]
    (output/'model_contrasts.json').write_text(json.dumps(contrasts,indent=2,allow_nan=False)+'\n')
    # Same physical PdD consumer at half P isolates discretization, not input.
    for T,P in ((300,16),(77,64)):
        m=full();name=f'PdD_{T}K_P{P}';path=output/(name+'_paths.npz')
        if any(r['case']==name for r in rows):continue
        if path.exists():
            with np.load(path) as z:row=summarize_saved(z['beads_A'],z['neighbor_field_eV_A'],m,'D',T)
        else:
            x,field,w,row=trajectory(m,'D',T,P)
            np.savez_compressed(path,beads_A=x,neighbor_field_eV_A=field,**{'log_weights_'+k:v for k,v in w.items()})
        row.update(case=name,model='Pd');rows.append(row);save()
    return rows
