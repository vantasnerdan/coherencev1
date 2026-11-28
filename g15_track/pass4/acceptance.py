"""Targeted physical closed limits for the single settled joined workflow."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import ast,json,math
import numpy as np
from scipy import sparse
import fields as f
import sources as s
import run as m


def reciprocal_packet():
    """One signed packet: identical covariance Hamiltonian, retarded bath ports."""
    import native,response,heapq
    arm=s.NOMINAL;p=s.inputs()['lattice'];k=native.primitive_basis(p)
    _,_,bare,_=native.dispersion(p,k,arm['allocation']);C=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4'];A=C/3
    theta=s.KB*300.;mass=p['mass_eVc2']*s.QE/299792458.**2*1e10/s.QE
    _,quiet,_,_=native.thermal_pole(p,bare,300.);q=np.zeros(len(k));moment=q.copy();w=np.sqrt(quiet);pw=q.copy();det=quiet*theta/mass
    dt=2e-15;t=np.arange(0.,30e-12+dt/2,dt);age=t-2e-12
    epsilon=np.where((age>=0)&(age<=20e-12),1e-4*np.sin(2*np.pi*250e9*age)*np.sin(np.pi*age/20e-12)**2,0.)
    K=bare.copy();g=p['collective_strain_slope_eV_A2'];pending=[];heat=0.;bath=0.;work=0.;minimum_loss=np.inf;error=0.;records=[];stress=[]
    initial=s.RHO*4460.**2*(1e-4)**2/(3.783e28*s.QE);wave=initial
    def free():
        Q=w*w;cov=w*pw/mass;V=pw*pw/mass**2+det/Q
        _,eq,_,_=native.thermal_pole(p,K,300.)
        return float(3*np.mean(response.relative_free_energy(Q,cov,V,eq,theta/mass,theta,C)))
    for j in range(1,len(t)):
        before=free();K=bare+g*epsilon[j]/3;parameter=free()-before;wave-=parameter;work+=parameter
        before=free();q,moment,w,pw=native.gaussian_step(q,moment,w,pw,det,K,A,mass,dt*1e15)
        error=max(error,abs(free()-before));before=free()
        q,moment,w,pw,det=native.damp_population(q,moment,w,pw,det,K,A,mass,dt*1e15,1/(arm['opt']*1e15),theta)
        loss=before-free();minimum_loss=min(minimum_loss,loss)
        for fraction,tau in ((.8,arm['pd']),(.2,arm['te'])):
            energy=fraction*loss;heapq.heappush(pending,(t[j]+tau,energy));bath+=energy
        while pending and pending[0][0]<=t[j]:
            _,energy=heapq.heappop(pending);bath-=energy;heat+=energy
        error=max(error,abs(wave+free()+bath+heat-initial))
        _,eq,_,_=native.thermal_pole(p,K,300.);stress.append(float(g/2*np.mean(w*w-eq)))
        if j%10==0:records.append(dict(age_s=t[j],signed_strain=epsilon[j],dynamic_stress_eV_per_D=stress[-1],wave_eV_per_D=wave,
            D_bound_availability_eV_per_D=free(),traveling_bath_eV_per_D=bath,heat_eV_per_D=heat,balance_eV_per_D=wave+free()+bath+heat-initial))
    frequency=np.fft.rfftfreq(len(t),dt);kernel=response.covariance(arm,frequency)['bulk_per_D_J']/s.QE
    predicted=np.fft.irfft(np.fft.rfft(epsilon)*np.conj(kernel),n=len(t))
    expected=float(np.trapezoid(predicted*np.gradient(epsilon,dt),t))
    m.table('signed_wave_packet.csv',records)
    return dict(packet_initial_eV_per_D=initial,packet_signed_stress_peak=max(stress),packet_signed_stress_hole=min(stress),
        packet_minimum_wave_energy_eV_per_D=wave,packet_minimum_dissipation_eV_per_D=minimum_loss,
        packet_energy_defect_eV_per_D=error,packet_measured_work_eV_per_D=work,packet_retarded_spectral_work_eV_per_D=expected,
        packet_response_relative_difference=abs(work-expected)/max(abs(expected),1e-30),
        packet_phase_convention='exp(-i omega t); rFFT convolution uses conjugate response; loss=-omega Im(B)>0')


def main():
    result={};rng=np.random.default_rng(1981)
    # Boring positive AST query in the actual campaign Python sources.
    result['python_function_definitions']=sum(sum(isinstance(n,ast.FunctionDef) for n in ast.walk(ast.parse(p.read_text()))) for p in s.HERE.glob('*.py'))
    manifest=json.loads((s.HERE/'source_index.json').read_text())
    result['frozen_input_hashes_match']=all(s.digest(s.HERE/r['path'])==r['sha256'] for r in manifest['source_files'])
    assert result['frozen_input_hashes_match']
    # Native exact forced oscillator plus its same interaction potential.
    w=np.array([2.,3.,4.]);q=rng.normal(size=3);v=rng.normal(size=3);force=rng.normal(size=3)
    h0=.5*np.sum(v*v+(w*q)**2)-force@q
    qq,vv,work=f.oscillator(q,v,force,w,np.zeros(3),.37)
    h1=.5*np.sum(vv*vv+(w*qq)**2)-force@qq
    result['reciprocal_closed_H_relative_defect']=float(abs(h1-h0)/max(abs(h0),1.))
    assert result['reciprocal_closed_H_relative_defect']<1e-13
    # Closed thermal box: conservative conduction, no boundary heat sink.
    heat=f.Heat();heat.G=heat.G-sparse.diags(heat.exit.sum(axis=0));heat.exit[:]=0.
    u=np.zeros_like(f.VOL);u[8,8,0]=1.;new,T,ex,defect=heat.step(u,np.zeros_like(u),.001)
    result['closed_heat_energy_defect_J']=float(new.sum()-1)
    result['closed_heat_minimum_temperature_K']=float(T.min())
    assert abs(new.sum()-1)<1e-10 and T.min()>=0 and ex.sum()==0
    # Uniform-cell causality: a front moves only vg*t, including a neighbour crossing.
    a=np.zeros_like(f.VOL);a[8,8,0]=1.;delta=np.array([f.DX/4,0.,0.]);out,ex=f.continuous_shift(a,delta,1.)
    result['causal_neighbour_fraction']=float(out[9,8,0]);result['transport_closed_energy_defect']=float(out.sum()+ex-1)
    assert abs(out[9,8,0]-.25)<1e-14 and abs(out.sum()-1)<1e-14
    # Free reflection changes signed amplitude; absorption removes positive energy.
    a[:]=0.;a[-1,8,0]=1.;out,ex,signed=f.shift(a,np.array([2*f.DX,0.,0.]),.8,-1)
    result['reflected_positive_energy']=float(out.sum());result['reflected_signed_amplitude_sum']=float(signed.sum())
    assert 0<out.sum()<1 and signed.sum()<0 and abs(out.sum()+ex-1)<1e-14
    # Phase-null comparison retains energy. Loss of first moment is not loss of population.
    age=np.linspace(0,10e-12,101);kap=2e12;gamma_phi=1e12
    total=np.exp(-kap*age);coherent=total*np.exp(-2*gamma_phi*age);incoherent=total*(-np.expm1(-2*gamma_phi*age))
    result['dephasing_population_identity_defect']=float(abs(coherent+incoherent-total).max())
    assert result['dephasing_population_identity_defect']<1e-15
    # Zero excitation is exactly the once-counted Maxwell background.
    lab=np.array([.5,1.,5.]);baseline=m.moving_rate(lab,.02585,'p',.04);pump=m.pump_rate(lab,np.zeros(3),'p',.02585,.04)
    result['quiet_velocity_kernel_relative_defect']=float(np.max(abs(pump/baseline-1)))
    assert result['quiet_velocity_kernel_relative_defect']<1e-12
    # The actual G3 alias symbol and reciprocal gradients at k/-k.
    p=s.inputs()['lattice'];k=np.array([[1.1e9,2e9,-3.2e9]])
    wp,vp,Dp=s.optical(p,k);wm,vm,Dm=s.optical(p,-k)
    result['G3_reciprocal_frequency_defect']=float(abs(wp-wm).max());result['G3_reciprocal_group_defect']=float(abs(vp+vm).max())
    assert np.array_equal(Dp,Dm) and np.max(abs(vp+vm))<1e-12
    import native
    anchor=native.Material(16).reference
    for allocation in ('onsite','coupling'):
        for U in native.Material(16).u:
            _,_,D,_=native.dispersion(p,np.zeros((1,3)),allocation,U)
            assert abs(float(D[0])-(anchor+float(native.Material(16).value(U))))<1e-12
    pole,Q,_,_=native.thermal_pole(p,anchor,300.)
    C=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4']
    assert abs(Q*(anchor+C*Q)-s.KB*300)<1e-14
    result['native_bare_Gamma_K']=anchor;result['native_300K_evolved_pole_K']=float(pole)
    result.update(reciprocal_packet())
    assert result['packet_energy_defect_eV_per_D']<1e-10
    assert result['packet_response_relative_difference']<.02
    assert result['packet_minimum_dissipation_eV_per_D']> -1e-13
    result['passed']=True;m.emit('acceptance.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
