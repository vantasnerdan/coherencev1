"""One settled changed-consumer proof, including one selected dt refinement."""
from pathlib import Path
import json
import time
import numpy as np
from numpy.polynomial import Polynomial
import engine as E
import rate_mapping as R

HERE=Path(__file__).resolve().parent


def polynomial_minimum(coefficients,length):
    p=Polynomial(coefficients)
    candidates=[0.,length]
    candidates += [float(r.real) for r in p.deriv().roots() if abs(r.imag)<1e-12 and 0<r.real<length]
    return min(float(p(x)) for x in candidates)


def verify():
    result=json.loads((HERE/'results.json').read_text())
    assert result['complete']
    assert len(result['rows'])==72
    target=[r for r in result['rows'] if r['scope']=='target']
    assert len(target)==36
    assert {(r['n_D'],r['energy_per_D_eV'],r['pattern']) for r in target}=={
        (n,e,p) for n in (32,256,2048) for e in (.0006,.0048,.025)
        for p in ('converging_wave','reverse_spreading','collision','direct_optimized')}
    proofs={}
    spectral=[]
    rng=np.random.default_rng(14025)
    for m in (2,4,8):
        lattice=E.Lattice(m)
        x=rng.normal(size=(lattice.n,3))
        f=np.empty_like(x)
        energy=E.energy_force(x,f,lattice.edges,lattice.weights,lattice.near,lattice.directions,True)
        grid=np.zeros(lattice.symbol.shape+(3,));grid[tuple(lattice.grid.T)]=x
        fft_f=-np.fft.ifftn(np.fft.fftn(grid,axes=(0,1,2))*lattice.symbol[...,None],axes=(0,1,2)).real
        error=float(abs(f-fft_f[tuple(lattice.grid.T)]).max())
        assert error<3e-14
        spectral.append(dict(n_D=lattice.n,min_curvature_eV_A2=float(lattice.symbol.min()),max_curvature_eV_A2=float(lattice.symbol.max()),fft_force_error_eV_A=error))
        if m==2:
            J=np.zeros((32,32))
            for (i,j),w in zip(lattice.edges,lattice.weights):J[i,j]=J[j,i]=w
            assert np.array_equal(J,E.BASE['J'])
    proofs['spectral']=spectral
    proofs['folded_32_graph_max_error']=0.
    lattice=E.Lattice(2)
    small=rng.normal(0,1e-6,(32,3));fs=np.empty_like(small)
    enew=E.energy_force(small,fs,lattice.edges,lattice.weights,lattice.near,lattice.directions)
    eold,fold=E.M.potential(small,E.BASE)
    assert abs(enew-eold)<1e-18 and abs(fs-fold).max()<1e-13
    proofs['small_signal_old_match']={'energy_difference_eV':float(enew-eold),'max_force_difference_eV_A':float(abs(fs-fold).max())}
    one=np.zeros((32,3));one[E.TI]=.2*lattice.target_axis
    fo=np.empty_like(one);eo=E.energy_force(one,fo,lattice.edges,lattice.weights,lattice.near,lattice.directions)
    el,fl=E.WELL.energy_force(one,'wide_quadratic')
    assert abs(eo-el.sum())<1e-15 and abs(fo-fl-J@one).max()<1e-14
    proofs['one_site_limit']={'correction_energy_eV':float(eo-el.sum()),'correction_force_max_eV_A':float(abs(fo-fl-J@one).max()),
        'C1_reason':'Origin correction is A[u+]²[v+]²/r0^4+B[u+]³[v+]³/r0^6; its first derivatives vanish at either one-site limit.'}
    x=rng.normal(size=(32,3));x*=rng.uniform(.02,.75,(32,1))/np.linalg.norm(x,axis=1)[:,None]
    local_f=np.empty_like(x);ee=E.local(x,local_f,False)
    er,fr=E.WELL.energy_force(x,'wide_quadratic')
    local_error=max(abs(ee-er.sum()),float(abs(local_f-fr).max()))
    assert local_error<1e-13
    f=np.empty_like(x);value=E.energy_force(x,f,lattice.edges,lattice.weights,lattice.near,lattice.directions)
    direction=rng.normal(size=x.shape);direction/=np.linalg.norm(direction)
    h=1e-6
    fd=(lattice.energy(x+h*direction)-lattice.energy(x-h*direction))/(2*h)
    directional_error=abs(fd+np.sum(direction*f))
    assert directional_error<1e-7
    proofs['updated_well_equivalence_max_error']=local_error
    proofs['total_directional_gradient_error_eV_A']=directional_error
    cut=[]
    for r,e,slope in zip(E.QR,E.QE,E.QS):
        x=np.zeros((32,3));x[E.TI]=r*lattice.target_axis;x[E.TJ]=-r*lattice.target_axis
        f=np.empty_like(x);total=E.energy_force(x,f,lattice.edges,lattice.weights,lattice.near,lattice.directions)
        local,lf=E.WELL.energy_force(x,'wide_quadratic')
        connected=total-local.sum()
        deriv=-np.dot(f[E.TI]-f[E.TJ],lattice.target_axis)
        local_deriv=-np.dot(lf[E.TI]-lf[E.TJ],lattice.target_axis)
        assert abs(connected-e)<1e-14 and abs(deriv-local_deriv-slope)<1e-13
        cut.append(dict(q_A=float(r),connected_eV=float(connected),connected_dE_du_eV_A=float(deriv-local_deriv)))
    proofs['symmetric_cut']=cut
    # Nonlinear energy bound on the retained support: Schur proves w100>=0;
    # q>=9p proves w110>=0; w111>=0 and sum weights=1. Each radial
    # Hermite cubic is compared to K*r²/2 at all derivative roots/endpoints.
    radial_min=[]
    for a in range(3):
        for j in range(len(E.RX)-1):
            c=E.RC[a,:,j][::-1].copy();r=E.RX[j]
            c[0]-=.5*E.K*r*r;c[1]-=E.K*r;c[2]-=.5*E.K
            radial_min.append(polynomial_minimum(c,E.RX[j+1]-r))
    correction_min=[]
    for j in range(len(E.QR)-1):
        c=E.PC[:,j][::-1].copy();r=E.QR[j]
        c[0]-=E.JN*r*r;c[1]-=2*E.JN*r;c[2]-=E.JN
        correction_min.append(polynomial_minimum(c,E.QR[j+1]-r))
    low_min=polynomial_minimum([0,0,E.C4,E.C6],1.)
    assert min(radial_min)>0 and min(correction_min)>0 and low_min>=-1e-15
    assert min(E.WELL.K4)>0
    proofs['energy_bound']={
        'domain':'all site |x|<=.80 Å, nearest bond R>=1.28782409436586 Å, N=32/256/2048',
        'minimum_radial_V_minus_Kr2_over2_eV':min(radial_min),
        'minimum_native_connected_W_minus_Jq2_eV':min(correction_min),
        'minimum_origin_connected_correction_eV':low_min,
        'conclusion':'V>=lambda_min*sum(|x|²)/2>=0 and Krel<=Ktotal<=initial N e in conservative supported motion',
        'evidence_tier':'analytic angular inequalities plus binary64 polynomial extrema, not outward-rounded formal enclosure'}
    launches=np.load(HERE/'launches.npz')
    for row in target:
        assert row['candidate_supported'] and row['peak_before_support_exit']
        assert row['peak_relative_kinetic_eV'] <= row['n_D']*row['energy_per_D_eV']*(1+max(.002,row['relative_energy_drift']))
        if row['pattern']=='direct_optimized':
            key=f"n{row['n_D']}_e{row['energy_per_D_eV']:g}_{row['pattern']}"
            q,v=launches[key+'_q'],launches[key+'_v']
            assert not q.any() and not v[[E.TI,E.TJ]].any()
            max_speed=2*np.sqrt(2*row['energy_per_D_eV']/E.MASS)*np.sqrt(row['n_D']/(row['n_D']-2))
            assert np.linalg.norm(v,axis=1).max()<=max_speed*(1+2e-12)
    lattice=E.Lattice(8);key='n2048_e0.025_direct_optimized'
    q,v=launches[key+'_q'],launches[key+'_v']
    started=time.monotonic()
    refined=E.trajectory(lattice,q,v,dt=.25,horizon=result['horizon_fs'])
    rr=E.summarize(lattice,q,v,refined,.025,'direct_optimized','single_dt_refinement')
    original=next(r for r in target if r['n_D']==2048 and r['energy_per_D_eV']==.025 and r['pattern']=='direct_optimized')
    proofs['single_limiting_step_refinement']={
        'n_D':2048,'energy_per_D_eV':.025,'dt_fs':[.5,.25],
        'coarse_peak_eV':original['peak_relative_kinetic_eV'],'fine_peak_eV':rr['peak_relative_kinetic_eV'],
        'fractional_peak_change':rr['peak_relative_kinetic_eV']/original['peak_relative_kinetic_eV']-1,
        'coarse_peak_time_fs':original['peak_time_fs'],'fine_peak_time_fs':rr['peak_time_fs'],
        'coarse_hold_fs':original['half_peak_hold_fs'],'fine_hold_fs':rr['half_peak_hold_fs'],
        'coarse_relative_energy_drift':original['relative_energy_drift'],'fine_relative_energy_drift':rr['relative_energy_drift'],
        'coarse_first_exit_fs':original['first_site_or_bond_exit_fs'],'fine_first_exit_fs':rr['first_site_or_bond_exit_fs'],
        'elapsed_s':time.monotonic()-started}
    tr=refined[2] if refined[3]<0 else refined[2][:-1]
    mask=tr[:,4]>.03
    ratios=tr[mask,10]/tr[mask,11]
    at_max=int(np.argmax(tr[:,4]))
    wall=dict(sample='actual most-displaced site at every .25-fs refinement step',
        minimum_sampled_radius_A=float(tr[mask,4].min()),maximum_sampled_radius_A=float(tr[:,4].max()),
        local_wall_over_old_harmonic_min=float(ratios.min()),local_wall_over_old_harmonic_max=float(ratios.max()),
        updated_local_eV_at_largest_radius=float(tr[at_max,10]),old_harmonic_eV_at_largest_radius=float(tr[at_max,11]),
        min_bond_A=float(tr[:,5].min()))
    proofs['sampled_wall_comparison']=wall
    proofs['wall_interpretation']=(
        f"The updated local wall costs {wall['local_wall_over_old_harmonic_min']:.3g}–"
        f"{wall['local_wall_over_old_harmonic_max']:.3g} times the old harmonic wall at the actual most-displaced site's sampled directions "
        f"above .03 Å (maximum radius {wall['maximum_sampled_radius_A']:.4f} Å). "
        f"The .25-fs trajectory {'reaches the finite horizon without a support exit' if rr['first_site_or_bond_exit_fs'] is None else 'hits the support boundary before the finite horizon'}; "
        "neither a near-end peak nor a first boundary crossing establishes saturation. The matched high-energy harmonic trajectory exits early, so its smaller supported peak is clipping, not a same-horizon nonlinear advantage.")
    proofs['refined_selected_row']=rr
    coarse=np.load(HERE/'trajectories.npz')
    ct=coarse[key]
    coarse_time=original['peak_time_fs'];fine_time=rr['peak_time_fs']
    near_coarse=tr[np.abs(tr[:,0]-coarse_time)<=30]
    near_fine=ct[np.abs(ct[:,0]-fine_time)<=30]
    proofs['single_limiting_step_refinement'].update(
        fine_peak_in_coarse_peak_window_eV=float(near_coarse[:,1].max()) if len(near_coarse) else None,
        coarse_peak_in_fine_peak_window_eV=float(near_fine[:,1].max()) if len(near_fine) else None,
        fine_min_bond_A=rr['min_bond_A_before_exit'],coarse_min_bond_A=original['min_bond_A_before_exit'],
        fine_max_site_radius_A=rr['max_site_radius_A_before_exit'],coarse_max_site_radius_A=original['max_site_radius_A_before_exit'],
        window_half_width_fs=30.,
        interpretation='one same-launch refinement only; exact rows retained, not a convergence claim')
    trace_payload={name:coarse[name] for name in coarse.files}
    trace_payload[key+'_dt025']=refined[2][::8]
    coarse.close()
    np.savez_compressed(HERE/'trajectories.npz',**trace_payload)
    spectra=R.load_contact_spectra()
    for row,stored in zip(target,result['conditional_readouts']):
        assert R.conditional_readout(row,spectra)==stored['readout']
    proofs['changed_consumer']={'target_rows':36,'whole_cell_rows':36,'all_36_readouts_match_retained_spectra':True,
        'candidate_score_only_pre_exit':True,'direct_target_rows_zero_and_site_cap_checked':True}
    (HERE/'consumer_proof.json').write_text(json.dumps(proofs,indent=2,allow_nan=False)+'\n')
    print(json.dumps(proofs,indent=2))

if __name__=='__main__':verify()
