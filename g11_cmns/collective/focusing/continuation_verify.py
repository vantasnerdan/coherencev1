"""One settled consumer check; no previous contact/rate calculation is rerun."""
import json
import numpy as np
import engine as E
import exterior as X
from continuation import HERE,load,save


def force(lattice,x,stiff=None):
    f = np.empty_like(x)
    args = (x,f,lattice.edges,lattice.weights,lattice.near,lattice.directions)
    energy = E.energy_force(*args) if stiff is None else X.energy_force(*args,stiff)
    return energy,f


def gradient_defect(lattice,state,stiff,rows):
    _,f = force(lattice,state,stiff)
    defect = 0.
    for i in sorted(rows):
        for a in range(3):
            xp,xm = state.copy(),state.copy(); h = 2e-6
            xp[i,a] += h; xm[i,a] -= h
            derivative = (force(lattice,xp,stiff)[0]-force(lattice,xm,stiff)[0])/(2*h)
            defect = max(defect,abs(derivative+f[i,a]))
    return float(defect)


def verify():
    data,traces = load()
    for key,case in data['cases'].items():
        if key.startswith('outside_'):
            case['trace_columns'] = data['trace_columns'][:10]
            for row in case['candidates']:
                name = key+'__'+row['candidate']
                traces[name] = traces[name][:,:10]
    expected = {'size_N16384':14,'size_N131072':14,'fine_N2048_25meV':11,
                'outside_soft_4.8meV':11,'outside_stiff_4.8meV':11,
                'outside_soft_25meV':11,'outside_stiff_25meV':11}
    for key,count in expected.items():
        case = data['cases'][key]
        assert len(case['candidates'])==count
        for row in case['candidates']:
            trace = traces[key+'__'+row['candidate']]
            assert trace.shape[1]==len(case.get('trace_columns',data['trace_columns']))
            assert np.isfinite(trace).all()
            assert row['target_kinetic_preload_fraction']==0
            assert abs(row['launch_total_energy_eV']/(case['n_D']*case['energy_per_D_eV'])-1)<1e-10
            if not row['candidate'].startswith('wave_'):
                assert row['launch_max_site_kinetic_eV']<=4*case['energy_per_D_eV']*case['n_D']/(case['n_D']-2)*(1+1e-12)
        for row in case['selected']:
            options = [r for r in case['candidates'] if r['pattern']==row['pattern']]
            valid = [r for r in options if r.get('model_valid_full_horizon',True)]
            eligible = [r for r in valid if r['headline_eligible']]
            assert row['candidate']==max(eligible or valid or options,key=lambda r:r['peak_relative_kinetic_eV'])['candidate']
    lattice = E.Lattice(2)
    rng = np.random.default_rng(20261125)
    x = rng.normal(size=(lattice.n,3))
    x *= (.72*rng.random(lattice.n)/np.linalg.norm(x,axis=1))[:,None]
    baseline,fbase = force(lattice,x)
    equality = []
    for stiff in (False,True):
        energy,f = force(lattice,x,stiff)
        assert energy==baseline and np.array_equal(f,fbase)
        equality.append(dict(model='stiff' if stiff else 'soft',energy_bitwise_equal=True,force_bitwise_equal=True))
    q = np.zeros_like(x); v = rng.normal(size=x.shape)
    v *= np.sqrt(2*lattice.n*.0006/(E.MASS*np.sum(v*v)))
    original = E.trajectory(lattice,q,v,.25,50.)
    trajectory_errors = []
    for stiff in (False,True):
        extended = X.trajectory(lattice,q,v,.25,50.,stiff)
        assert np.array_equal(original[0],extended[0]) and np.array_equal(original[1],extended[1])
        error = float(np.max(abs(original[2][:,:10]-extended[2][:,:10])))
        assert error<1e-12
        trajectory_errors.append(error)
    # C1 checks on the full angular energy across a single-site boundary and
    # across the connected q boundary. These are declared synthetic probes.
    boundary_errors = []
    for stiff in (False,True):
        for direction in (np.array([1.,0.,0.]),np.array([1.,2.,3.])/np.sqrt(14)):
            xm = np.zeros_like(x); xp = xm.copy()
            xm[0] = (.8-1e-8)*direction; xp[0] = (.8+1e-8)*direction
            em,fm = force(lattice,xm,stiff); ep,fp = force(lattice,xp,stiff)
            boundary_errors.append(dict(stiff=stiff,energy_jump=float(abs(ep-em)),force_jump=float(np.max(abs(fp-fm)))))
            assert abs(ep-em)<1e-5 and np.max(abs(fp-fm))<1e-5
        xm = np.zeros_like(x); xp = xm.copy()
        xm[list(lattice.target)] = np.array([1.,-1.])[:,None]*(.8-1e-8)*lattice.target_axis
        xp[list(lattice.target)] = np.array([1.,-1.])[:,None]*(.8+1e-8)*lattice.target_axis
        em,fm = force(lattice,xm,stiff); ep,fp = force(lattice,xp,stiff)
        boundary_errors.append(dict(stiff=stiff,connected_pair=True,
            energy_jump=float(abs(ep-em)),force_jump=float(np.max(abs(fp-fm)))))
        assert abs(ep-em)<1e-5 and np.max(abs(fp-fm))<1e-5
    synthetic_exterior = []
    synthetic = np.zeros_like(x)
    synthetic[0] = 1.05*np.array([1.,2.,3.])/np.sqrt(14)
    synthetic[list(lattice.target)] = np.array([1.,-1.])[:,None]*.92*lattice.target_axis
    for stiff in (False,True):
        defect = gradient_defect(lattice,synthetic,stiff,{0,*lattice.target})
        assert defect<2e-6
        synthetic_exterior.append(dict(stiff=stiff,site_radius_A=1.05,target_q_A=.92,
            force_gradient_max_abs_eV_A=defect,scope='declared synthetic model check, not an observed trajectory state'))
    lattice = E.Lattice(8)
    states = dict(np.load(HERE/'continuation_states.npz'))
    outside = [c['selected'][0] for k,c in data['cases'].items() if k.startswith('outside_')]
    worst = max(outside,key=lambda r:r['max_site_radius_A_before_exit'])
    key = next(k for k,c in data['cases'].items() if k.startswith('outside_') and c['selected'][0] is worst)
    k = int(np.argmax([g['largest_site']['radius_A'] for g in worst['simultaneous_geometry']]))
    state = states[key+'__x'][k]
    stiff = '_stiff_' in key
    geometry = worst['simultaneous_geometry'][k]
    rows = {geometry['largest_site']['row'],geometry['minimum_NN_pair']['first']['row'],
            geometry['minimum_NN_pair']['second']['row'],*lattice.target}
    defect = gradient_defect(lattice,state,stiff,rows)
    assert defect<2e-6
    sanity = dict(case=key,time_fs=geometry['time_fs'],max_radius_A=geometry['largest_site']['radius_A'],
        finite_difference_force_max_abs_eV_A=defect,
        outside_site_extent=geometry['largest_site']['radius_A']>.8)
    actual_outer = []
    for key,case in data['cases'].items():
        if not key.startswith('outside_'):
            continue
        event = case['search_extrema'][0]
        g = event['geometry']
        if g['largest_site']['radius_A']<=.8:
            continue
        xx = states[key+'__largest_site_over_search__x']
        selected_rows = {g['largest_site']['row'],g['minimum_NN_pair']['first']['row'],
                         g['minimum_NN_pair']['second']['row'],*lattice.target}
        defect = gradient_defect(lattice,xx,'_stiff_' in key,selected_rows)
        assert defect<2e-6
        candidate = next(r for r in case['candidates'] if r['candidate']==event['candidate'])
        actual_outer.append(dict(case=key,candidate=event['candidate'],time_fs=g['time_fs'],
            radius_A=g['largest_site']['radius_A'],force_gradient_max_abs_eV_A=defect,
            relative_energy_drift=candidate['relative_energy_drift'],
            model_failure=candidate['model_failure']))
    soft = traces['outside_soft_25meV__direct_1500fs']
    stiff_trace = traces['outside_stiff_25meV__direct_1500fs']
    first_exit = next(r['first_site_or_bond_exit_fs'] for r in
        data['cases']['outside_soft_25meV']['candidates'] if r['candidate']=='direct_1500fs')
    prefix = int(round(first_exit/.25))
    assert np.array_equal(soft[:prefix],stiff_trace[:prefix])
    assert np.max(abs(soft[prefix:,4]-stiff_trace[prefix:,4]))>1e-4
    result = dict(scope='settled ordered continuation consumers only; no legacy trajectory or rate replay',
        cases_checked=expected, on_domain_energy_force=equality,
        on_domain_final_x_v_bitwise_equal=True,on_domain_trace_max_abs=trajectory_errors,
        synthetic_C1_boundary_errors=boundary_errors,selected_state_sanity=sanity,
        synthetic_outside_force_checks=synthetic_exterior,
        actual_outside_state_checks=actual_outer,
        actual_soft_stiff_prefix_bitwise_equal_before_first_exit=True,
        actual_soft_stiff_post_exit_max_radius_difference_A=float(np.max(abs(soft[prefix:,4]-stiff_trace[prefix:,4]))),
        no_claim_of_timestep_convergence=True)
    data['verification'] = result
    old = json.loads((HERE/'results.json').read_text())
    anchors = [r for r in old['rows'] if r['scope']=='target' and
               r['energy_per_D_eV']==.0006 and r['pattern'] in ('converging_wave','direct_optimized')]
    anchors += data['cases']['size_N16384']['selected']+data['cases']['size_N131072']['selected']
    scaling = {}
    for pattern in ('converging_wave','direct_optimized'):
        rows = sorted((r for r in anchors if r['pattern']==pattern),key=lambda r:r['n_D'])
        fields = ('n_D','peak_relative_kinetic_eV','pair_gain','fraction_total_energy_at_peak',
                  'peak_time_fs','half_peak_hold_fs','signed_inward_speed_A_per_fs_at_peak',
                  'first_site_or_bond_exit_fs','max_site_radius_A_before_exit','min_bond_A_before_exit',
                  'participation_at_launch','target_kinetic_preload_fraction')
        scaling[pattern] = [dict((k,r[k]) for k in fields) for r in rows]
        for first,second in zip(scaling[pattern],scaling[pattern][1:]):
            second['adjacent_finite_peak_exponent'] = float(
                np.log(second['peak_relative_kinetic_eV']/first['peak_relative_kinetic_eV'])/
                np.log(second['n_D']/first['n_D']))
    data['finite_size_scaling'] = scaling
    data['resources']['N131072'].update(
        approved_parallel_numerical_workers=2,
        actual_full_search_seconds=data['cases']['size_N131072']['runtime_seconds'],
        decision='full unchanged search carried by two one-thread workers after measured supported probe',
        harmonic_interpretation=dict(source='parent independent sampled 64^3 exact-symbol grid',
            sampled_max_group_speed_A_fs=.0386909,sampled_max_110_projection_A_fs=.0339444,
            sampled_fastest_2400fs_distance_A=92.86,
            caveat='sampled harmonic speed is not a nonlinear causal bound; no horizon sweep'))
    data['scaling_interpretation'] = 'old anchors reused without dynamics; finite capped launch classes, finite2400fs arrival/horizon and original support constraints do not prove saturation'
    save(data,traces)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    verify()
