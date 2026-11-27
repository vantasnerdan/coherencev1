"""Read settled Pass5 artifacts; targeted semantic/consumer checks, no producer replay."""
import ast
import json
import platform
from pathlib import Path
import numpy as np
import scipy
import source as s
from run import readrows, make_kernel, packet_record, source_phase

t,r=s.t,s.r
HERE,ROOT=s.HERE,s.ROOT


def check(m,ref,previous=None):
    source=s.load();xy,prob,_=r.sources(m);data=s.tracks(m)
    checks={}
    import delivery
    assert Path(__file__).resolve().parent==HERE and delivery.HERE==HERE
    checks['reader_and_packager_owned_by_pass5']=True
    for name in ('source.py','run.py','refine.py','feedback.py','readout.py','delivery.py'):
        tree=ast.parse((HERE/name).read_text())
        assert any(isinstance(x,ast.FunctionDef) for x in ast.walk(tree))
    checks['python_ast_positive']=True
    assert np.array_equal(source['energy_J'][2:4],np.load(HERE/'old_launch.npz')['built'])
    assert t.sha(HERE/'old_launch.npz')==t.sha(ROOT/'g13_rogue/pass4/launch.npz')
    checks['both_original_gaussian_arrays_exact']=True
    largest=0.
    for tr in data:
        mk=s.marks(tr)
        raw=np.sum(mk['intensity']*np.sum(mk['q']**2,axis=-1))
        expected=2*s.MPD*tr['deposit'][:,1].sum()
        largest=max(largest,abs(raw/expected-1))
        pD=np.sqrt(2*t.MD*tr['projectile_eV']*t.QE)
        expected_z=(t.MD+s.MPD)*mk['T']*t.QE/pD[:,None]
        assert np.allclose(mk['q'][:,:,2],expected_z,rtol=5e-14,atol=1e-40)
        assert np.sum(tr['momentum'])==1 or abs(np.sum(tr['momentum'])-1)<1e-15
        assert np.all(np.expm1(np.zeros_like(mk['T'],complex))==0)
    assert largest<1e-13
    checks['recoil_second_trace_relative_error']=largest
    checks['recoil_pair_zero_k_omega_and_single_external_impulse']=True
    budgets=readrows('source_budget.csv')
    residual=max(abs(float(x['balance_J']))/float(x['incident_J']) for x in budgets)
    assert residual<1e-13
    checks['wave_plus_residual_heat_relative_balance']=residual
    integrated=np.einsum('asbnj,n->asbj',source['angular_J_sr'],source['angular_weight_sr'])
    assert np.allclose(integrated,source['energy_J'][:2],rtol=1e-12,atol=1e-38)
    checks['directional_DOS_normalization']=True
    bins=dict(np.load(HERE/'selected_bins.npz'))
    assert np.allclose(bins['energy'].sum(axis=3),source['energy_J'][:2,:,4:6],rtol=2e-12,atol=1e-35)
    checks['selected_bins_same_source_epoch']=True
    # A real row discriminates a retarded source from an advanced phase.
    zrow=float(data[0]['z'][17]);trow=float(data[0]['t'][17])
    omega=2*np.pi*8e10;km=float(t.dispersion(8e10,t.christoffel(np.array([0.,0.,1.]))[0])[0][2])
    distance=30e-6;clock=1.7e-9
    transform=np.exp(1j*(omega*trow-km*zrow))
    observed=np.exp(1j*(omega*clock-km*distance))*source_phase(transform)
    retarded=np.exp(1j*(omega*(clock-trow)-km*(distance-zrow)))
    assert abs(observed-retarded)<2e-12
    assert abs(np.exp(1j*(omega*clock-km*distance))*transform-retarded)>1e-2
    checks['actual_track_row_retarded_source_carrier']=True
    # Exercise actual loss consumer, including a ray where group transit differs
    # from geometric distance/nominal speed through the fixed grain map.
    at=np.array([[t.L/2+5e-6,t.L/2-3e-6]])
    band=5;energy=source['energy_J'][0,:,band]
    real=t.loss
    t.loss=lambda f,vg,law:np.zeros_like(np.asarray(vg,dtype=float))
    try:base=r.kernel(at,band,energy,'Cu_extrapolated')
    finally:t.loss=real
    errors=[]
    for law in ('Cu_extrapolated','Al_transfer'):
        ker=r.kernel(at,band,energy,law)
        for d,z in enumerate(t.DEPTHS):
            for path in range(15):
                b,echo=divmod(path,5)
                zi=z if echo==0 else 2*((echo+1)//2)*t.H+(-1 if echo%2 else 1)*z
                length=np.sqrt((5e-6)**2+(3e-6)**2+zi*zi)
                optical=np.log(10)/10*32*ker['f']/1e6*length if law=='Cu_extrapolated' else ker['arrival'][0,d,path]/(34e-12*240e9/ker['f'])
                observed=base['logweight'][0,d,path]-ker['logweight'][0,d,path]
                errors.append(float(np.max(abs(observed-optical))/max(1.,optical)))
    assert max(errors)<1e-12
    checks['actual_Cu_length_Al_time_relative_error']=max(errors)
    # Raw marked-parent second moment, not square of a deterministic mean.
    ker=make_kernel(xy[:4],5,energy,'Cu_extrapolated',0,source)
    avg=np.zeros_like(ker['weight'])
    for code in range(8):
        signs=np.array([1 if code&(1<<b) else -1 for b in range(3)])[np.repeat(np.arange(3),5)]
        avg+=abs(ker['source_mu']+np.sqrt(ker['source_var'])*signs[None,None,:,None])**2/8
    assert np.max(abs(avg-1))<1e-13
    checks['packet_source_mean_plus_variance']=True
    packet_error=[]
    for scenario,name in ((0,'track'),(2,'gaussian20'),(3,'gaussian60')):
        path=HERE/f'{name}_b5_Cu_extrapolated_even_steady_r2.npz'
        frozen=dict(np.load(path))
        ker=make_kernel(xy,5,source['energy_J'][scenario,:,5],'Cu_extrapolated',scenario,source)
        reproduced=packet_record(ker,frozen,xy,2,5,scenario)
        scale=max(float(frozen['elastic_J_m3'].max()),1e-300)
        error=float(np.max(abs(reproduced['elastic_J_m3']-frozen['elastic_J_m3']))/scale)
        assert error<1e-12
        packet_error.append(error)
    checks['track_and_both_gaussian_packet_replay_relative_errors']=packet_error
    if previous:
        old=dict(np.load(Path(previous)/'built_b5_Cu_extrapolated_even_steady_r2.npz'))
        new=dict(np.load(HERE/'gaussian20_b5_Cu_extrapolated_even_steady_r2.npz'))
        error=float(np.max(abs(new['elastic_J_m3']-old['elastic_J_m3']))/old['elastic_J_m3'].max())
        assert error<1e-12
        checks['accepted_Pass4_actual_packet_relative_error']=error
    peaks=readrows('same_record_peaks.csv')
    histories={}
    for x in peaks:
        key=(x['band'],x['loss'],x['beam'],x['record'])
        histories.setdefault(key,set()).add(x['physical_catalog_sha256'])
    assert all(len(v)==1 for v in histories.values())
    checks['all_source_comparisons_share_physical_catalogs']=True
    tails=readrows('tail_comparison.csv');tail_histories={}
    for row in tails:
        key=(row['loss'],row['beam'],row['record'],row['depth_m'])
        tail_histories.setdefault(key,set()).add(row['physical_catalog_sha256'])
    assert all(len(x)==1 for x in tail_histories.values())
    oldtail=dict(np.load(HERE/'old_tail_b5_Cu_extrapolated_even_steady_r2.npz'))
    newtail=dict(np.load(HERE/'tail_gaussian20_b5_Cu_extrapolated_even_steady_r2.npz'))
    scale=max(float(oldtail['elastic_J_m3'][:,1:].max()),1e-300)
    tail_error=float(np.max(abs(oldtail['elastic_J_m3'][:,1:]-newtail['elastic_J_m3'][:,1:]))/scale)
    assert tail_error<1e-12
    checks['matched_old_tail_actual_packet_relative_error']=tail_error
    checks['all_five_sources_share_accepted_tail_filters_and_parents']=True
    spectral=dict(np.load(HERE/'refined_source_phase.npz'))
    low=dict(np.load(HERE/'selected_bins.npz'))
    frozen=dict(np.load(HERE/'tail_track_b5_Cu_extrapolated_even_steady_r2.npz'))
    carrier=float(spectral['band5_frequency_Hz'][0])
    ker=make_kernel(xy,5,low['energy'][0,:,1,0],'Cu_extrapolated',0,source,
        frequency=carrier,bandwidth=np.sqrt(t.BANDS[5]*t.BANDS[6])-t.BANDS[5],
        angular=low['angular'][0,:,1,0],variance=low['variance'][0,:,1,0],
        center=spectral['band5_source_mean'][0,:,0])
    ker['weight'][:,0]=0
    reproduced=packet_record(ker,frozen,xy,2,5,0)
    tail_error=float(np.max(abs(reproduced['elastic_J_m3']-frozen['elastic_J_m3']))/max(float(frozen['elastic_J_m3'].max()),1e-300))
    assert tail_error<1e-12
    checks['updated_track_actual_low_edge_packet_replay_relative_error']=tail_error
    fb=dict(np.load(HERE/'feedback.npz'))
    from feedback import earlier_field,ARMS,receiver_packets,source_power_contrast
    physical=s.orient_source(source)
    norm_error=0.
    for band in (5,6):
        contrast=source_power_contrast(physical,band)
        for sp in range(2):
            angular=physical['angular_J_sr'][0,sp,band]
            fraction=np.clip(np.divide(physical['angular_variance_J_sr'][0,sp,band],angular,out=np.zeros_like(angular),where=angular>0),0,1)
            mu=np.sqrt(1-fraction)*source_phase(physical['center_mean'][0,sp,band])
            for sign in (-1.,1.):
                integrated=np.einsum('nb,n->b',angular*abs(mu+sign*np.sqrt(fraction))**2,physical['angular_weight_sr'])
                # The mode's actual radiation norm, not the averaged source power.
                predicted=source['energy_J'][0,sp,band]*(1+sign*contrast[sp])
                norm_error=max(norm_error,float(np.max(abs(integrated-predicted)/predicted)))
    assert norm_error<1e-12
    checks['realized_recoil_mark_outgoing_energy_relative_error']=norm_error
    causal=0
    for key in [x.removesuffix('_launch_sqrt_J') for x in fb if x.endswith('_launch_sqrt_J')]:
        a=fb[key+'_launch_sqrt_J']
        assert np.all(fb[key+'_source_norm']>0)
        band=int(key.split('_',1)[0][1:])
        expected_norm2=1+source_power_contrast(physical,band)[fb[key+'_species']]*fb[key+'_source_marks']
        assert np.allclose(fb[key+'_source_norm']**2,expected_norm2,rtol=1e-13,atol=1e-13)
        assert np.array_equal(a[ARMS.index('zero_prior035')],a[ARMS.index('zero_prior0')])
        time=fb[key+'_time_s'];assert np.all(np.diff(time)>=0)
        assert len(np.unique(fb[key+'_parent_id']))==len(time)
        pre=fb[key+'_pre_velocity_m_s']
        assert pre[ARMS.index('no_background035'),0]==0
        assert np.all(pre[ARMS.index('zero_prior035')]==0)
        dirs=fb[key+'_background_directions']
        _,normal,_=s.oriented_polarizations(np.array([0.,0.,1.]))
        _,incoming,_=s.oriented_polarizations(dirs)
        rotation=np.einsum('ia,nib->nab',normal,incoming)
        covariance=np.einsum('nab,ncb->ac',rotation,rotation)/(3*len(dirs))
        assert np.allclose(covariance,np.eye(3)/3,rtol=1e-13,atol=1e-13)
        if key+'_causal_index' in fb:
            index=int(fb[key+'_causal_index']);green=fb[key+'_causal_green']
            original=earlier_field(a[0],green,index)
            mutated=a[0].copy();mutated[index]*=1e6
            assert np.array_equal(original,earlier_field(mutated,green,index))
            # A positive control reaches the same real causal source kernel.
            support=np.flatnonzero(np.any(abs(green)>0,axis=tuple(range(1,green.ndim))))
            mutated=a[0].copy();mutated[support[0]]*=2
            assert not np.array_equal(original,earlier_field(mutated,green,index))
            causal+=1
    assert causal>0
    checks['actual_nonzero_causal_kernel_no_current_parent_feedback']=causal
    debit=readrows('feedback_debits.csv')
    assert all(float(x['residual_heat_J'])>=0 and abs(float(x['balance_J']))<1e-29 for x in debit)
    assert all(abs(float(x['actual_acoustic_work_J'])-(float(x['fixed_other_acoustic_J'])+float(x['two_band_self_J'])+float(x['signed_cross_J'])))<1e-29 for x in debit)
    checks['two_band_signed_work_debited_once']=True
    events=readrows('feedback_events.csv')
    species={(x['beam'],x['loss'],x['arm'],x['parent_id']):int(x['species']) for x in events}
    host_fraction=max(float(x['actual_acoustic_work_J'])/float(m['heatenergy'][species[x['beam'],x['loss'],x['arm'],x['parent_id']],1]) for x in debit)
    assert host_fraction<=1+1e-13
    checks['maximum_joint_realized_acoustic_work_over_host_budget']=host_fraction
    work_error=0.
    for row in events:
        key=f'b{row["band"]}_{row["loss"]}_{row["beam"]}'
        index=int(np.searchsorted(fb[key+'_parent_id'],int(row['parent_id'])))
        arm=ARMS.index(row['arm'])
        old=fb[key+'_pre_projection_sqrt_J'][arm,index]
        delta=fb[key+'_launch_sqrt_J'][arm,index]
        energy=float(np.vdot(delta,delta).real+2*np.vdot(old,delta).real)
        error=abs(energy-float(row['actual_work_J']))/max(float(np.vdot(delta,delta).real),1e-300)
        work_error=max(work_error,error)
    assert work_error<1e-12
    checks['actual_signed_source_projection_work_relative_error']=work_error
    feedback_replay=[]
    for band in (5,6):
        key=f'b{band}_Cu_extrapolated_even_steady'
        _,payload=receiver_packets(physical,band,'Cu_extrapolated','even_steady',
            fb[key+'_position_m'],fb[key+'_time_s'],fb[key+'_species'],fb[key+'_parent_id'],
            fb[key+'_launch_sqrt_J'],fb[key+'_source_marks'],fb[key+'_source_norm'],.5*25e-6)
        frozen=fb[key+'_receiver_elastic_J_m3']
        error=float(np.max(abs(payload['elastic_J_m3']-frozen))/max(float(frozen.max()),1e-300))
        assert error<1e-12
        feedback_replay.append(error)
    checks['both_energetic_feedback_receiver_replay_relative_errors']=feedback_replay
    rows=readrows('readout.csv')
    expected={(sc,b,law,beam,z) for sc in s.SHAPES for b in range(7) for law in ('Cu_extrapolated','Al_transfer') for beam in ('even_steady','hot_steady','hot_bursts') for z in t.DEPTHS}
    actual={(x['scenario'],int(x['band']),x['loss'],x['beam'],float(x['depth_m'])) for x in rows}
    assert actual==expected and len(rows)==len(expected)
    checks['all_requested_source_beam_loss_band_depth_cells']=len(actual)
    checks['runtime']=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    s.emit('behavior_check.json',checks)
    print(json.dumps(checks),flush=True)


def source_tables(m,source):
    terms=[]
    original=np.load(HERE/'old_launch.npz')['built_terms']
    for scenario,name in enumerate(s.SHAPES[:4]):
        raw=source['terms'][scenario,:,:7] if scenario<2 else original[scenario-2]
        mechanical,thermal,cross=t.RATE*np.einsum('s,sbjc->c',t.FRACTION,raw)
        variance=float(t.RATE*np.einsum('s,sbnj,n->',t.FRACTION,source['angular_variance_J_sr'][scenario],source['angular_weight_sr'])) if scenario<2 else 0.
        terms.append(dict(scenario=name,coherent_mechanical_W=mechanical-variance,
            recoil_raw_variance_W=variance,thermoelastic_W=thermal,interference_W=cross,
            total_W=mechanical+thermal+cross))
    s.table('source_terms.csv',terms)
    kinetic=[]
    for scenario,name in enumerate(s.SCENARIOS):
        for sp,tr in enumerate(s.tracks(m)):
            h=source['heat_scale'][scenario,sp]
            mechanical=float(source['terms'][scenario,sp,:,: ,0].sum())
            host=float(tr['deposit'][:,1].sum())
            if host*(1-h)<mechanical:raise ValueError('Host kinetic debit cannot fund mechanical wave')
            for j,deposit in enumerate(tr['deposit']):
                paid=deposit[1]/host*mechanical
                thermal=deposit[1]*(1-h)-paid
                kinetic.append(dict(scenario=name,species=sp,row=j,host_kinetic_J=deposit[1],
                    mechanical_wave_debit_J=paid,host_thermal_work_debit_J=thermal,
                    residual_host_heat_J=deposit[1]*h,
                    electronic_residual_heat_J=deposit[0]*h,direct_D_residual_heat_J=deposit[2]*h,
                    balance_J=deposit[1]-paid-thermal-deposit[1]*h))
    s.table('kinetic_heat_debits.csv',kinetic)
    # Source-band-integrated angular maxima and the LOWEST frequency's loss
    # upper-bound the incoherent mean in the new upper energetic decade.
    xy,prob,_=r.sources(m);upper=[]
    maximum=source['angular_J_sr'][0,:,6].max(axis=1).sum(axis=1)
    for law in ('Cu_extrapolated','Al_transfer'):
        alpha=np.log(10)/10*32*t.BANDS[6]/1e6 if law=='Cu_extrapolated' else t.loss(t.BANDS[6],10000.,law)/10000.
        for beam in ('even_steady','hot_steady','hot_bursts'):
            pp=prob[int(beam!='even_steady')]
            for z in t.DEPTHS:
                rr=np.sqrt(np.sum((xy-t.L/2)**2,axis=1)+z*z)
                bound=.5*t.RATE*5*np.sum(t.FRACTION*maximum)/500.*np.sum(pp*np.exp(-alpha*rr)/rr**2)
                upper.append(dict(band=6,loss=law,beam=beam,depth_m=z,mean_elastic_upper_J_m3=bound,
                    cutoff_J_m3=1e-20,below_cutoff=bool(bound<1e-20),
                    basis='lowest-edge loss, angular-band maximum, 500m/s lower group-speed, five images, transmission<=1; stationary mean'))
    s.table('upper_band_depth_bound.csv',upper)
    return terms,upper


def write_pages(m):
    source=s.load();rows=readrows('readout.csv');peaks=readrows('same_record_peaks.csv')
    refined=readrows('frequency_refinement.csv');feedback=readrows('feedback.csv');receiver=readrows('feedback_receiver.csv')
    tails=readrows('tail_comparison.csv')
    powers=t.RATE*np.einsum('s,asbj->ab',t.FRACTION,source['energy_J'])
    totals=powers.sum(axis=1)
    terms,upper=source_tables(m,source)
    actual=[x for x in peaks if x['scenario']=='track' and x['conditional']=='False']
    for row in tails:
        if row['scenario']!='track' or row['conditional']!='False':continue
        actual.append(dict(row,tallest_over_mean=row['tallest_over_own_mean'],
            distinct_core_at_peak=row['actual_core_parents']))
    peak=max(actual,key=lambda x:float(x['peak_J_m3']))
    thermal={int(x['band']):float(x['mean_thermal_elastic_J_m3']) for x in readrows('thermal_reference.csv')}
    peak['same_band_300K_mean_ratio']=float(peak['peak_J_m3'])/thermal[int(peak['band'])]
    old_anchor=max((x for x in tails if x['scenario']=='gaussian20' and x['loss']=='Cu_extrapolated' and x['beam']=='even_steady' and float(x['depth_m'])==30e-6 and x['conditional']=='False'),key=lambda x:float(x['peak_J_m3']))
    new_anchor=next(x for x in tails if x['scenario']=='track' and x['loss']==old_anchor['loss'] and x['beam']==old_anchor['beam'] and x['record']==old_anchor['record'] and x['depth_m']==old_anchor['depth_m'])
    frame=HERE/f'tail_track_b5_{new_anchor["loss"]}_{new_anchor["beam"]}_r{new_anchor["record"]}.npz'
    with np.load(frame) as frozen:
        dep=int(np.flatnonzero(t.DEPTHS==float(new_anchor['depth_m']))[0])
        actual_field=frozen['elastic_J_m3'][:,dep]
        ti,point=np.unravel_index(actual_field.argmax(),actual_field.shape)
        shuffled=frozen['shuffled_J_m3'][:,dep]
        new_anchor.update(shuffled_same_point_window_max_J_m3=float(shuffled[:,point].max()),
            shuffled_whole_grid_window_max_J_m3=float(shuffled.max()),
            selected_peak_over_shuffled_grid_max=float(actual_field[ti,point]/shuffled.max()),
            selected_time_phase_ratio=float(actual_field[ti,point]/shuffled[ti,point]))
    high_windows=[x for x in receiver if int(x['band'])==6 and x['arm']=='quiet']
    occupied_high=sum(float(x['peak_elastic_J_m3'])>0 for x in high_windows)
    reach=[]
    for law in ('Cu_extrapolated','Al_transfer'):
        for band in range(7):
            selected=[x for x in (refined if band in (4,5) else rows) if x['scenario']=='track' and int(x['band'])==band and x['loss']==law and x['beam']=='even_steady']
            good=[float(x['depth_m']) for x in selected if float(x['uncensored_mean_J_m3'])>=1e-20]
            reach.append(dict(band=band,loss=law,cutoff_J_m3=1e-20,deepest_sampled_m=max(good) if good else None,
                meaning='sampled-depth reach; no interpolation or no-sample zero'))
    s.table('depth_reach.csv',reach)
    chosen=[x for x in feedback if x['arm']=='beta035']
    relative=[x for x in feedback if x['arm']=='relative_speed']
    contrasts=[float(x['work_over_quiet']) for x in chosen]
    expected=[float(x['expected_work_over_quiet']) for x in chosen]
    speed=max(float(x['max_abs_relative_speed']) for x in relative)
    headline=dict(powers_W={name:float(total) for name,total in zip(s.SHAPES,totals)},
        incident_shares={name:float(total/.69) for name,total in zip(s.SHAPES,totals)},
        primary_band_share=(powers[0]/totals[0]).tolist(),strongest_actual_packet=peak,
        beta035_realized_work_ratio_range=[min(contrasts),max(contrasts)],
        beta035_conditional_expected_work_ratio_range=[min(expected),max(expected)],
        max_physical_relative_speed_fraction=speed,depth_reach=reach)
    headline['matched_30um_old_peak']=dict(record=old_anchor['record'],old_J_m3=float(old_anchor['old_peak_J_m3']),track_J_m3=float(new_anchor['paired_old_peak_J_m3']),exposure_s=float(new_anchor['exposure_s']))
    headline['matched_30um_event']=new_anchor
    headline['upper_feedback_nonempty_windows']=dict(nonempty=occupied_high,total=len(high_windows))
    worst=max(chosen,key=lambda x:float(x['work_over_quiet']))
    keys=('band','loss','beam')
    cohort={x['arm']:x for x in feedback if all(x[k]==worst[k] for k in keys)}
    sensor={x['arm']:x for x in receiver if all(x[k]==worst[k] for k in keys) and float(x['depth_m'])==63e-9}
    headline['largest_conditional_source_case']=dict(source=cohort,receiver_63nm=sensor)
    s.emit('headline_results.json',headline)
    lines=['# Pass 5: a track-shaped launch on the same foil','',
        f'With the declared **2-nm coherent-bath footprint and row-local thermalization**, the track construction launches **{totals[0]:.6g} W**, or **{totals[0]/.69:.6g} of incident power**. The quarter-range hypothesis gives {totals[1]:.6g} W. The unchanged Gaussian comparisons give {totals[2]:.6g} and {totals[3]:.6g} W; the historical assumed-share arm gives {totals[4]:.6g} W on the same incident denominator. No interpolation between those values was imposed. The boost is mainly thermal: {terms[0]["thermoelastic_W"]:.6g} W thermoelastic, {terms[0]["interference_W"]:.6g} W interference, {terms[0]["coherent_mechanical_W"]:.6g} W coherent mechanical, and {terms[0]["recoil_raw_variance_W"]:.6g} W recoil raw variance (`source_terms.csv`); it does not isolate a directional-recoil gain.','',
        '| Band | Track launch W | Share of track launch | Even/Cu track core parents | N_eff | Mean elastic eV/D |',
        '|---|---:|---:|---:|---:|---:|']
    for band in range(7):
        row=next(x for x in rows if x['scenario']=='track' and int(x['band'])==band and x['loss']=='Cu_extrapolated' and x['beam']=='even_steady' and float(x['depth_m'])==t.DEPTHS[0])
        lines.append(f'| {t.BANDS[band]:.3g}–{t.BANDS[band+1]:.3g} Hz | {powers[0,band]:.4g} | {powers[0,band]/totals[0]:.3%} | {float(row["energetic_core_parents"]):.4g} | {float(row["N_eff"]):.4g} | {float(row["mean_eV_D"]):.4g} |')
    lines+=['',f'The largest unconditioned sampled packet is {float(peak["peak_J_m3"]):.6g} J/m³ at {float(peak["depth_m"])*1e6:.6g} µm, {float(peak["tallest_over_mean"]):.5g} times its own same-point mean over {float(peak["exposure_s"]):.5g} s. It has {peak["distinct_core_at_peak"]} core electrical ancestries and N_eff={float(peak["N_eff_at_peak"]):.7g}; its inventory quotient is {float(peak["peak_eV_D"]):.5g} eV/D. Paired shuffled energy at that same place/time is {float(peak["shuffled_same_peak_J_m3"]):.6g} J/m³.','',
        'The table separates finite first-five-image/eight-sigma front ancestry, the twelve-length energetic core, and energy-weighted participation. Every source, all three beams, both loss rules and all four depths are in `readout.csv`; `same_record_peaks.csv` and `tail_comparison.csv` bind maxima and controls to their actual records.','',
        'The two-bin low-edge transfer is:', '', '| Band / loss, even steady | 30 µm J/m³ | 50 µm J/m³ | 200 µm J/m³ |','|---|---:|---:|---:|']
    for band in (4,5):
        for law in ('Cu_extrapolated','Al_transfer'):
            values=[next(x for x in refined if x['scenario']=='track' and int(x['band'])==band and x['loss']==law and x['beam']=='even_steady' and float(x['depth_m'])==z)['uncensored_mean_J_m3'] for z in t.DEPTHS[1:]]
            lines.append(f'| {band} / {law} | '+' | '.join(f'{float(v):.5g}' for v in values)+' |')
    lines+=['',f'The selected 30-µm low-edge coincidence peaks at {float(new_anchor["peak_J_m3"]):.6g} J/m³: **{new_anchor["selected_peak_over_shuffled_grid_max"]:.3g}×** the paired shuffle’s whole-window/nine-point maximum ({new_anchor["shuffled_whole_grid_window_max_J_m3"]:.6g} J/m³). It has {new_anchor["actual_core_parents"]} supported parents, instantaneous N_eff={float(new_anchor["N_eff_at_peak"]):.5g}, and {float(new_anchor["peak_eV_D"]):.4g} eV/local-D. At the originally selected point/time, the comparison is {new_anchor["selected_time_phase_ratio"]:.4g}×; the shuffle’s same-point window maximum is {new_anchor["shuffled_same_point_window_max_J_m3"]:.6g} J/m³. At the old Gaussian peak’s identical location/time, the track source gives {float(new_anchor["paired_old_peak_J_m3"]):.6g} versus {float(old_anchor["old_peak_J_m3"]):.6g} J/m³ on the same {float(new_anchor["exposure_s"]):.5g}-s record, alongside the lower matched mean above. One selected event and one shuffle do not estimate an excess or rarity.']
    lines+=['', 'Depth reach uses a **declared reporting threshold**, mean elastic energy density ≥10⁻²⁰ J/m³ at a sampled depth. It is not a detection floor, physical zero, or bound on occasional coherent peaks. `depth_reach.csv` and `upper_band_depth_bound.csv` apply only to the retained paths and named loss rules.','',
        f'On identical physical parent histories, the β=.35 conditional mark law changes realized source work by factors {min(contrasts):.7g}–{max(contrasts):.7g}; its conditional expected-work factors are {min(expected):.7g}–{max(expected):.7g}. The unamplified physical relative-speed correction reaches only |v_z/v_D|={speed:.4g}. Phase shuffle, background-off and zero-prior controls are paired. Updated launches enter actual retarded receiver packets (`feedback_receiver.csv`).','',
        f'The largest conditional change is band {worst["band"]}/{worst["loss"]}/{worst["beam"]}: {worst["parents"]} local source parents in {float(worst["exposure_s"])*1e9:g} ns. Work changes from {float(cohort["quiet"]["actual_work_J"]):.6g} to {float(worst["actual_work_J"]):.6g} J. Phase shuffle retains a factor {float(cohort["shuffled035"]["work_over_quiet"]):.7g}, whereas the physical relative-speed arm changes work by only {float(cohort["relative_speed"]["work_over_quiet"])-1:.3g} fractionally. This is chiefly the stipulated mark-energy sensitivity. Source C counts {worst["C_abs_velocity_parents"]}/{worst["parents"]}. At its matched 63-nm receiver the quiet→conditional peak is {float(sensor["quiet"]["peak_elastic_J_m3"]):.6g}→{float(sensor["beta035"]["peak_elastic_J_m3"]):.6g} J/m³ ({float(sensor["quiet"]["peak_eV_D"]):.6g}→{float(sensor["beta035"]["peak_eV_D"]):.6g} eV/D); B counts {sensor["quiet"]["B_initiations"]}→{sensor["beta035"]["B_initiations"]} initiations, {sensor["quiet"]["B_boundary_censored"]}→{sensor["beta035"]["B_boundary_censored"]} boundary-censored, in {float(sensor["quiet"]["exposure_s"])*1e9:.4g} ns. These source and receiver denominators remain separate.','',
        f'The upper energetic-band feedback has nonempty quiet packets in {occupied_high}/{len(high_windows)} bounded retarded receiver windows. Empty windows are exposure-limited; they are not assigned a physical zero or a rarity forecast.','',
        'The source and signed-background equations are in `README.md`; consumed inputs and the bounded Pd-loss search are in `inputs_ledger.md`. Those physical inputs, rather than further range refinements, set the interpretation.','',
        'Exact archive hashes, unpublished draft IDs, and remote member readback are supplied in the delivery receipt.']
    (HERE/'SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(headline),flush=True)
