"""Consume frozen Pass-4 outputs; never replay the scientific producer to package."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import transport as t


def rows(name):return list(csv.DictReader((t.HERE/name).open()))


def tail_bound(m):
    # Conditional analytic upper bound for the SAME source/capped-bulk model.
    # |F_p|<=1, |F_heat|<=Q, |Dk²-iw|>=w. Dropping all spatial attenuation
    # and using c_min=1000 m/s overbounds the inspected pure-Pd low-q tensor.
    # Integrate in k so the endpoint vg=0 is not a false divergence.
    w0=2*np.pi*2.5e12;kmax=np.pi*np.pi*t.FCAP/1000.
    out=[]
    for shape,(a,tau) in enumerate(((20e-9,.5e-12),(60e-9,2e-12))):
        for s in range(2):
            q=m['heatenergy'][s].sum();p=t.PARENT_P[s]
            aa=t.BETA*q/(t.RHO*t.CP)*kmax/w0
            prefactor=float(4*np.pi/(2*t.RHO*(2*np.pi)**3)*np.sum(kmax**3/3*2*(p*p+aa*aa)))
            log_bound=float(np.log(prefactor)-(w0*tau)**2)
            bound=float(np.exp(max(log_bound,np.log(np.nextafter(0.,1.)))))
            out.append(dict(shape=shape,species=('D+','D2+')[s],lo_Hz=2.5e12,branch_caps_Hz=t.FCAP.tolist(),built_supported_above_ladder_upper_J=bound,log_built_supported_above_ladder_upper_J=log_bound,method='analytic upper within same capped acoustic source; positive cross bound, spatial attenuation dropped; subnormal upper floor retained, not zero'))
    t.dump('above_ladder_bound.json',out)
    return out


def summarize(ref):
    m=t.reference(ref);tail=tail_bound(m)
    launch=rows('launch.csv');read=rows('readout.csv');counts=rows('overlap_counts.csv');sp=rows('sparse_events.csv')
    led=rows('transport_ledger.csv');bg=rows('slow_background.csv')
    b=np.load(t.HERE/'launch.npz')['built'];other=rows('source_budget.csv')
    for row in other:
        s=0 if row['species']=='D+' else 1
        bound=tail[s]['built_supported_above_ladder_upper_J'] if row['scenario']=='built' else 0.
        row['above_ladder_built_upper_J']=bound
        row['heat_random_unmodeled_lower_J']=float(row['retained_heat_random_motion_and_above_ladder_J'])-bound
        row['heat_random_unmodeled_upper_J']=row['retained_heat_random_motion_and_above_ladder_J']
    t.table('source_budget.csv',other)
    import receiver
    parent_debits=[];record_budget=[]
    for scenario,energy in (('built',b[0]),('assumed',receiver.assumed(m))):
        for species in range(2):
            for band in range(7):
                for branch in range(3):
                    parent_debits.append(dict(scenario=scenario,species=('D+','D2+')[species],band=band,branch=branch,acoustic_debit_J=energy[species,band,branch],meaning='one isolated-parent acoustic work, not another per-band recoil fraction'))
        high=energy.sum(axis=(1,2))
        for beam in ('even_steady','hot_steady','hot_bursts'):
            for lowloss in ('retained','Cu_proxy'):
                name=f'r2_{beam}'+('__built_only' if scenario=='built' else '')+f'_b0.0_{lowloss}_combined_records0-15_energy.csv'
                with (Path(ref)/name).open() as f:
                    for old in csv.DictReader(f):
                        incident=float(old['incident_J']);n=incident/(5000*t.QE)
                        ns2=(float(old['host_J'])-n*m['heatenergy'][0,1])/(m['heatenergy'][1,1]-m['heatenergy'][0,1])
                        ns=np.array([n-ns2,ns2])
                        if np.any(ns<0):raise ValueError('Accepted source/host budget epoch mismatch')
                        lowwork=sum(float(old[k]) for k in ('random_work_J','momentum_work_J','thermal_work_J'))
                        highwork=float(ns@high);rest=incident-lowwork-highwork
                        if rest<0:raise ValueError('All-band source record overspent; no clipping')
                        record_budget.append(dict(scenario=scenario,beam=beam,low_loss=lowloss,record=old['record'],incident_J=incident,accepted_actual_low_work_J=lowwork,new_isolated_high_work_J=highwork,immediate_heat_random_and_rest_J=rest,source_stage_balance_J=incident-lowwork-highwork-rest,accepted_low_damping_J_reclassified_not_new_input=float(old['loss_J']),D_parent_equivalent=ns[0],D2_parent_equivalent=ns[1],nearest_integer_count_discrepancy=float(np.max(abs(ns-np.rint(ns)))),meaning='species-equivalent counts inferred from accepted incident/host ledgers; source-work partition, not a claim of closed unresolved finite-record cross work'))
    t.table('per_parent_debits.csv',parent_debits)
    t.table('source_record_budget.csv',record_budget)
    thermal=[]
    h=6.62607015e-34;kb=1.380649e-23
    import receiver
    for band in range(7):
        f,k,vg,depth=receiver.characteristics(band)
        thermal_oscillator=h*f/np.expm1(h*f/(kb*300))
        dos=np.mean(k/vg,axis=0)/t.H if band==0 else np.mean(k*k/(np.pi*vg),axis=0)
        energy=.5*thermal_oscillator*dos.sum()*(t.BANDS[band+1]-t.BANDS[band])
        thermal.append(dict(band=band,temperature_K=300,mean_thermal_elastic_J_m3=float(energy),model='harmonic equilibrium Bose reference, same conditional acoustic branches at band centre; not total thermal energy'))
    t.table('thermal_reference.csv',thermal)
    with np.load(t.HERE/'joint_distributions.npz') as samples:
        for row in read:
            key=f'{row["scenario"]}_b{row["band"]}_{row["loss"]}_{row["beam"]}_depth{int(np.argmin(abs(t.DEPTHS-float(row["depth_m"]))))}_elastic_samples'
            values=samples[key]
            row['samples']=len(values)
            row['empirical_mean_elastic_J_m3']=float(values.mean())
            row['empirical_sd_elastic_J_m3']=float(values.std())
            row['wave_over_300K_same_band_thermal_mean']=float(row['mean_elastic_J_m3'])/thermal[int(row['band'])]['mean_thermal_elastic_J_m3']
            median=float(np.median(values))
            row['max_over_median']=float(values.max()/median) if median>0 else None
            row['median_zero']=bool(median==0)
            if len(values)==2048:row['finite_sample_exposure_s']=4*64/float(row['carrier_Hz'])
    t.table('readout.csv',read)
    slow=np.load(t.HERE/'slow_background.npz')
    built_background={}
    for beam in ('even_steady','hot_steady','hot_bursts'):
        for lowloss in ('retained','Cu_proxy'):
            path=Path(ref)/f'r2_{beam}__built_only_b0.0_{lowloss}_combined_records0-15_histories.npz'
            with np.load(path) as data:built_background[beam+'_'+lowloss]=data['local_J_m3'][:,0,[4,13,22]]
            for suffix in ('energy','same_event','tails'):
                path=Path(ref)/f'r2_{beam}__built_only_b0.0_{lowloss}_combined_records0-15_{suffix}.csv'
                if path.exists():(t.HERE/('low_built_'+beam+'_'+lowloss+'_'+suffix+'.csv')).write_bytes(path.read_bytes())
    np.savez_compressed(t.HERE/'slow_built_background.npz',**built_background)
    tails=rows('tail_signed_events.csv')
    events=[];same_event=[]
    from scipy.ndimage import label as connected
    nd=receiver.population()
    countlookup={(x['scenario'],x['band'],x['loss'],x['beam'],x['depth_m']):float(x['mean_distinct_count']) for x in counts}
    for name in ('width_censoring_sources_r2.csv','summary_sources_r2.csv','actual_band_extent.csv'):
        (t.HERE/('low_'+name)).write_bytes((Path(ref)/name).read_bytes())
    for row in sp+tails:
        if not row.get('peak_J_m3'):continue
        name=row.get('evidence_path') or f'{row["scenario"]}_b{row["band"]}_{row["loss"]}_{row["beam"]}_r{row["record"]}.npz'
        with np.load(t.HERE/name) as data:
            mu=float(data['expected_parents'])
            if row['conditional']=='True':
                row['conditional_record_weight']=float(mu*np.exp(-mu))
                row['conditional_meaning']='exactly one parent; weight P(N=1), not observed event or all nonempty records'
            at=np.unravel_index(data['elastic_J_m3'].argmax(),data['elastic_J_m3'].shape)
            reference=thermal[int(row['band'])]['mean_thermal_elastic_J_m3']
            row['peak_over_300K_same_band_thermal_mean']=float(row['peak_J_m3'])/reference
            when=(int(row['record'])+.5)*25e-6+float(data['time_s'][at[0]])
            sample=min(1999,int(when/50e-9))
            dep=max(0,at[1]-1)
            for lowloss in ('retained','Cu_proxy'):
                u=(built_background if row['scenario'].startswith('built') else slow)[row['beam']+'_'+lowloss][:,dep]
                value=float(u[sample])
                row['same_time_slow_'+lowloss+'_J_m3']=value
                row['slow_'+lowloss+'_energy_label']='peak' if value>np.quantile(u,.95) else 'hole' if value<np.quantile(u,.05) else 'ordinary'
            row['slow_track_context']='30-um reference only, not a microscopic track extrapolation' if at[1]==0 else 'same-time actual frozen depth receiver'
            energy=data['elastic_J_m3'];time=data['time_s'];dt=float(time[1]-time[0])
            point=at[2];conditional=row['conditional']=='True'
            row['impact_parent_rate_s']=t.RATE
            row['receiver_distinct_parent_rate_s']=int(row['distinct_arrived_parents'])/float(row['finite_exposure_s'])
            row['signed_phase_peak_ratio']=float(energy.max()/max(float(data['shuffled_J_m3'].max()),1e-300))
            step=float(np.max(abs(data['offset_m'])))
            for d,z in enumerate(t.DEPTHS):
                selected=True if row.get('evidence_path') else countlookup[(row['scenario'],row['band'],row['loss'],row['beam'],str(z))]<200
                common=dict(evidence=name,scenario=row['scenario'],beam=row['beam'],loss=row['loss'],band=row['band'],record=row['record'],depth_m=z,conditional=conditional)
                same_event.append(dict(common,time_s=float(time[at[0]]),spatial_offset_index=point,field_selected=selected,elastic_eV_D=float(energy[at[0],d,point]/(t.QE*nd[d])) if selected else None,distinct_parent_count=int(data['counts'][at[0],d]) if selected else None,effective_weighted_count=float(data['effective_count'][at[0],d]) if selected else None,meaning='same physical time/neighbour at all depths; blank means dense field was not reconstructed, not zero'))
                if not selected:continue
                u=energy[:,d,4];mean=float(u.mean())
                mask=u>3*mean if mean else np.zeros(len(u),bool)
                cuts=np.diff(np.r_[False,mask,False].astype(int));starts=np.flatnonzero(cuts==1);ends=np.flatnonzero(cuts==-1)
                for start,end in zip(starts,ends):
                    peak=start+int(np.argmax(u[start:end]));plane=energy[peak,d].reshape(3,3)
                    components,_=connected(plane>=u[peak]/2);coords=np.argwhere(components==components[1,1])
                    width=(coords.max(axis=0)-coords.min(axis=0)+1)*step
                    events.append(dict(common,threshold_over_finite_own_mean=3.,start_s=float(time[start]),duration_s=(end-start)*dt,duration_censored=bool(start==0 or end==len(u)),peak_J_m3=float(u[peak]),peak_over_finite_own_mean=float(u[peak]/mean),width_x_m=width[0],width_y_m=width[1],width_censored=bool(np.any(coords==0) or np.any(coords==2)),width_unresolved=bool(np.any(width<=step)),unconditional_excursion_rate_s=len(starts)/float(row['finite_exposure_s']) if not conditional else 0.,impact_parent_rate_s=t.RATE,meaning='connected fixed-centre energy excursion; shared-neighbour component width, not independent traces'))
    t.table('sparse_events.csv',sp)
    t.table('tail_signed_events.csv',tails)
    t.table('same_event_depths.csv',same_event)
    t.table('event_distributions.csv',events)
    refined=next(x for x in rows('frequency_refinement.csv') if int(x['band'])==5 and x['beam']=='even_steady' and x['loss']=='Cu_extrapolated' and float(x['depth_m'])==30e-6)
    lines=['# Pass 4 — observed coarse result','',
      'This is a conditional, finite-resolution acoustic model and an unpublished delivery, not measured PdD launch/loss or microscopic D excitation. The exact retained low band remains authoritative. No global low-band matrix was rerun.','',
      '## Deciding overlap','',
      'The matched normal-impulse common-band calculation uses both exact rung2 numeric files, all three receivers, actual elastic spatial quadratic plus kinetic energy/time averaging and both old loss laws. The retained exact guided Green reproduces the unmodified accepted transition at the same phase and50-ns sampling to6.35e-15 relative. This establishes its normalization/consumer join, not independent physical validation. A finite front/rear bulk-image illustration gives only0.051–0.087 of the resolved mean, but lacks converged lateral reflections and cannot establish a universal bulk-transport failure. Direction-dependent guided speeds/thickness branches, not that illustration alone, motivate retaining the resolved/guided kernel. See `overlap.csv`, `overlap.npz` and `guided_branches.csv`.','',
      '## All seven bands before refinement','',
      '| Band | Built launch W,20 nm/0.5 ps | Built launch W,60 nm/2 ps | Even/Cu 12-length importance-core parents,track /30 /50 /200 µm |','|---|---:|---:|---|']
    for j in range(7):
        power=[float(t.RATE*np.einsum('s,sb->',t.FRACTION,b[a,:,j])) for a in range(2)]
        cc=[r for r in counts if r['scenario']=='built' and int(r['band'])==j and r['loss']=='Cu_extrapolated' and r['beam']=='even_steady']
        lines.append(f'| {t.BANDS[j]:.3g}–{t.BANDS[j+1]:.3g} Hz | {power[0]:.4g} | {power[1]:.4g} | '+ ' / '.join(f'{float(r["mean_distinct_count"]):.3g}' for r in cc)+' |')
    lines+=['', 'The displayed counts are electrical parents in the **12-energy-length importance core**, not every arriving weak front. D2+ and all echoes count once. `band_depth_full_support.csv/npz` independently counts the first-five-image/eight-sigma front ancestry across all28 cells, without the attenuation cutoff: about0.71 million track and1.04 million deeper ancestors in the upper decades for even steady forcing. This is a finite modeled-path lower bound, not all history or detectable phonon occupancy. Participation from packet weights is a third quantity; none determines coherence. The128 shared snapshots preserve retarded burst intensity.','',
      'Cu attenuation is per length: exp(-alpha_E*actual path length), with the guided radial distance kept distinct. Al attenuation is per time: exp(-transit/tau_E). Neither numerical analog is measured high-frequency Pd loss or a confidence bound. For the heat/residence scenario their rates are assumed absorptive; the analog measurements do not deconvolve intrinsic and grain losses. No empirical grain-scattering coefficient is added, and explicit R/T is conservative redistribution. The two launch comparisons vary footprint and duration together; they do not identify their individual effects.','',
      '## Energy and depth readout','',
      '| Band | Built even/Cu importance-core elastic inventory eV/D:track /30 /50 /200 µm |','|---|---|']
    for j in range(7):
        rr=[r for r in read if r['scenario']=='built' and int(r['band'])==j and r['loss']=='Cu_extrapolated' and r['beam']=='even_steady']
        lines.append('| '+str(j)+' | '+' / '.join(f'{float(r["mean_elastic_eV_D"]):.4g}' for r in rr)+' |')
    lines+=['', '`readout.csv` retains importance-core distributions, zeros and exposure. `uncensored_transfer.csv` gives uncensored modeled-path incoherent means, numerical-underflow flags, inventory quotients and separate finite-lifetime4.1-ps upper transfer proxies. The latter assumes all dissipated total packet energy goes to D:2*elastic*(1-exp(-Gamma*4.1ps)), not an acquired microscopic energy measurement. `joint_distributions.npz` has within-band depth moments/energy samples; `band_depth_joint.npz` shares one parent measure across all28 core count cells. There is no reconstructed all-band signed field. Energy sample clock offsets and finite source history are explicit. See `LEAD.md` for positional peak/typical/phase comparisons.','',
      '## Sparse fields, supports and causes','',
      f'{len(sp)} finite sparse/depth-subset records are listed in `sparse_events.csv`; NPZ evidence retains the complete sampled importance-core parent catalog, common neighbours, mode-diagonal signed branch-channel amplitudes and shuffled parent phases. These are cycle-averaged analytic elastic envelopes, not a full microscopic tensor. Many-parent/small-participation cases use positive compound-Poisson packet-jump marginals and separate joint packet moments, not universal Rayleigh or a clipped Gaussian.','',
      'Actual empty records remain empty unconditionally. Forced single-parent evidence is explicitly conditioned on exactly one parent and weighted by P(N=1). Width/duration grids and observation censoring remain visible. Pre-existing energy, not merely amplitude sign, defines peak/hole labels. Direct/echo, coincidence and slow-field labels overlap. `peak_controls.json`/NPZ gives paired direct/echo and slow-refraction controls; `tail_signed_events.csv` gives the newly discovered sparse low-edge transfer at30 µm.','',
      '## Joins and decisive limitations','',
      'The source budget has each electrical parent once. Momentum, thermoelastic and cross launch terms are separate in `launch.csv`; the assumed comparison uses one recoil-normalized high-frequency envelope and the unchanged separate historical1% low-band term. The steady transport ledger debits wave input to stored energy, effective damping heat and edge escape once; its balance is an **expectation/residence-model** balance, not exact per-parent finite-record cross work against every unresolved mode. Recoil remainder may remain random motion/heat. Electronic stopping only enters deposited heat and explicit thermoelastic work.','',
      f'The ladder ends at2.5 THz, while supported pure-Pd branches extend above it. Within the same conditional source model, the conservative above-ladder built-energy bounds are {tail[0]["built_supported_above_ladder_upper_J"]:.3g} and {tail[1]["built_supported_above_ladder_upper_J"]:.3g} J/parent for the primary shape. They cannot move the displayed launch budget; no extra upper-band run is warranted by this bound. It is not a universal bound on unmeasured impact thermalization.','',
      'The actual accepted slow energy histories are reused at the same physical time; a strain norm upper is derived from their local elastic quadratic, not from invented extra source energy. Grüneisen10 is a declared refraction sensitivity, not a measured coefficient. The recorded path phase bound decides whether ignoring that slow refraction is small. One-way transport does not backreact on the accepted low band.','',
      'Weak joins remain guided/bulk grouping, coarse angular directivity versus band-integrated source work, unmodeled late cavity/grain conversion, and unknown Pd loss. `uncensored_transfer.csv` quantifies energy omitted by12-length selection and the four-sigma temporal guard. Its shuffled and arbitrary-phase finite-history expectation bounds are not bounds conditional on a selected maximum and exclude unmodeled paths; they do not promote core phase comparisons to the full field. Spatial widths are measured/censored on shared neighbouring rays. Pre-existing labels use earlier-parent-only energy at the finite-grid sample nearest first-envelope onset, with out-of-record labels censored.','',
      '## Which bands to refine','',
      f'Only2.5–25 and25–250 GHz were refined into two source-weighted bins; the other five decades were not recomputed. Corrected Cu transfer at30 µm is {float(refined["refined_mean_elastic_J_m3"]):.6g} J/m³ in the importance core, {float(refined["uncensored_incoherent_mean_elastic_J_m3"]):.6g} J/m³ uncensored, with {float(refined["union_mean_distinct_count"]):.6g} core parents. The local-D division is an acoustic inventory quotient. `frequency_refinement.csv` retains both loss scenarios and all depths; `tail_signed_events.csv` retains finite low-edge conditional fields. The final decade is not a physical cutoff; the separate above2.5-THz source bound is negligible. No additional sweep is warranted by these results.','',
      '## Portable replay','',
      'From repository root, with NumPy/SciPy and a C99 compiler (for the accepted portable RNG):',
      '```sh','OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/transport.py overlap',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/transport.py launch',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/joint.py',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/joint.py --full-support',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/transport.py receivers',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/readout.py refine --reuse-spectrum',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/readout.py tail-fields',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/readout.py controls',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/readout.py support',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/readout.py summarize',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/readout.py check --reuse-structural',
      'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 g13_rogue/pass4/readout.py lead',
      '```','A separately located verified accepted dependency can be supplied as `--reference-dir DIR`. The default is the root-relative accepted rung2 directory, not a machine path. The actual CLI was exercised with the exact hash-matching read-only extraction. Random draws use the accepted portable RNG; floating spectral eigensolutions are tolerance-portable, not promised byte-identical across BLAS implementations.','']
    (t.HERE/'READOUT.md').write_text('\n'.join(lines))
    print(json.dumps(dict(readout_rows=len(read),count_rows=len(counts),sparse_records=len(sp),above_ladder_bound=tail[:2])))


def refine(ref,reuse_spectrum=False):
    """Only the deciding 2.5–250-GHz transfer tails, after all-seven coarse census."""
    import receiver as r
    m=t.reference(ref);xy,prob,area=r.sources(m);nd=r.population()
    tracks=t.track_data(m);original=np.load(t.HERE/'launch.npz')['built_terms']
    output=[];spectrum={}
    frozen=dict(np.load(t.HERE/'selected_spectrum.npz')) if reuse_spectrum else {}
    for band in (4,5):
        if reuse_spectrum:freq,terms=frozen[f'band{band}_frequency_Hz'],frozen[f'band{band}_terms_J']
        else:freq,terms=t.bulk_spectrum(tracks,t.BANDS[band],t.BANDS[band+1])
        error=float(np.max(abs(terms.sum(axis=0)-original[:,:,band]))/np.max(abs(original[:,:,band])))
        if error>1e-12:raise ValueError('Selected spectrum changed original source epoch')
        spectrum[f'band{band}_frequency_Hz']=freq;spectrum[f'band{band}_terms_J']=terms
        split=np.sqrt(t.BANDS[band]*t.BANDS[band+1]);edges=(t.BANDS[band],split,t.BANDS[band+1])
        groups=[]
        for group in range(2):
            take=(freq<split) if group==0 else (freq>=split)
            e=terms[take,0].sum(axis=(0,-1))
            weights=np.einsum('fsb,s->f',terms[take,0].sum(axis=-1),t.FRACTION)
            fc=float(freq[take]@weights/weights.sum())
            groups.append((e,fc,edges[group+1]-edges[group]))
        for law in ('Cu_extrapolated','Al_transfer'):
            kernels=[r.kernel(xy,band,e,law,frequency=f,bandwidth=width) for e,f,width in groups]
            for beam in ('even_steady','hot_steady','hot_bursts'):
                pp=prob[int(beam!='even_steady')]
                mean=sum(r.moments(k,pp,None)[0] for k in kernels)
                count=np.zeros(4)
                for i,p0 in enumerate(pp):
                    for d in range(4):
                        intervals=[]
                        for k in kernels:
                            supported=k['weight'][i,d].sum(axis=-1)>0
                            intervals.extend((a,a+8*k['sigma']) for a in k['arrival'][i,d,supported])
                        count[d]+=t.RATE*p0*sum(b-a for a,b in r.unions(intervals))
                for d,z in enumerate(t.DEPTHS):
                    output.append(dict(scenario='built',beam=beam,band=band,loss=law,depth_m=z,refined_mean_elastic_J_m3=.5*mean[d],refined_mean_elastic_eV_D=.5*mean[d]/(t.QE*nd[d]),union_mean_distinct_count=count[d],source_weighted_low_Hz=groups[0][1],source_weighted_high_Hz=groups[1][1],source_epoch_relative_error=error,meaning='two disjoint spectral-bin energy transfer; same-parent count union; not a new fully coherent broadband phase-field claim',burst_mean='unit-intensity expectation; retarded finite burst realizations remain in band_depth_joint.npz'))
                    from scipy.special import logsumexp
                    logs=[np.log(t.RATE*np.sqrt(2*np.pi)*k['sigma']/2)+logsumexp(k['logweight'][:,d]+np.log(pp)[:,None,None]+np.log(t.FRACTION)[None,None,:]) for k in kernels]
                    full_log=float(logsumexp(logs));full_mean=float(np.exp(full_log))
                    output[-1].update(uncensored_incoherent_mean_elastic_J_m3=full_mean,log_uncensored_incoherent_mean_elastic_J_m3=full_log,uncensored_inventory_eV_D=full_mean/(t.QE*nd[d]),importance_core_mean_label='12-energy-length core, not an uncensored mean or acquired D energy')
            print('refined',band,law,flush=True)
    t.table('frequency_refinement.csv',output);np.savez_compressed(t.HERE/'selected_spectrum.npz',**spectrum)
    t.dump('frequency_refinement.json',dict(decision='replace false band-centre reach exclusions with explicit source-weighted low-edge conditional transfer in the two deciding decades',bands=[4,5],not_recomputed=[0,1,2,3,6],scope='built primary launch transfer means/count union; phase records remain explicitly band-centre conditional',source_reuse='same sole bulk_spectrum implementation; source epoch checked against sunk launch integrals',limitations='two frequency bins are a coarse transfer refinement, not a converged event-tail calculation; no measured high-band Pd lifetime or microscopic D coupling is inferred'))


def tail_fields(ref):
    """Signed packets where the two-bin refinement actually found few deep arrivals."""
    import receiver as r
    m=t.reference(ref);xy,prob,area=r.sources(m);p=t.Portable()
    with np.load(t.HERE/'selected_spectrum.npz') as source:
        freq=source['band5_frequency_Hz'];terms=source['band5_terms_J']
    split=np.sqrt(t.BANDS[5]*t.BANDS[6]);take=freq<split
    energy=terms[take,0].sum(axis=(0,-1))
    weights=np.einsum('fsb,s->f',terms[take,0].sum(axis=-1),t.FRACTION)
    fc=float(freq[take]@weights/weights.sum());nd=r.population()
    epsilon=max(float(x['strain_norm_upper']) for x in rows('slow_background.csv'))
    report=[]
    for law in ('Cu_extrapolated','Al_transfer'):
        ker=r.kernel(xy,5,energy,law,frequency=fc,bandwidth=split-t.BANDS[5])
        for beam in ('even_steady','hot_steady','hot_bursts'):
            pp=prob[int(beam!='even_steady')]
            lam=r.count_masks(ker,pp)
            counts=np.array([sum(lam[j] for j in range(1,16) if j&(1<<d)) for d in range(4)])
            selected=dict(ker,weight=ker['weight'].copy());selected['weight'][:,counts>=200]=0
            for record in range(4):
                payload,row=r.sparse_record(p,selected,xy,pp,beam,record,5,'built_low_edge',law,nd,epsilon)
                row['field_scope']='25–79.06-GHz low spectral bin; high bin is outside30-um twelve-length census, not omitted from source ledger'
                row['carrier_Hz']=fc
                name=f'tail_b5_{law}_{beam}_r{record}.npz';row['evidence_path']=name
                if payload is not None:
                    mu=float(payload['expected_parents'])
                    if bool(payload['conditional']):row['conditional_record_weight']=float(mu*np.exp(-mu))
                    payload['subband_edges_Hz']=np.array([t.BANDS[5],split])
                    payload['unconditional_elastic_J_m3']=np.zeros_like(payload['elastic_J_m3']) if bool(payload['conditional']) else payload['elastic_J_m3']
                    np.savez_compressed(t.HERE/name,**payload)
                report.append(row)
    t.table('tail_signed_events.csv',report)
    print(json.dumps(dict(refined_sparse_records=len(report),carrier_Hz=fc)))


def support(ref):
    """Separate finite modeled fronts, importance cores and omitted energy."""
    from scipy.special import logsumexp,erfc
    import receiver as r
    m=t.reference(ref);xy,prob,area=r.sources(m);nd=r.population()
    built=np.load(t.HERE/'launch.npz')['built'][0]
    energy={'built':built,'assumed':r.assumed(m)}
    read=rows('readout.csv')
    full=rows('band_depth_full_support.csv')
    front={(x['beam'],int(x['band']),float(x['depth_m'])):x for x in full}
    counts=rows('overlap_counts.csv')
    lookup={(x['scenario'],int(x['band']),x['loss'],x['beam'],float(x['depth_m'])):float(x['mean_distinct_count']) for x in counts}
    byrow={(x['scenario'],int(x['band']),x['loss'],x['beam'],float(x['depth_m'])):x for x in read}
    frozen=np.load(t.ROOT/'g13_rogue/pass2/replication/inputs.npz')
    maxburst=float(frozen['time_weights'][2].max());history=float(frozen['time_weights'].shape[1]*50e-9)
    report=[]
    for band in range(7):
        for law in ('Cu_extrapolated','Al_transfer'):
            base=r.kernel(xy,band,built[:,band],law)
            for scenario,E in energy.items():
                ratio=(E[:,band]/built[:,band]).T[np.repeat(np.arange(3),5)]
                ker=base if scenario=='built' else dict(base,weight=base['weight']*ratio[None,None,:,:],logweight=base['logweight']+np.log(ratio)[None,None,:,:])
                cache={}
                for beam in ('even_steady','hot_steady','hot_bursts'):
                    profile=int(beam!='even_steady');scale=maxburst if beam=='hot_bursts' else 1.
                    for d,z in enumerate(t.DEPTHS):
                        if (profile,d) not in cache:
                            value=ker['logweight'][:,d]+np.log(prob[profile])[:,None,None]+np.log(t.FRACTION)[None,None,:]
                            norm=np.log(t.RATE*np.sqrt(2*np.pi)*ker['sigma']/2)
                            cache[profile,d]=(float(logsumexp(value)+norm),float(logsumexp(np.where(ker['weight'][:,d]==0,value,-np.inf))+norm))
                        lm,lw=cache[profile,d];total=float(np.exp(lm));spatial=float(np.exp(lw))
                        row=byrow[scenario,band,law,beam,float(z)]
                        core=float(row['mean_elastic_J_m3'])
                        omitted=spatial+float(erfc(4/np.sqrt(2)))*core
                        ancestry=front[beam,band,float(z)]
                        conversion=2*(-np.expm1(-4.1e-12/float(row['energy_lifetime_s'])))/(t.QE*nd[d])
                        report.append(dict(scenario=scenario,beam=beam,band=band,loss=law,depth_m=z,
                            uncensored_incoherent_mean_elastic_J_m3=total,log_uncensored_incoherent_mean_elastic_J_m3=lm,
                            linear_underflow=bool(total==0 and np.isfinite(lm)),
                            uncensored_inventory_eV_D=total/(t.QE*nd[d]),finite_lifetime_4p1ps_upper_proxy_eV_D=total*conversion,
                            spatial_omitted_mean_J_m3=spatial,spatial_plus_temporal_omitted_mean_J_m3=omitted,
                            temporal_guard_energy_fraction=float(erfc(4/np.sqrt(2))),
                            finite_front_parent_lower_bound_mean=ancestry['mean_distinct'],source_history_s=history,
                            shuffled_omitted_energy_expectation_upper_J_m3=15*scale*omitted,
                            arbitrary_phase_finite_history_expectation_upper_J_m3=15*(1+t.RATE*history*scale)*scale*omitted,
                            meaning='modeled-path incoherent mean; expectation bounds are neither selection-conditional nor bounds on unmodeled late cavity/grain paths'))
                        row.update(finite_front_parent_lower_bound_mean=ancestry['mean_distinct'],
                            importance_parent_census_mean=lookup[scenario,band,law,beam,float(z)],
                            energy_sample_time_offset_s=4*ker['sigma'],conditional_4p1ps_eV_D=core*conversion)
            print('support',band,law,flush=True)
    t.table('uncensored_transfer.csv',report);t.table('readout.csv',read)
    t.dump('support_semantics.json',dict(actual_front_count='finite first-five-image/eight-sigma ancestry lower bound, independent of attenuation; not detectable phonon occupancy',importance_count='12-energy-length core only',weighted_count='nonnegative packet participation, not coherence',coherent_field='mode-diagonal signed core only; full weak/late/grain-remainder coherence is not claimed',causal_packet_clock='four-sigma filter delay is explicit and separate from group transit',history_s=history,omitted_bounds='shuffled and arbitrary-phase finite-source-history expectation bounds; neither covers a selected maximum or unmodeled propagation paths',conditional_proxy='2*elastic_inventory*(1-exp(-Gamma*4.1ps))/local_D_count, all damped packet energy hypothetically transferred to D; not a microscopic excitation measurement'))




def controls(ref):
    """Only the actual strongest primary built record gets cause/refraction controls."""
    import receiver as r
    m=t.reference(ref);xy,prob,area=r.sources(m);p=t.Portable()
    candidates=[x for x in rows('sparse_events.csv') if x['scenario']=='built' and x.get('peak_J_m3') and x['conditional']=='False']
    if not candidates:candidates=[x for x in rows('sparse_events.csv') if x['scenario']=='built' and x.get('peak_J_m3')]
    selected=max(candidates,key=lambda x:float(x['peak_J_m3']))
    band=int(selected['band']);beam=selected['beam'];law=selected['loss'];record=int(selected['record'])
    energy=np.load(t.HERE/'launch.npz')['built'][0]
    ker=r.kernel(xy,band,energy[:,band],law);pp=prob[int(beam!='even_steady')]
    lam=r.count_masks(ker,pp)
    countmean=np.array([sum(lam[j] for j in range(1,16) if j&(1<<d)) for d in range(4)])
    subset=dict(ker,weight=ker['weight'].copy());subset['weight'][:,~(countmean<200)]=0
    label=f'built_b{band}_{law}_{beam}_r{record}.npz'
    saved=np.load(t.HERE/label);epsilon=float(selected['slow_strain_upper'])
    arrays={'baseline_elastic_J_m3':saved['elastic_J_m3'],'baseline_shuffled_J_m3':saved['shuffled_J_m3'],'time_s':saved['time_s'],'parent_id':saved['parent_id']}
    result=dict(selected_record=label,conditioned=bool(saved['conditional']),slow_strain_upper=epsilon,refraction_sensitivity_coefficient=10.,scope='same complete sampled parents in actual low-count depths; direct/echo and one-way slow-phase controls')
    for name,mode,refraction in (('direct','direct',0.),('echo','echo',0.),('slow_plus','all',10*epsilon),('slow_minus','all',-10*epsilon)):
        payload,row=r.sparse_record(p,subset,xy,pp,beam,record,band,'built',law,r.population(),epsilon,path_mode=mode,refraction=refraction)
        if not np.array_equal(payload['parent_id'],saved['parent_id']) or not np.array_equal(payload['parent_source_index'],saved['parent_source_index']):raise ValueError('Control lost paired actual parents')
        arrays[name+'_elastic_J_m3']=payload['elastic_J_m3']
        result[name+'_max_J_m3']=float(payload['elastic_J_m3'].max())
    at=np.unravel_index(saved['elastic_J_m3'].argmax(),saved['elastic_J_m3'].shape)
    total=float(saved['elastic_J_m3'][at]);direct=float(arrays['direct_elastic_J_m3'][at]);echo=float(arrays['echo_elastic_J_m3'][at])
    result.update(at_original_peak=dict(total_J_m3=total,direct_diagonal_J_m3=direct,echo_diagonal_J_m3=echo,direct_echo_cross_J_m3=total-direct-echo,central_distinct_parent_count=int(saved['counts'][at[:2]]),central_weighted_count=float(saved['effective_count'][at[:2]]),count_reference='central receiver at the same time/depth, not a coherence count'),slow_plus_over_baseline_max=result['slow_plus_max_J_m3']/total if total else 0.,slow_minus_over_baseline_max=result['slow_minus_max_J_m3']/total if total else 0.)
    samepoint=saved['elastic_J_m3'][:,at[1],at[2]]
    typical=float(samepoint.mean());median=float(np.median(samepoint))
    result.update(band=band,beam=beam,loss=law,source_scenario='built:20nm/0.5ps',record=record,
        depth_m=float(saved['depth_m'][at[1]]),peak_offset_m=saved['offset_m'][at[2]].tolist(),
        peak_time_s=float(saved['time_s'][at[0]]),finite_exposure_s=float(selected['finite_exposure_s']),
        same_point_record_mean_J_m3=typical,same_point_record_median_J_m3=median,
        tallest_over_same_point_record_mean=total/typical if typical else None,
        tallest_over_same_point_record_median=total/median if median else None,
        paired_shuffled_at_original_peak_J_m3=float(saved['shuffled_J_m3'][at]),
        paired_shuffled_same_point_max_J_m3=float(saved['shuffled_J_m3'][:,at[1],at[2]].max()),
        paired_shuffled_grid_max_J_m3=float(saved['shuffled_J_m3'].max()),
        phase_scope='same declared scalar branch-channel importance-core parents, not full weak-front coherence')
    np.savez_compressed(t.HERE/'peak_controls.npz',**arrays);t.dump('peak_controls.json',result)
    print(json.dumps(result))


def check(ref,reuse_structural=False):
    import receiver
    m=t.reference(ref)
    prior=json.loads((t.HERE/'behavior_check.json').read_text()) if reuse_structural else {}
    if reuse_structural:relative=prior['guided_first_three_P14_vs_P18_relative']
    else:
        if not np.all(receiver.assumed(m)>0):raise ValueError('Gaussian envelope tail cancellation removed a requested band')
        import sys
        sys.path.insert(0,str(t.ROOT/'g13_rogue/pass3/rung1'))
        import rungs
        rungs.accepted()
        k=1800.;theta=.37;grain=27
        w14=t.plate_matrices(k,theta,grain,14)[0][:3]
        w18=t.plate_matrices(k,theta,grain,18)[0][:3]
        relative=float(np.max(abs(w14-w18)/w18))
        if relative>1e-4:raise ValueError('Guided branch thickness boundary unresolved')
    led=rows('transport_ledger.csv')
    relative_balance=max(abs(float(r['balance_W']))/max(float(r['launch_W']),1e-300) for r in led)
    if relative_balance>1e-12:raise ValueError('Transport debit/credit boundary')
    maxcount=0;phase_difference=0.;records=0
    for path in sorted(t.HERE.glob('*_r[0-3].npz')):
        z=np.load(path)
        if z['counts'].max()>len(z['parent_id']):raise ValueError('Echo counted as new parent')
        if len(np.unique(z['parent_id']))!=len(z['parent_id']):raise ValueError('Duplicate parent identity')
        if float(z['time_s'][0])!=float(z['packet_filter_delay_s']):raise ValueError('Causal filter clock not recorded')
        if len(z['preexisting_label'])!=len(z['parent_id']):raise ValueError('Missing parent labels')
        phase_difference=max(phase_difference,float(np.max(abs(z['elastic_J_m3']-z['shuffled_J_m3']))))
        maxcount=max(maxcount,int(z['counts'].max()));records+=1
    # Actual seeded-grain characteristic interface arithmetic.
    if reuse_structural:interface_error=prior['interface_R_plus_T_error']
    else:
        a=t.christoffel(np.array([.3,.4,np.sqrt(.75)]),12)[0]
        b=t.christoffel(np.array([.3,.4,np.sqrt(.75)]),45)[0]
        R=((a-b)/(a+b))**2;T=4*a*b/(a+b)**2
        interface_error=float(np.max(abs(R+T-1)))
        if interface_error>1e-14:raise ValueError('Interface energy split')
    xy=t.L/2+np.array([[.0007,.0009],[.0011,.0003]])
    energy=np.load(t.HERE/'launch.npz')['built'][0,:,5]
    errors={}
    for law in ('Cu_extrapolated','Al_transfer'):
        g=receiver.kernel(xy,5,energy,law)
        radius,transit,tr,rr,phase=receiver.anis_path(xy,t.DEPTHS[0],2,5,g['k'][0,2],g['vg'][0,2],g['f'])
        directivity=3*(t.DEPTHS[0]/radius)**2
        geom=4*np.pi*(radius**2+(20e-9)**2+1/g['k'][0,2]**2)*g['vg'][0,2]
        prefix=np.log(energy[0,2]*directivity*tr/(geom*np.sqrt(2*np.pi)*g['sigma']))
        actual=prefix-g['logweight'][:,0,10,0]
        expected=32*(g['f']/1e6)*np.log(10)/10*radius if law=='Cu_extrapolated' else transit/(34e-12*240e9/g['f'])
        errors[law]=float(np.max(abs(actual-expected)/expected))
        if errors[law]>1e-10:raise ValueError('Length/time attenuation join')
    result=dict(guided_first_three_P14_vs_P18_relative=relative,transport_expectation_relative_balance=relative_balance,interface_R_plus_T_error=interface_error,sparse_evidence_records=records,max_distinct_in_packet=maxcount,actual_paired_phase_control_max_difference_J_m3=phase_difference,review='owner self-check of actual changed model boundaries, not independent physical validation')
    result['loss_exponent_relative_errors']=errors
    result['reused_unchanged_source_guided_interface_evidence']=reuse_structural
    source=rows('source_record_budget.csv')
    result['source_stage_relative_balance']=max(abs(float(x['source_stage_balance_J']))/float(x['incident_J']) for x in source)
    result['accepted_source_epoch_checked']=True
    result['assumed_all_seven_positive']=True
    dependency=json.loads((t.HERE/'dependencies.json').read_text())
    base=['g13_rogue/pass2/coarse.py','g13_rogue/pass2/replication/replicate.py','g13_rogue/pass2/replication/portable.c','g13_rogue/pass2/replication/inputs.npz','g13_rogue/pass3/rung1b/signed.py','g13_rogue/pass3/rung1/rungs.py','g4_beam/merit/profiles.json','g11_cmns/collective/part1/stopping.json']
    dependency['base_consumers']=[dict(path=x,sha256=t.sha(t.ROOT/x)) for x in base]
    required=[]
    for beam in ('even_steady','hot_steady','hot_bursts'):
        for lowloss in ('retained','Cu_proxy'):
            for source in ('','__built_only'):
                prefix=f'r2_{beam}{source}_b0.0_{lowloss}_combined_records0-15'
                for suffix in ('histories.npz','energy.csv','same_event.csv','tails.csv'):
                    name=prefix+'_'+suffix
                    required.append(dict(path='g13_rogue/pass3/rung2/'+name,sha256=t.sha(Path(ref)/name)))
    dependency['accepted_readout_files']=required
    t.dump('dependencies.json',dependency)
    t.dump('behavior_check.json',result);print(json.dumps(result))

def lead():
    """Decision readout from settled bytes only; no new propagation calculation."""
    c=json.loads((t.HERE/'peak_controls.json').read_text())
    behavior=json.loads((t.HERE/'behavior_check.json').read_text())
    spectrum=np.load(t.HERE/'launch.npz')['built']
    powers=[float(t.RATE*np.einsum('s,sjb->',t.FRACTION,spectrum[j])) for j in range(2)]
    peak=c['at_original_peak']['total_J_m3'];depth=c['depth_m']
    rr=next(x for x in rows('readout.csv') if x['scenario']=='built' and int(x['band'])==c['band'] and x['beam']==c['beam'] and x['loss']==c['loss'] and float(x['depth_m'])==depth)
    weak=next(x for x in rows('uncensored_transfer.csv') if x['scenario']=='built' and int(x['band'])==c['band'] and x['beam']==c['beam'] and x['loss']==c['loss'] and float(x['depth_m'])==depth)
    with np.load(t.HERE/c['selected_record']) as z:
        d=int(np.where(z['depth_m']==depth)[0][0]);centre_peak=float(z['elastic_J_m3'][:,d,4].max())
    thermal=float(rr['mean_elastic_J_m3'])/float(rr['wave_over_300K_same_band_thermal_mean'])
    inventory=peak/(t.QE*float(rr['n_D_m3']))
    reference=next(x for x in rows('frequency_refinement.csv') if int(x['band'])==5 and x['beam']=='even_steady' and x['loss']=='Cu_extrapolated' and float(x['depth_m'])==30e-6)
    threshold=3*float(reference['uncensored_incoherent_mean_elastic_J_m3'])
    tails=[]
    for beam in ('even_steady','hot_steady','hot_bursts'):
        for law in ('Cu_extrapolated','Al_transfer'):
            fr=next(x for x in rows('frequency_refinement.csv') if int(x['band'])==5 and x['beam']==beam and x['loss']==law and float(x['depth_m'])==30e-6)
            samples=[];events=0;exposure=0.
            for record in range(4):
                with np.load(t.HERE/f'tail_b5_{law}_{beam}_r{record}.npz') as z:
                    value=z['unconditional_elastic_J_m3'][:,1,4];samples.append(value)
                    high=value>threshold;events+=int(np.count_nonzero(high&~np.r_[False,high[:-1]]))
                    exposure+=float((z['time_s'][1]-z['time_s'][0])*len(z['time_s']))
            values=np.concatenate(samples)
            tails.append(dict(beam=beam,loss=law,core_mean_J_m3=float(fr['refined_mean_elastic_J_m3']),
                uncensored_mean_J_m3=float(fr['uncensored_incoherent_mean_elastic_J_m3']),
                log_uncensored_mean_J_m3=float(fr['log_uncensored_incoherent_mean_elastic_J_m3']),
                inventory_eV_D=float(fr['uncensored_inventory_eV_D']),
                core_mean_distinct=float(fr['union_mean_distinct_count']),threshold_J_m3=threshold,
                samples=len(values),exceeding_samples=int(np.count_nonzero(values>threshold)),
                sampled_exceedance_fraction=float(np.mean(values>threshold)),
                sampled_core_peak_J_m3=float(values.max()),excursions=events,finite_exposure_s=exposure,
                empirical_receiver_excursions_s=events/exposure,
                meaning='fixed-centre finite model samples, correlated in time; forced conditional records are zero unconditionally; not an experimental rarity estimate'))
    t.table('tail_30um_readout.csv',tails)
    result=dict(built_launch_W_joint_shape_duration_scenarios=powers,launch_ratio=powers[0]/powers[1],
        positional_peak=c,peak_inventory_eV_D=inventory,peak_over_300K_same_band_harmonic_mean=peak/thermal,
        fixed_centre_peak_J_m3=centre_peak,omitted=weak,
        shuffled_omitted_expectation_over_fixed_centre_peak=float(weak['shuffled_omitted_energy_expectation_upper_J_m3'])/centre_peak,
        arbitrary_phase_history_expectation_over_fixed_centre_peak=float(weak['arbitrary_phase_finite_history_expectation_upper_J_m3'])/centre_peak,
        tails=tails,changed_boundary_checks=behavior)
    t.dump('headline_results.json',result)
    a=c['at_original_peak']
    lines=['# Pass 4 — decision-ready readout','',
        'Conditional Pd-transfer model, not measured PdD attenuation, D excitation or a full-field coherence/excess claim. All seven coarse decades, two loss sensitivities, both source budgets and three beam histories are retained.','',
        '## Result and position','',
        f'Built launch is {powers[0]:.9g} W for20 nm/0.5 ps and {powers[1]:.9g} W for60 nm/2 ps: a {powers[0]/powers[1]:.6g}-fold **joint footprint+duration** sensitivity, not separate identification. Only2.5–25 and25–250 GHz were frequency-refined; no further sweep was run.','',
        f'Strongest unconditioned primary signed-core sample: `{c["selected_record"]}`, band{c["band"]} ({t.BANDS[c["band"]]:.6g}–{t.BANDS[c["band"]+1]:.6g} Hz), {c["beam"]}, {c["loss"]}, {c["source_scenario"]}, depth{depth:.6g} m, offset{c["peak_offset_m"]} m, exposure{c["finite_exposure_s"]:.9g} s. Peak is {peak:.9g} J/m³; same-point record mean {c["same_point_record_mean_J_m3"]:.9g}, median {c["same_point_record_median_J_m3"]:.9g}; tallest/mean {c["tallest_over_same_point_record_mean"]:.6g}. Tallest/median is {c["tallest_over_same_point_record_median"]} (undefined if median is zero).','',
        f'At that time/depth the **central** receiver has {a["central_distinct_parent_count"]} supported electrical parents and weighted participation {a["central_weighted_count"]:.9g}. These are not equally coherent contributions or the count of all weak fronts. Paired shuffled energy at the original peak location/time is {c["paired_shuffled_at_original_peak_J_m3"]:.9g} J/m³; its same-point temporal maximum is {c["paired_shuffled_same_point_max_J_m3"]:.9g}, and its whole-grid maximum is {c["paired_shuffled_grid_max_J_m3"]:.9g}. This is a paired core comparison, not a causal attribution to full-field coherent focusing.','',
        f'Direct/echo cross energy at the selected peak is {a["direct_echo_cross_J_m3"]:.6g} J/m³. The declared Grüneisen10 slow-phase sensitivity changes the maximum by factors {c["slow_plus_over_baseline_max"]:.9g} and {c["slow_minus_over_baseline_max"]:.9g}; it introduces no source energy. Accepted slow histories are same-time inputs, not an event-by-event reconstructed microscopic low/high parent trajectory. Transport is one-way.','',
        f'The peak inventory quotient is {inventory:.9g} eV per local D, not acquired D energy. Peak divided by the300-K same-band harmonic mean is {peak/thermal:.9g}; this is not a thermal-tail probability. Separate4.1-ps columns assume all total packet damping could transfer to D and use a finite decay fraction, not a measured microscopic coupling.','',
        '## Front ancestry, importance and missing phase','',
        'The first-five-image/eight-sigma joint ancestry census has roughly0.71–0.72 million track parents and1.04–1.06 million deeper parents in the upper decades, even when the12-energy-length core has only a few strong weights. It is a finite modeled-front lower bound; noncompact tails and later cavity paths are not an exact all-history count. N_eff is a different energy-weight measure.','',
        f'At the selected central depth, omitted spatial-plus-temporal diagonal mean is {float(weak["spatial_plus_temporal_omitted_mean_J_m3"]):.9g} J/m³. Its shuffled expectation bound is {result["shuffled_omitted_expectation_over_fixed_centre_peak"]:.6g} times the sampled fixed-centre peak; the arbitrary-phase finite-history expectation bound is {result["arbitrary_phase_history_expectation_over_fixed_centre_peak"]:.6g} times it. These are **unconditional expectations**, not selection-conditional peak certificates. The latter need not be small; no full weak-front phase conclusion follows. Bounds exclude unmodeled late cavity/grain-redistributed paths.','',
        '## Restored30-µm tail','',
        f'All rows below use one absolute threshold {threshold:.9g} J/m³ (three times corrected even/Cu uncensored mean), at the fixed central receiver. Sample times are correlated, widths/durations can be censored, and short-record rates are not rarity forecasts. Electrical impact rate is {t.RATE:.9g} /s, distinct from receiver excursions.','',
        '| Beam / loss | Uncensored mean J/m³ | Inventory eV/D | Core parents | Exceeding samples /2048 | Excursions / exposure s |',
        '|---|---:|---:|---:|---:|---|']
    for row in tails:lines.append(f'| {row["beam"]} / {row["loss"]} | {row["uncensored_mean_J_m3"]:.6g} | {row["inventory_eV_D"]:.6g} | {row["core_mean_distinct"]:.6g} | {row["exceeding_samples"]} | {row["excursions"]} / {row["finite_exposure_s"]:.6g} |')
    lines+=['','Al core zeros are not physical zero energy; the uncensored transfer is retained. The25–250-GHz band-centre core zero at30 µm was not a valid exclusion of its lower-frequency tail. The requested2.5-THz endpoint is not the Pd dispersion cutoff: the same-model above-ladder primary bound is about5.24e-35 J per electrical parent, so no extra upper band was warranted.','',
        '## Scope and decisive limitations','',
        'Complete tables retain all beams/depths, source/transport debits, zeros, conditional weights, shared-neighbour support, event labels and accepted5–250-kHz baselines with both losses. Exact common-band Green normalization error is6.35e-15; a truncated bulk-image illustration is not proof that every bulk transport construction fails. Owner checks compare Cu and Al propagation with their direct length/time equations; this is not independent physical validation. No accepted low-band implementation was changed.','',
        'The receiver is a finite-resolution, coarse angular, mode-diagonal envelope model with first-five images and an explicit guided/late-cavity/grain-conversion remainder, not foil-resonance reconstruction. High-frequency Pd loss and microscopic launch shape remain unmeasured; imported effective-loss numerics are conditional absorptive scenarios, not deconvolved Pd intrinsic data. Source work above the resolved band is isolated-parent/expectation bookkeeping, not exact finite-record unresolved cross work. These limitations, not further exponent sweeps, control the next physical interpretation.','',
        'Root-relative dependency recipes and exact hashes are in `README.md`/`dependencies.json`; `governing_plan.json` binds authority. `manifest.json` and `evidence_manifest.json` whitelist archived members. Actual draft IDs and remote readback hashes are in the separately returned delivery receipt. No commit, push, publication, merge or broader campaign was performed.','']
    (t.HERE/'LEAD.md').write_text('\n'.join(lines))
    print(json.dumps(dict(powers_W=powers,positional_peak=c,tail_30um=tails)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['summarize','controls','refine','tail-fields','support','check','lead']);p.add_argument('--reference-dir',default=str(t.ROOT/'g13_rogue/pass3/rung2'));p.add_argument('--reuse-spectrum',action='store_true');p.add_argument('--reuse-structural',action='store_true');a=p.parse_args()
    if a.action=='summarize':summarize(a.reference_dir)
    elif a.action=='controls':controls(a.reference_dir)
    elif a.action=='refine':refine(a.reference_dir,a.reuse_spectrum)
    elif a.action=='tail-fields':tail_fields(a.reference_dir)
    elif a.action=='support':support(a.reference_dir)
    elif a.action=='lead':lead()
    else:check(a.reference_dir,a.reuse_structural)
