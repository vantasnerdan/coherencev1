"""Pass 5: signed whole-foil conditional receiver -> resolved medium -> G12.

--case is the bounded joined consumer. --all --stage runs the same code over
all selected cases. --consume reads only selected saved arrays and inputs.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from common import HERE,BASE,MD,QE,KB,ENERGIES,inherited,calendar,emit,digest


def first_case(path4):
    from receiver import Green
    from live import Strain
    from sampling import marks,distribution,frequency_readout
    from medium import Medium
    from yield_consumer import difference
    started=time.monotonic();g=Green(HERE,path4,2.5);clock=calendar(path4)
    point=np.array([.0025,-.0025,1e-8]);mark=marks(g,point,'hot_bursts');pdf=distribution(mark,clock,'hot_bursts',samples=64)
    sampler=Strain(HERE,path4);xyz=np.tile(point,(32,1));xyz[:,0]+=np.linspace(-1e-8,1e-8,32);ts=np.linspace(.1,1.99,32)
    start=time.monotonic();field=sampler.bands(2.5,'hot_bursts',xyz,time_s=ts);sampler_time=time.monotonic()-start
    low_mean=sampler.low.mean(2.5,'hot_bursts',xyz,ts)
    strain=field.sum(axis=1).real+low_mean;medium=Medium(HERE);T,U,gradT,gradU=medium.sample(2.5,'hot_bursts',xyz,time_s=ts)
    omega,vg,kdot=medium.optical(2.5,'hot_bursts',xyz,np.array([1e8,0.,0.]),time_s=ts,strain=strain/3)
    quiet=medium.optical(2.5,'hot_bursts',xyz,np.array([1e8,0.,0.]),time_s=ts,strain=0.)[0]
    theta=[]
    for j,t in enumerate(T):
        coefficient=np.array([g.response.covariance(g.arm,ch['frequency'],float(t))['theta_per_linear_strain_eV'] for ch in g.channels])
        theta.append(float(np.real(coefficient@field[j])))
    theta=np.array(theta);zero=np.zeros((32,3));kernel=difference(np.full(32,5.),KB*T,zero,KB*T+theta,zero,'p',.04)
    proof=json.loads((HERE/'passive_controls.json').read_text())
    low_pdf=sampler.low.mean(2.5,'hot_bursts',point,pdf['time_s'])
    pdf_total=pdf['strain']+low_pdf
    gamma_fraction=frequency_readout(medium,medium.sample(2.5,'hot_bursts',np.broadcast_to(point,(64,3)),time_s=pdf['time_s']),pdf_total,'onsite')[0]
    import meetings
    from yield_consumer import build_acoustic_probe,consume
    encounter=meetings.build(HERE,path4,medium.sample,sampler.sample,probe_samples=16,event_samples=128,event_response_samples=16,
        case_selection={(2.5,'nominal_retained','onsite','hot_bursts')},prefix='first_')
    build_acoustic_probe(HERE,path4,prefix='first_')
    yield_rows=consume(HERE,path4,prefix='first_')
    identity={p.name:digest(p) for p in sorted(HERE.glob('*.py'))}
    result=dict(base=BASE,source_identity=identity,physical_seconds=2.,source_nodes=len(mark['h']),corrected_AC_filter=True,
        filter_L2=proof['isolated_Hann_L2'],filter_zero_DC_abs=proof['isolated_Hann_zero_DC_quadrature_abs'],
        isolated_port_real_flux_quadratures=proof['isolated_Hann_real_flux_quadratures'],isolated_port_pseudocovariance_abs=proof['isolated_Hann_pseudocovariance_abs'],
        force_to_strain_phase='bulk +i; outgoing guide +pi/4; source/receiver modal gauge checked, residual guide source phase is a scenario',
        source_terminal_range_m=[r-2e-9 for r in g.source_extent],
        signed_strain_min=float(pdf_total.min()),signed_strain_max=float(pdf_total.max()),
        centered_sample_rms=float(np.std(pdf_total-pdf['mean']-low_pdf)),Campbell_high_variance=float(mark['weight']@(mark['h']**2)),
        low_prompt_pressure_innovation_rms_bound=float(sampler.low.innovation_rms_bound(2.5,'hot_bursts',point)),
        low_prompt_pressure_Cox_mean_abs_max=float(abs(low_pdf).max()),
        actual_high_Cox_mean_abs_max=float(abs(pdf['mean']).max()),sampler_points=32,sampler_seconds=sampler_time,
        fractional_Gamma_frequency_min=float(gamma_fraction.min()),fractional_Gamma_frequency_max=float(gamma_fraction.max()),
        encounter_cases=len(encounter['cases']),actual_G12_case_rows=len(yield_rows),
        encounter_array_sha256=digest(HERE/'first_meetings.npz'),actual_G12_sha256=digest(HERE/'first_yield_outcomes.csv'),
        maximum_fractional_strain_Omega=float(abs(omega/quiet-1).max()),
        signed_velocity_covariance_delta_min_eV=float(theta.min()),signed_velocity_covariance_delta_max_eV=float(theta.max()),
        target_axis_velocity_rms_min_m_s=float(np.sqrt((KB*T+theta).min()*QE/MD)),target_axis_velocity_rms_max_m_s=float(np.sqrt((KB*T+theta).max()*QE/MD)),
        connected_G12_kernel_min_m2=float(kernel.min()),connected_G12_kernel_max_m2=float(kernel.max()),
        total_seconds=time.monotonic()-started,conditional_model=True,
        limitations='Conditional source-interior/Fresnel and within-band representation; low thermoelastic continuation, exterior Gaussian phase/return covariance, and microscopic contact remain unclosed. No physical tail absence or contact probability inferred.')
    medium.close();emit(HERE/'first_integrated.json',result);print(json.dumps(result),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare',metavar='PASS4')
    parser.add_argument('--case',action='store_true')
    parser.add_argument('--native-case',action='store_true')
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--stage',choices=('receivers','meetings','probe','native','native-observers','assemble','controls','maps','report'))
    parser.add_argument('--energy',type=float,choices=ENERGIES)
    parser.add_argument('--samples',type=int,default=256)
    parser.add_argument('--consume',action='store_true')
    parser.add_argument('--package',action='store_true');parser.add_argument('--cold',action='store_true');parser.add_argument('--upload',action='store_true')
    args=parser.parse_args()
    import delivery
    if args.prepare:delivery.prepare(args.prepare)
    path4=HERE/'inherited'
    if args.case:
        from receiver import build_spectra
        build_spectra(HERE,path4,energies=(2.5,),closures=('nominal_retained',))
        first_case(path4)
    if args.native_case:
        import native_gammax
        from yield_consumer import build_acoustic_probe,consume
        native_gammax.controls(HERE,path4);native_gammax.build(HERE,path4,2.5,narrow=True)
        build_acoustic_probe(HERE,path4,prefix='first_native_')
        rows=consume(HERE,path4,prefix='first_native_')
        print(json.dumps(dict(native_case=True,actual_G12_rows=len(rows),source='frozen onsite incoming launch, original pressure',maximum_abs_fraction=max(abs(r['fraction_of_G12_residual']) for r in rows))),flush=True)
    if args.all:
        if args.stage is None:parser.error('--all requires a bounded --stage; every physical job has a two-hour limit')
        if args.stage=='receivers':
            import sampling
            energies=(args.energy,) if args.energy is not None else ENERGIES
            prefix=f'E{args.energy:g}_' if args.energy is not None else ''
            sampling.build(HERE,path4,energies=energies,samples=args.samples,prefix=prefix)
        elif args.stage=='meetings':
            from medium import Medium
            from live import Strain
            import meetings
            from native_gammax import CachedStrain
            medium=Medium(HERE);sampler=Strain(HERE,path4)
            selected=None if args.energy is None else {(args.energy,c,a,b) for c in meetings.ARMS for a in meetings.ALLOCATIONS for b in meetings.BEAMS}
            prefix='' if args.energy is None else f'E{args.energy:g}_'
            cached=CachedStrain(sampler,HERE,prefix)
            result=meetings.build(HERE,path4,medium.sample,cached.sample,case_selection=selected,prefix=prefix)
            medium.close();print(json.dumps(dict(meetings_ready=True,cases=len(result['cases']))),flush=True)
        elif args.stage=='probe':
            from yield_consumer import build_acoustic_probe
            build_acoustic_probe(HERE,path4,prefix='' if args.energy is None else f'E{args.energy:g}_')
        elif args.stage=='native':
            if args.energy is None:parser.error('--stage native requires a bounded --energy')
            import native_gammax
            native_gammax.build(HERE,path4,args.energy)
        elif args.stage=='assemble':
            import meetings,analyze
            from yield_consumer import merge_acoustic_probe
            import native_gammax
            result=meetings.merge(HERE);merge_acoustic_probe(HERE);analyze.merge_receivers(HERE)
            native_gammax.assemble(HERE);native_gammax.receiver_projection(HERE,path4)
            print(json.dumps(dict(assembled_cases=len(result['cases']))),flush=True)
        elif args.stage=='controls':
            import analyze
            analyze.sensitivities(HERE,path4)
        elif args.stage=='maps':
            import analyze
            analyze.maps(HERE,path4)
        elif args.stage=='report':
            import analyze
            analyze.source_ledger(HERE,path4);analyze.limits(HERE,path4)
            print(json.dumps(analyze.report(HERE)),flush=True)
        elif args.stage=='native-observers':
            import analyze,native_gammax
            native_gammax.tagged_sensitivity(HERE,path4);analyze.native_landscape(HERE,path4)
    if args.consume:
        import native_gammax
        rows=native_gammax.consume_all(HERE,path4)
        selected=json.loads((HERE/'selected_inputs.json').read_text())
        for row in selected['files']:
            if digest(path4/row['path'])!=row['sha256']:raise ValueError('Selected inherited input changed')
        low_inputs=json.loads((HERE/'low_inputs.json').read_text())
        if digest(HERE/'low_modes.npz')!=low_inputs['cache_sha256']:raise ValueError('Selected low-mode projection changed')
        native_inputs=json.loads((HERE/'native_Gamma_X_inputs/provenance.json').read_text())
        for row in native_inputs['files']:
            if digest(HERE/row['local'])!=row['sha256']:raise ValueError('Selected native Gamma/X input changed')
        import scipy
        emit(HERE/'cold_consume.json',dict(success=True,base=BASE,rows=len(rows),selected_inputs_verified=len(selected['files']),native_input_commit=native_inputs['commit'],
            consumer_identity={p.name:digest(p) for p in sorted(HERE.glob('*.py'))},runtime=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
            meeting_sha256=digest(HERE/'meetings.npz'),native_meeting_sha256=digest(HERE/'nativeGX_meetings.npz'),medium_sha256=digest(HERE/'medium.npz'),acoustic_probe_sha256=digest(HERE/'acoustic_probe.npz')))
        print(json.dumps(dict(consume=True,rows=len(rows),selected_only=True)),flush=True)
    if args.package:print(json.dumps(delivery.package()),flush=True)
    if args.cold:print(json.dumps(delivery.cold(sys.executable)),flush=True)
    if args.upload:print(json.dumps(delivery.upload()),flush=True)


if __name__=='__main__':main()
