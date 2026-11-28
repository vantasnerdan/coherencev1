"""One joined entrypoint, with disjoint case writers and a single reducer."""
import concurrent.futures,json
import numpy as np
import sources as s
import fields as f
import native


def case_job(arguments):
    import run
    arms,beams,ec=arguments;calendar=f.Calendar();tracks=json.loads((s.HERE/'frozen/tracks.json').read_text())
    merged={key:[] for key in ('energy','bands','checks','travel','fine','diagnostics','optical','moments','shuffled')}
    for arm in arms:
        for beam in beams:
            controls=('temperature_off','strain_off','frozen','native_e_off') if arm['closure']=='nominal_retained' and ec==2.5 else ()
            result=run.run_case(arm,beam,ec,calendar,tracks,controls)
            for key in merged:merged[key].extend(result[key])
    return merged

def shared_band_rows(rows):
    """Lossless compaction of identical physical rows, retaining every case."""
    groups={}
    for row in rows:
        body={k:v for k,v in row.items() if k not in ('arm','allocation','beam','control')}
        key=json.dumps(body,sort_keys=True,separators=(',',':'))
        case=f'E{row["E_cm_keV"]:g}:{row["arm"]}:{row["beam"]}:{row["control"]}'
        if key not in groups:groups[key]=(body,[])
        if case not in groups[key][1]:groups[key][1].append(case)
    return [dict(body,applies_to_cases=';'.join(cases)) for body,cases in groups.values()]


def all_fields(model,workers=3,closure=None):
    arms=[a for a in s.ARMS if closure is None or a['closure']==closure]
    # Source preparation is an explicit single-writer stage (--native), never a
    # side effect of a field worker or a replay of already accepted packets.
    native.load()
    # One writer creates shared paths before disjoint workers read them.
    for arm in arms:model.carrier_kernel(arm)
    calendar=f.Calendar();clock=np.linspace(0,2,100001)
    np.savez_compressed(s.HERE/'calendar.npz',edges_s=calendar.edges,burst_on=calendar.state,phase=calendar.phase,normalization=calendar.norm,
        sample_time_s=clock,intensity_hot_bursts=calendar.value(clock,'hot_bursts'),seconds=2.,physical_parent_seed=5151)
    jobs=[([a for a in arms if a['closure']==name],beams,ec)
        for name in dict.fromkeys(a['closure'] for a in arms)
        for ec in s.ENERGIES for beams in (('even_steady',),('hot_steady','hot_bursts'))]
    merged={key:[] for key in ('energy','bands','checks','travel','fine','diagnostics','optical','moments','shuffled')}
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(case_job,jobs):
            for key in merged:merged[key].extend(result[key])
    if closure is not None:
        print(json.dumps(dict(completed_closure=closure,cases=len(arms)*len(s.BEAMS)*len(s.ENERGIES),maximum_relative_energy_defect=max(r['relative_balance'] for r in merged['checks']))),flush=True)
        return
    merged['bands']=shared_band_rows(merged['bands'])
    for name,key in (('energy_history.csv','energy'),('band_transport.csv','bands'),('transport_scales.csv','travel'),
                     ('phase_extrema.csv','fine'),('medium_controls.csv','diagnostics'),('conditioned_optical.csv','optical'),
                     ('parent_moments.csv','moments'),('phase_shuffle.csv','shuffled')):model.table(name,merged[key])
    model.emit('conservation.json',merged['checks']);mixture,consumer=model.consume();model.figures()
    result=dict(model='joined whole-foil conservative FV enthalpy, reciprocal 303-mode Hamiltonian, retarded directional carriers and seven-band ports, physical conditioned phase windows',
        seconds=2.,beams=list(s.BEAMS),arms=[a['name'] for a in s.ARMS],energies_keV=list(s.ENERGIES),
        conditioned_Poisson_phase_window_E_cm_keV=2.5,source_covariance_histories_E_cm_keV=list(s.ENERGIES),
        maximum_relative_macro_energy_defect=max(r['relative_balance'] for r in merged['checks']),
        minimum_signed_wave_strain=min(r['signed_strain_hole'] for r in merged['fine']),maximum_signed_wave_strain=max(r['signed_strain_peak'] for r in merged['fine']),
        maximum_primary40_fraction_of_residual=max(r['fraction_of_G12_residual'] for r in mixture if r['branch']=='p' and r['comparator']=='static40'),
        cold_consume='run.py --consume reads only shipped projections, source covariances and built whole-foil fields',
        primary_omissions=['native q calibration, onsite/neighbour partition and Gaussian/focus closure',
                          'measured loaded-PdD transport/dephasing and microscopic launch phase',
                          'native neighbour-shell to nuclear contact continuation'])
    model.emit('results.json',result);print(json.dumps(result),flush=True)
