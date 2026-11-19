"""Consume matched trajectories, native target quadrature and conditional pair law.

Run after trajectories.py. Outputs are deterministic finite-window estimates,
not error bars on a nuclear tail or first-principles material contact rates.
"""
from pathlib import Path
import csv
import importlib.util
import json
import sys
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from g4_beam.forward import model as M
from g4_beam.forward.consumer import forward_yields
from g4_beam.merit.api import reference
from response import cross_section_adapter, LAB, MD, QE


def save_csv(name, rows):
    with (HERE/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def folded_cross_section(lab_keV, velocity, u_keV, channel, reaction, direction=2):
    """Empirical <sigma(Erel)*vrel/vbeam>, mb; DD target isotope.

    Direction is a Cartesian axis; both pump-parallel x and transverse z are used.
    No starting-distance tunnelling multiplier accompanies a cross section
    already defined for an incoming scattering state.
    """
    speed = np.sqrt(2*np.asarray(lab_keV)*1000*QE/MD)
    v = np.asarray(velocity, float)
    # Keep joint v_parallel and v_perpendicular, rather than independent tails.
    parallel = v[:,direction]
    transverse2 = np.sum(v*v, axis=1)-parallel*parallel
    result = np.empty(len(speed))
    for first in range(0, len(speed), 24):
        speeds = speed[first:first+24,None]
        rel2 = (speeds-parallel[None,:])**2 + transverse2[None,:]
        relative = np.sqrt(rel2)
        ec = MD*rel2/(4*1000*QE)
        result[first:first+24] = np.mean(M.reaction_sigma(ec, u_keV, channel, reaction)*relative/speeds, axis=1)
    return result


def project_motion(velocity, channel, reaction, *, direction=2, kernel_grid=None):
    target = reference()
    grid = np.geomspace(.1, 16., 250) if kernel_grid is None else kernel_grid
    raw = M.reaction_sigma
    average = folded_cross_section(grid, velocity, .04, channel, reaction, direction)
    ratio = average/raw(grid/2, .04, channel, reaction)
    def transform(original, ec, u, ch, re):
        if ch != channel or re != reaction:
            raise ValueError('Kernel branch mismatch')
        return original(ec, u, ch, re)*np.interp(2*np.asarray(ec), grid, ratio)
    with cross_section_adapter(transform):
        return forward_yields(LAB/2, target.profiles, Ue_keV=.04,
                              measurement={'channel':channel},
                              stopping_fn=target.stopping_fn, reaction=reaction)


def pair_summary(item, data, p):
    """Resolve a connected radial width only above accumulator roundoff."""
    s = item['statistics']
    mean = float(np.mean(data['pair_mean_r_A']))
    second = float(np.mean(data['pair_mean_r2_A2']))
    variance = second-mean*mean
    resolution = 8*s['observation_samples']*np.finfo(float).eps*max(second,mean*mean)
    if variance < -resolution:
        raise ArithmeticError('Pair second moment is inconsistent beyond accumulator roundoff')
    resolved = variance > resolution
    support_min = p['a_A']/np.sqrt(2)-s['maximum_observed_site_displacement_A']
    turning_max = 14.399645478/50
    assert support_min > turning_max, 'Observed displacement reaches the finite-start barrier; unity fold no longer justified'
    return dict(case=item['summary']['label'],reference_nn_A=p['a_A']/np.sqrt(2),
        record_ps=s['record_ps'],sampled_joint_pairs=s['sampled_joint_pairs'],
        minimum_A=s['min_pair_A'],median_pair_min_A=s['pair_min_quantiles_A'][3],
        mean_pair_distance_A=mean,
        pooled_radial_std_A=float(np.sqrt(variance)) if resolved else None,
        radial_variance_accumulator_resolution_A2=resolution,
        radial_std_status='resolved' if resolved else 'below accumulated-second-moment resolution',
        radial_std_observed_support_upper_A=s['maximum_observed_pair_displacement_A'],
        conservative_beam_start_at_nn_support_min_A=support_min,
        largest_unscreened_turning_radius_A_at_consumer_cutoff=turning_max,
        starting_separation_barrier_fold=1.,
        closest_approach_contact_rate_s=None,
        rate_reason='Classical short-record joint approaches do not identify a normalized quantum state or rare-contact tail.')


def main():
    manifest = json.loads((HERE/'trajectories.json').read_text())
    records = {}
    arrays = {}
    for name in manifest['cases']:
        item = json.loads((HERE/name).read_text())
        label = item['summary']['label']
        records[label] = item
        arrays[label] = np.load(HERE/item['arrays_file'])
    contrasts = [('NiH_cold_h0', 'NiH_cold_h0.08', 'NiH_H_motion_transplanted_to_DD'),
                 ('PdD_random_x1_h0', 'PdD_random_x1_h0.08', 'PdD_cold_x1_patched'),
                 ('PdD_thermal_x0.875_h0', 'PdD_thermal_x0.875_h0.08', 'PdD_300K_x0.875_unlocked')]
    beam_rows = []
    for off, on, case in contrasts:
        for direction_name, direction in (('parallel_x',0),('transverse_z',2)):
            for reaction in ('eq7', 'proton'):
                results = {}
                for label in (off,on):
                    velocity = arrays[label]['velocities_m_s']
                    results[label] = {ch:project_motion(velocity,ch,reaction,direction=direction) for ch in ('p','n')}
                for channel in ('p','n','total'):
                    ys = {}
                    for label in (off,on):
                        channels = ('p','n') if channel=='total' else (channel,)
                        ys[label] = sum(results[label][ch]['physical_yield_per_C'] for ch in channels)
                    for i, energy in enumerate(LAB):
                        beam_rows.append(dict(case=case,reaction=reaction,channel=channel,
                            E_lab_keV=float(energy), E_cm_keV=float(energy/2),
                            unpumped_physical_events_per_C=float(ys[off][i]),
                            pumped_physical_events_per_C=float(ys[on][i]),
                            pumped_over_unpumped=float(ys[on][i]/ys[off][i]),
                            normalized_R_factor=float((ys[on][i]/ys[on][-1])/(ys[off][i]/ys[off][-1])),
                            screening_eV=40., beam_direction=direction_name,
                            species='D on D; NiH velocities are explicitly transplanted, not NiD trajectories'))
    save_csv('beam.csv',beam_rows)

    budget = 1e-9*17.6e6
    power_rows = []
    pair_rows = []
    for label, item in records.items():
        s = item['statistics']; c = item['summary']['configuration']
        source = ROOT/item['parameters_file']
        p = json.loads(source.read_text())['materials'][c['material']]
        A = p['K4_100_eV_A4']; B=(3*p['K4_111_eV_A4']-A)/2
        drag = np.array(s['drag_eV_per_s_per_occupied_site_blocks'])
        pump = np.array(s['explicit_pump_work_eV_per_s_per_occupied_site_blocks'])
        mass = p['mass_eVc2']/(2.99792458e18)**2
        omega = s['omega_rad_ps']*1e12
        coherent_drag = s['gamma_per_s']*mass*omega**2*item['summary']['driven_demodulated_rms_A']**2/2
        power_rows.append(dict(case=label,temperature_K=c['temperature_K'],
            occupied_site_drag_eV_s=float(drag.mean()),
            explicit_pump_work_eV_s=float(pump.mean()),
            thermal_injection_eV_s=s['classical_thermal_injection_eV_per_s_per_occupied_site'],
            ratio4_omega_over_gamma=s['ratio4_omega_over_gamma'],
            ratio4_frequency_over_gamma=s['ratio4_omega_over_gamma']/(2*np.pi),
            required_ratio4_gross_drag_DT_budget=s['ratio4_omega_over_gamma']*drag.mean()/budget,
            required_ratio4_gross_drag_DD_budget=s['ratio4_omega_over_gamma']*drag.mean()/(1e-9*3.65e6),
            harmonic_Omega_drag_eV_s=coherent_drag,
            required_ratio4_harmonic_Omega_DT_budget=s['ratio4_omega_over_gamma']*coherent_drag/budget,
            target_DT_budget_eV_s=budget,
            ratio1_Jsum_over_Kii=sum(row['J_sum_eV_A2'] for row in p['shells'])/p['Kii_eV_A2'],
            ratio2_K4_100_variance_over_Kii=A*p['zero_point_axis_variance_A2']/p['Kii_eV_A2'],
            ratio2_K4_111_variance_over_Kii=p['K4_111_eV_A4']*p['zero_point_axis_variance_A2']/p['Kii_eV_A2'],
            ratio3_hbarOmega_over_kBT='infinite' if c['temperature_K']==0 else (6.582119569e-16*s['omega_rad_ps']*1e12)/(8.617333262e-5*c['temperature_K']),
            ratio5_h_over_linear_Mathieu_threshold=c['pump_h']*s['ratio4_omega_over_gamma']/2,
            ratio6_vacancy_fraction=1-item['summary']['actual_loading'],
            ratio6_vacancy_shift_over_coupling='spatial law, not one scalar; see native loading_response',
            conditional_ratio7_lambda_over_nn=.541/(p['a_A']/np.sqrt(2)),
            velocity_squared_attenuation=s['harmonic_velocity_squared_attenuation']))
        pair_rows.append(pair_summary(item,arrays[label],p))
    save_csv('power.csv',power_rows); save_csv('pairs.csv',pair_rows)

    # Proper conditional pair potential, not a rate at instantaneous r_min.
    path=ROOT/'g6_pairlaw/script/pairlaw.py'
    spec=importlib.util.spec_from_file_location('conditional_pairlaw',path)
    pairlaw=importlib.util.module_from_spec(spec);spec.loader.exec_module(pairlaw)
    conditional=[]
    for material, parameter_file in (('NiH','g3_lattice/rerun/part2/updated/ni_epoch/parameters.json'),
                                     ('PdD','g3_lattice/rerun/part2/updated/parameters.json')):
        p=json.loads((ROOT/parameter_file).read_text())['materials'][material]
        center=p['a_A']/np.sqrt(2)
        stiffness=p['Kii_eV_A2']/2
        for isotope in ('DD','DT'):
            value=pairlaw.harmonic(.541,center,stiffness,isotope)
            conditional.append(dict(material_geometry=material,isotope=isotope,
                screening_A=.541,external_center_A=center,stiffness_eV_A2=stiffness,
                achieved_mean_A=value['mean_A'],log10_rate_s=value['log10_rate_s'],
                log10_p_rate_s=value['log10_rate_s']+np.log10(55.576/(55.576+53.701)) if isotope=='DD' else None,
                log10_n_rate_s=value['log10_rate_s']+np.log10(53.701/(55.576+53.701)) if isotope=='DD' else None,
                pumped_over_unpumped=1.,ratio7_lambda_over_actual_nn=.541/center,
                resonance_rate_s=None,resonance_reason='No identified bound-state matrix element; positive-energy DD reaction_S cannot be evaluated at a bound energy.',
                preparation='Conditional normalized radial ground state; identical potential in both pump labels. Not inferred from lattice trajectories.',
                model='C exp(-r/lambda)/r + k(r-r_center)^2/2; k=Kii/2 for two independent fixed harmonic host wells; same .541A screening scenario for both geometries.'))
    save_csv('conditional_pair.csv',conditional)

    # One decisive consumer check after all producers/edits settle: the
    # zero-velocity adapter must recover the actual native target, and the
    # actual pump contrast must reach and change that very consumer.
    target=reference()
    native=forward_yields(LAB/2,target.profiles,Ue_keV=.04,
                          stopping_fn=target.stopping_fn,reaction='proton')
    resting=project_motion(np.zeros((1,3)),'p','proton')
    rest_error=float(np.max(abs(resting['physical_yield_per_C']/native['physical_yield_per_C']-1)))
    assert rest_error<1e-12, 'Zero target velocity changed native consumer'
    ni=next(row for row in beam_rows if row['case']=='NiH_H_motion_transplanted_to_DD' and row['reaction']=='proton' and row['channel']=='p' and row['E_lab_keV']==.5)
    assert abs(ni['pumped_over_unpumped']-1)>1e-7, 'Actual pump velocity never reached the target consumer'
    # Finite observed temporal blocks, not independent-site confidence bars.
    blocks=[]; yield_blocks=[]
    off,on,_=contrasts[0]
    for block in range(8):
        result=[]; yields=[]
        for label in (off,on):
            data=arrays[label]
            velocity=data['velocities_m_s'][data['velocity_block']==block]
            value=folded_cross_section(np.array([.5]),velocity,.04,'p','proton',direction=0)[0]
            result.append(value)
            yields.append(project_motion(velocity,'p','proton',direction=0)['physical_yield_per_C'][0])
        blocks.append(float(result[1]/result[0]))
        yield_blocks.append(float(yields[1]/yields[0]))
    evidence=dict(native_rest_consumer_relative_error=rest_error,
        named_contrast='NiH cold .08 vs0 pump, H velocities explicitly transplanted into DD scattering on unchanged transported PdD target',
        named_contrast_E_lab_keV=.5,named_contrast_pumped_over_unpumped=ni['pumped_over_unpumped'],
        named_contrast_projectile_axis='parallel_x',
        consumer_matched_temporal_block_yield_ratios=yield_blocks,
        entrance_cross_section_matched_temporal_block_ratios=blocks,
        block_interpretation='Eight contiguous finite trajectory blocks; spread is a temporal diagnostic, not iid error bars or a rare-tail enclosure.',
        pairlaw_script_SHA256=hashlib.sha256(path.read_bytes()).hexdigest(),
        nuclear_contact_status='Only conditional quantum potential rates are quantified; actual pumped/unpumped lattice pair rate remains unidentified.',
        pump_off_result='At zero bath temperature, dH/dt=-m gamma sum(v^2)<=0 when pump is off. Positive drag cannot sustain nonzero periodic motion; no dynamic pump-off lock is provided by this model. Static displacement is not oscillatory lock.',
        power_interpretation='No division by domain N. Thermal injection belongs to all3 velocity components. Finite-difference drag at300K is not an exact OU-work estimate; compare matched pump increments and explicit Hamiltonian pump work.',
        largest_omissions=['Displacement-induced host electronic/nonlinear collision-frequency response; awaits collision-scope paired frozen snapshot.',
                           'Prepared many-body quantum pair density/contact, isotope-specific pair state and bound-state threshold-resonance matrix element.',
                           'Pump delivery/efficiency and backreaction, long-time thermal/damage/transport evolution and high-barrier rare tails beyond short classical records.'])
    (HERE/'summary.json').write_text(json.dumps(evidence,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(consumer=evidence['named_contrast_pumped_over_unpumped'],
                         rest_error=rest_error,power=power_rows,pairs=pair_rows,conditional_pair=conditional),indent=2),flush=True)


if __name__=='__main__':
    main()
