"""Read saved ensembles and exact pulse output; never rerun their producers."""
from pathlib import Path
import csv
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'g13_rogue/output'


def rows(name):
    with (OUT/name).open() as f:
        return list(csv.DictReader(f))


def table(name, data):
    with (OUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def main():
    survival=rows('survivals.csv')
    configs=rows('packet_geometry.csv')
    stats=json.loads((OUT/'statistics.json').read_text())
    duration=stats['window_fs']*1e-15
    ambient={(r['model'],r['threshold_eV']):r for r in survival if r['zone']=='-1'}
    with np.load(OUT/'maxima.npz') as arrays:
        maxima=(arrays['local_eV'],arrays['closing_eV'])
    output=[]
    for model in ('linear','nonlinear'):
        for energy,eta in dict.fromkeys((r['E_cm_keV'],r['eta']) for r in configs):
            for current in (.1,1.,10.):
                selected=[r for r in survival if r['model']==model and r['E_cm_keV']==energy and r['eta']==eta and r['zone']!='-1']
                thresholds=list(dict.fromkeys(r['threshold_eV'] for r in selected))
                for threshold in thresholds:
                    base=ambient[model,threshold]
                    response=np.zeros((2,stats['ensemble_phases'],2))
                    resolution=0.
                    for i,config in enumerate(configs):
                        if config['E_cm_keV']!=energy or config['eta']!=eta:
                            continue
                        weight=current*float(config['flux_nm2_s'])*float(config['area_stratum_nm2'])/32
                        resolution=max(resolution,weight/stats['ensemble_phases'])
                        for kind,values in enumerate(maxima):
                            response[:,:,kind]+=weight*((values[:,:,i+1]>float(threshold)).astype(float)-
                                                         (values[:,:,0]>float(threshold)).astype(float))
                    branch=int(model=='nonlinear')
                    local,closing=response[branch].mean(axis=0)
                    sem=response[branch].std(axis=0,ddof=1)/np.sqrt(stats['ensemble_phases'])
                    difference=response[1]-response[0]
                    comparison=difference.mean(axis=0)
                    comparison_sem=difference.std(axis=0,ddof=1)/np.sqrt(stats['ensemble_phases'])
                    output.append(dict(model=model,E_cm_keV=energy,eta=eta,current_scale=current,
                        threshold_eV=threshold,tagged_window_duration_s=duration,
                        ambient_window_local_max_survival=base['max_local_survival'],
                        ambient_window_closing_max_survival=base['max_closing_survival'],
                        single_packet_local_max_excess_events_per_D_s=local,
                        single_packet_closing_max_excess_events_per_D_s=closing,
                        local_excess_cluster_SE_per_D_s=float(sem[0]),closing_excess_cluster_SE_per_D_s=float(sem[1]),
                        paired_nonlinear_minus_linear_local_excess_per_D_s=float(comparison[0]),
                        paired_nonlinear_minus_linear_closing_excess_per_D_s=float(comparison[1]),
                        paired_comparison_local_cluster_SE_per_D_s=float(comparison_sem[0]),
                        paired_comparison_closing_cluster_SE_per_D_s=float(comparison_sem[1]),
                        largest_stratum_one_count_resolution_per_D_s=resolution,
                        all_sampled_responses_zero=bool(not np.any(response)),
                        semantics='arrival-centered window maxima excess; not a stationary arbitrary-window maximum CDF or individual upcrossing count'))
    table('current_maxima.csv',output)
    contact=json.loads((OUT/'contact.json').read_text())
    proton=5.81220367098898e7
    sizes=[]
    for pulse in contact['pulses']:
        for variant in pulse['variants']:
            for band in ('lower_band','upper_band'):
                b=variant[band]
                sizes.append(dict(join_A=variant['join_A'],include_flagged=variant['include_flagged'],band=band,
                    log10_proton_probability_per_pair_pulse=b['log10_probability_proton_per_pair_pulse'],
                    log10_either_branch_probability_per_pair_pulse=b['log10_probability_either_DD_branch_per_pair_pulse'],
                    log10_total_contact_exposure_A_minus3_s=b['log10_total_Xi_A_minus3_s'],
                    log10_peak_proton_hazard_per_pair_s=b['cusp_peak_log10_rho_A_minus3']+math.log10(proton),
                    log10_max_snapshot_R_lt_1fm=b['log10_max_instantaneous_probability_R_lt_1fm'],
                    rms_quantum_classical_R_mismatch_A=b['rms_R_mismatch_A'],
                    rms_quantum_classical_velocity_mismatch_A_fs=b['rms_velocity_mismatch_A_fs'],
                    closing_fraction_of_contact_exposure=b['closing_fraction_of_Xi'],
                    upper_minus_lower_log10_exposure=variant['upper_minus_lower_log10_Xi'],
                    excess_sign=b['signed_excess_Xi']['sign'],
                    log10_excess_contact_exposure=b['signed_excess_Xi']['log10_magnitude_A_minus3_s']))
    table('contact_sizes.csv',sizes)
    growth=rows('growth.csv')
    thresholds=[]
    for case in dict.fromkeys(r['case'] for r in growth):
        candidates=[r for r in growth if r['case']==case and r['carrier']=='top' and r['amplitude_lifetime_fs']=='8200.0' and r['threshold_action_quanta']]
        row=min(candidates,key=lambda r:float(r['threshold_action_quanta']))
        thresholds.append(dict(case=case,wavelength_nm=float(row['wavelength_nm']),
            threshold_action_quanta=float(row['threshold_action_quanta']),
            threshold_meV_per_D=float(row['threshold_meV_per_D']),
            threshold_displacement_A=float(row['threshold_displacement_A'])))
    table('minimum_thresholds.csv',thresholds)
    summary=dict(thresholds=thresholds,
        contact_log10_proton_chance_range=[min(r['log10_proton_probability_per_pair_pulse'] for r in sizes),max(r['log10_proton_probability_per_pair_pulse'] for r in sizes)],
        contact_log10_peak_proton_hazard_range=[min(r['log10_peak_proton_hazard_per_pair_s'] for r in sizes),max(r['log10_peak_proton_hazard_per_pair_s'] for r in sizes)],
        contact_band_shift_dex=[min(r['upper_minus_lower_log10_exposure'] for r in sizes),max(r['upper_minus_lower_log10_exposure'] for r in sizes)],
        quantum_classical_R_mismatch_A=[min(r['rms_quantum_classical_R_mismatch_A'] for r in sizes),max(r['rms_quantum_classical_R_mismatch_A'] for r in sizes)],
        physical_pair_preparation_probability=None,physical_reaction_rate=None,
        reaction_scope='pumped prepared seed, normalized spatial l=0; two-transition Kerr continuation crossed, no material population inferred',
        runtime_readout='python g13_rogue/readout.py',
        rate_semantics='current_maxima.csv counts signed Palm excess of one arrival-centered maximum per packet window, distinct from rates.csv multiple upward crossings and occupancy')
    (OUT/'sizes.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':
    main()
