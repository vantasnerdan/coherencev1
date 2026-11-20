"""Source-bound finite-domain LiT return readout; NumPy only, no Native imports."""
import json
import numpy as np
from screen_io import manifest,require,write_json
from substrate_framework.lithium_fuel_return import (
    MASS_D,MASS_T,break_even_barns,dt_eq7_sigma_mb,
    lit_icru_electronic,lit_inherited_sqrt_electronic,lit_zbl_nuclear,
    thick_target_return)


def tabulated_bragg(energy_lab_keV, bridge):
    """Interpolate literal source rows at equal velocity; never extrapolate."""
    t = np.asarray(energy_lab_keV, dtype=float)/MASS_D
    def atom(key):
        rows = bridge[key]['rows']
        x = np.array([float(r['energy']) for r in rows])
        y = np.array([float(r['stopping_power']) for r in rows])
        return x, y
    lithium_x, lithium_y = atom('Li')
    require(t.min() >= lithium_x[0] and t.max() <= lithium_x[-1], 'Li observed domain')
    total = np.interp(t, lithium_x, lithium_y)
    switch = bridge['hydrogen_switch_keV_per_u']
    for key, selected in [('H_low', t < switch), ('H_high', t >= switch)]:
        if selected.any():
            x, y = atom(key)
            require(t[selected].min() >= x[0] and t[selected].max() <= x[-1], 'H observed domain')
            total[selected] += np.interp(t[selected], x, y)
    return total


def finite_return_with_jumps(energy, total, screening_keV, jumps):
    """Integrate each continuous segment with explicit one-sided limits.

    A rejected/conditional source join must not acquire a grid-dependent
    smoothing interval merely because a trapezoidal consumer was used.
    """
    if not jumps:
        return thick_target_return(energy, total, screening_keV)
    lookup = {row['energy_lab_keV']: row for row in jumps
              if energy[0] < row['energy_lab_keV'] < energy[-1]}
    bounds = [float(energy[0]), *sorted(lookup), float(energy[-1])]
    gain = np.zeros_like(energy)
    integrated_gain_times_energy = 0.
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        selected = (energy >= lo) & (energy <= hi)
        e = energy[selected]
        require(e[0] == lo and e[-1] == hi, 'Explicit source boundary on energy grid')
        s = total[selected].copy()
        if lo in lookup:
            s[0] = lookup[lo]['right_total_per_formula']
        if hi in lookup:
            s[-1] = lookup[hi]['left_total_per_formula']
        g = thick_target_return(e, s, screening_keV)
        gain[selected] = g+integrated_gain_times_energy/e
        integrated_gain_times_energy = float(gain[selected][-1]*hi)
    return gain


def return_curves(stopping, output, bridge_path=None):
    native_e = np.asarray(stopping['energy_lab_keV'])
    native_s = np.asarray(stopping['electronic_cross_section_per_formula'])
    bridge = json.loads(bridge_path.read_text()) if bridge_path is not None else None
    source_boundaries = [bridge['hydrogen_switch_keV_per_u']*MASS_D] if bridge else []
    # Exact requested fixed cut .1keV, plus twofold sensitivities .05 and .2.
    energy = np.unique(np.r_[np.geomspace(.05, 200, 2401), .1, .2, .5, 1., 2.,
                             5., 8., 13., 25., native_e[-1], 50., 100., 200., source_boundaries])
    icru = lit_icru_electronic(energy)
    inherited = lit_inherited_sqrt_electronic(energy)
    nuclear = lit_zbl_nuclear(energy)
    below = energy <= native_e[-1]
    native_values = np.interp(np.log(energy[below]), np.log(native_e), native_s)
    replacements = {}
    nuclear_at_join = float(lit_zbl_nuclear(native_e[-1]))
    baselines = [('source_ICRU49', icru, float(lit_icru_electronic(native_e[-1])), []),
                 ('inherited_sqrt', inherited, float(lit_inherited_sqrt_electronic(native_e[-1])), [])]
    bridge_details = None
    if bridge_path is not None:
        table = icru.copy()
        high = energy > native_e[-1]
        if high.any():
            table[high] = tabulated_bragg(energy[high], bridge)
        overlap = float(tabulated_bragg(np.array([native_e[-1]]), bridge)[0])
        switch = bridge['hydrogen_switch_keV_per_u']
        h_at_switch = {}
        for key in ('H_low', 'H_high'):
            rows = bridge[key]['rows']
            h_at_switch[key] = float(np.interp(switch,
                [float(r['energy']) for r in rows], [float(r['stopping_power']) for r in rows]))
        switch_lab = switch*MASS_D
        li_rows = bridge['Li']['rows']
        li_switch = float(np.interp(switch, [float(r['energy']) for r in li_rows],
                                    [float(r['stopping_power']) for r in li_rows]))
        nuclear_switch = float(lit_zbl_nuclear(switch_lab))
        high_jumps = [dict(energy_lab_keV=switch_lab,
                          left_total_per_formula=li_switch+h_at_switch['H_low']+nuclear_switch,
                          right_total_per_formula=li_switch+h_at_switch['H_high']+nuclear_switch,
                          reason='Literal H-source handoff, no continuity forcing')] if switch_lab > native_e[-1] else []
        baselines.append(('compiled_Bragg_high_bridge', table, overlap, high_jumps))
        bridge_details = dict(source=manifest(bridge_path), equal_velocity_keV_per_u=float(native_e[-1]/MASS_D),
                              compiled_Bragg_overlap=overlap, native_overlap=float(native_s[-1]),
                              relative_mismatch=float(native_s[-1]/overlap-1),
                              source_handoff_keV_per_u=switch, H_per_atom_at_handoff=h_at_switch,
                              H_handoff_relative_difference=h_at_switch['H_high']/h_at_switch['H_low']-1,
                              applicability=bridge['policy'], availability=bridge['availability'],
                              uncertainty=bridge['errors'], normalization_caveat=bridge['Li_normalization_caveat'])
        bridge_details['within_30_percent'] = abs(bridge_details['relative_mismatch']) <= .3
        bridge_details['first_tabulated_Li_equal_velocity_keV_per_u'] = float(bridge['Li']['rows'][0]['energy'])
        bridge_details['diagnostic_25keV'] = 'Native loss query retained; below the first Li table datum at equal velocity, so no measured-data extrapolation.'
    jumps_by_model = {}
    for label, baseline, right_electronic, high_jumps in baselines:
        replacements[label+'_baseline'] = baseline
        jumps_by_model[label+'_baseline'] = list(high_jumps)
        if label == 'compiled_Bragg_high_bridge':
            jumps_by_model[label+'_baseline'].append(dict(energy_lab_keV=float(native_e[-1]),
                left_total_per_formula=float(lit_icru_electronic(native_e[-1]))+nuclear_at_join,
                right_total_per_formula=right_electronic+nuclear_at_join,
                reason='Conditional source-model to compiled-Bragg boundary'))
        replacement = baseline.copy()
        replacement[below] = native_values
        replacements[label+'_native_replacement'] = replacement
        jumps_by_model[label+'_native_replacement'] = [*high_jumps,
            dict(energy_lab_keV=float(native_e[-1]), left_total_per_formula=float(native_s[-1])+nuclear_at_join,
                 right_total_per_formula=right_electronic+nuclear_at_join,
                 reason='Conditional native-range boundary; no calibrated matching')]
    result = dict(energy_lab_keV=energy.tolist(), inherited_E_cut_keV=.1,
                  cutoff_sensitivity_keV=[.05, .1, .2],
                  inherited_Nuclear_ZBL_per_formula=nuclear.tolist(),
                  source_ICRU49_electronic_per_formula=icru.tolist(),
                  inherited_sqrt_electronic_per_formula=inherited.tolist(),
                  native_join_lab_keV=float(native_e[-1]),
                  native_join_relative_to_source_ICRU49=float(native_s[-1]/lit_icru_electronic(native_e[-1])-1),
                  join_policy='No matching factor or smoothing; explicit finite-domain substitution. Eligibility requires independent30percent overlap and numerical controls.',
                  below_Bosch_Hale_fit_domain_lab_keV=.5*(MASS_D+MASS_T)/MASS_T,
                  cross_section_caveat='Below .5keV centre-of-mass, Bosch-Hale is labelled model extrapolation, not measured cross section.',
                  isotope_caveat='Electronic LiD Born-Oppenheimer curve reused for LiT at the same calculated volume; no phonon/zero-point isotope correction.',
                  approximation='Source H+Li Bragg sum above native range; source_ICRU49 is full Geant4 table law. inherited_sqrt retains the old forward low-velocity-law extrapolation.',
                  curves={})
    result['compiled_measured_Bragg_bridge'] = bridge_details
    result['ICRU49_overlap_within_30_percent'] = abs(result['native_join_relative_to_source_ICRU49']) <= .3
    for name, electronic in replacements.items():
        total = electronic+nuclear
        entry = dict(electronic_per_formula=electronic.tolist(), total_per_formula=total.tolist(),
                     break_even_barns=break_even_barns(total).tolist(),
                     explicit_one_sided_boundaries=jumps_by_model[name], gains={})
        entry['applicability'] = (
            'Conditional replacement arithmetic. Above the native range, bridge admission and numerical controls are required; no physical G verdict is licensed by this curve alone.'
            if name.endswith('native_replacement') else
            'Finite-domain baseline model arithmetic with explicit inherited/source extrapolations; not measured LiT yield.')
        for cutoff in (.05, .1, .2):
            keep = energy >= cutoff
            gains = {str(u):finite_return_with_jumps(energy[keep], total[keep], u,
                                                    jumps_by_model[name]).tolist()
                     for u in (0., .84, 1.4)}
            entry['gains'][str(cutoff)] = dict(energy_lab_keV=energy[keep].tolist(), G=gains)
        result['curves'][name] = entry
    # Nuclear reference samples are original printed TableV (p621), not fit data.
    ecm = np.array([3., 5., 10., 20., 50., 100.])
    expected = np.array([.009808, .5383, 27.02, 407.7, 4219., 3427.])
    predicted = dt_eq7_sigma_mb(ecm*(MASS_D+MASS_T)/MASS_T)
    result['Bosch_Hale_printed_TableV_check'] = dict(Ecm_keV=ecm.tolist(), sigma_mb=predicted.tolist(),
                                                   printed_mb=expected.tolist(), relative_difference=(predicted/expected-1).tolist())
    write_json(output/'return-curves.json', result)
    return result
