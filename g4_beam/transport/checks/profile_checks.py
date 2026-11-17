"""Seven-field SI profile windows and independent adaptive CSDA yields."""
import numpy as np
from scipy.constants import elementary_charge
from scipy.integrate import solve_ivp


def window(record, depth_m):
    metal = 0.0
    deuterium = 0.0
    previous = 0.0
    for layer in record['layers']:
        if len(layer) != 7:
            raise ValueError('Expected seven-field SI layers, including metal density')
        top, bottom, host, density_D, Z, mass, density_metal = layer
        if top < previous-1e-15 or bottom <= top or min(host, density_D, density_metal) < 0:
            raise ValueError('Invalid physical layer ordering or densities')
        previous = bottom
        width = max(0., min(bottom, depth_m)-max(top, 0.))
        metal += width*density_metal
        deuterium += width*density_D
    if metal <= 0:
        raise ValueError('Window has no metal inventory')
    return dict(D_per_metal=deuterium/metal,
                D_density_nm3=deuterium/depth_m*1e-27,
                D_column_m2=deuterium, metal_column_m2=metal,
                adsorbate_excluded=True)


def probability(record, E_lab_eV, stopping, sigma_mb):
    E = float(E_lab_eV)
    probability = record['nads_m2']*sigma_mb(E/2000, 0.)*1e-31
    cutoff = .03*E_lab_eV
    for layer in record['layers']:
        top, bottom, host, density_D, Z, mass, density_metal = layer
        if E <= cutoff:
            break
        width_nm = (bottom-top)*1e9
        if width_nm <= 0:
            raise ValueError('Nonpositive layer width')

        def derivative(z_nm, state):
            energy = max(state[0], 1.)
            return [-float(stopping(energy, host, Z, mass, density_D)),
                    density_D*1e-27*float(sigma_mb(energy/2000, 0.))*1e-13]

        def stop(z_nm, state):
            return state[0]-cutoff

        stop.terminal = True
        stop.direction = -1
        solution = solve_ivp(derivative, [0., width_nm], [E, probability],
                             rtol=1e-7, atol=[1e-6, 1e-80],
                             max_step=min(width_nm, 1.), events=stop)
        if not solution.success:
            raise ValueError('Independent CSDA integration failed: '+solution.message)
        E, probability = solution.y[:, -1]
    if not np.isfinite(probability) or not 0 <= probability < 1:
        raise ValueError('Invalid physical per-deuteron fusion probability')
    return float(probability)


def reactions_per_true_C(record, stopping, sigma_mb):
    E_lab = 2000*record['case']['E_cm_keV']
    return (.85*probability(record, E_lab, stopping, sigma_mb)
            +.30*probability(record, E_lab/2, stopping, sigma_mb))/elementary_charge


def check_profiles(records, stopping, sigma_mb):
    windows = []
    yields = {}
    max_error = 0.
    for record in records:
        case = record['case']
        key = (case['condition'], case['E_cm_keV'])
        row = dict(case=case, windows={})
        for nm in [5, 20, 70]:
            result = window(record, nm*1e-9)
            expected = record['end']['DM'+str(nm)]
            error = abs(result['D_per_metal']-expected)/max(abs(expected), 1e-30)
            max_error = max(max_error, error)
            if error > 1e-10:
                raise ValueError('Saved window differs from independent SI integral')
            row['windows'][str(nm)+'nm'] = result
        value = reactions_per_true_C(record, stopping, sigma_mb)
        yields[key] = value
        row['bare_DDp_reactions_per_true_incident_C'] = value
        windows.append(row)
    ti = [dict(E_cm_keV=energy,
               Ti_reactions_per_true_C=yields[('Ti', energy)],
               TiD_reactions_per_true_C=yields[('TiD', energy)],
               TiD_over_Ti=yields[('TiD', energy)]/yields[('Ti', energy)])
          for energy in sorted(energy for condition, energy in yields if condition == 'Ti')]
    return dict(profile_results=windows, Ti_yield_comparison=ti,
                max_saved_window_relative_error=max_error,
                yield_units='physical reactions per TRUE incident charged C; no detector/cup correction',
                method='independent adaptive RK45 CSDA; same declared stopping/sigma model, not independent nuclear data')
