"""Per-metal normalisation and conditional, non-causal screening/d-band slopes."""
import itertools
import numpy as np
from scipy.stats import linregress
from screen_io import require


def scalar(row):
    dos = row['dos']
    metals = sum(s != 'H' for s in row['geometry']['symbols'])
    d_atoms = row['geometry']['symbols'].count('H')
    result = dict(material_id=row['material_id'], metal_atoms_per_cell=metals,
        D_per_metal=d_atoms/metals,
        N_EF_tetrahedron_states_per_eV_per_cell=dos['n_EF_tetrahedron_states_per_eV_per_cell'],
        N_EF_tetrahedron_states_per_eV_per_metal=dos['n_EF_tetrahedron_states_per_eV_per_cell']/metals,
        occupied_metal_d_centre_relative_EF_eV=dos['occupied_d_centre_relative_EF_eV'],
        occupied_metal_d_rms_width_eV=dos['occupied_d_rms_width_eV'],
        by_metal_element=dos['by_metal_element'],
        D_projected_N_EF_states_per_eV_per_D=dos.get('D_projected_n_EF_tetrahedron_states_per_eV_per_D'),
        ground_state_settings=row['ground_state_settings'], geometry=row['geometry'],
        occupied_DOS_settings=dos['settings'], nspins=dos['nspins'],
        spin_degeneracy=dos['spin_degeneracy'], model=row['model'])
    if 'screening' in row:
        result['screening'] = row['screening']
    return result


def slope(rows, allowance=.1, coordinate='Pd'):
    require(len(rows) >= 2, 'A slope needs two or more independent hosts')
    require(coordinate in ['Pd','all_metal'],'Unknown d-band coordinate')
    x = np.array([(r['dos']['by_metal_element']['Pd'] if coordinate=='Pd' else r['dos'])[
                  'occupied_d_centre_relative_EF_eV'] for r in rows])
    y = np.array([r['screening']['Ue_static_Yukawa_proxy_eV'] for r in rows])
    require(len(np.unique(x)) == len(x) and np.all(y > 0), 'Degenerate slope inputs')
    weights = (x-x.mean())/np.sum((x-x.mean())**2)
    value = float(weights@y)
    oracle = float(linregress(x, y).slope)
    require(np.isclose(value, oracle, rtol=1e-12, atol=1e-12), 'OLS independent route mismatch')
    # This is an interval under a stated per-case lambda allowance, not a fitted CI.
    lower_U, upper_U = y/(1+allowance), y/(1-allowance)
    lower = float(np.sum(weights*np.where(weights >= 0, lower_U, upper_U)))
    upper = float(np.sum(weights*np.where(weights >= 0, upper_U, lower_U)))
    vertices = [float(weights@np.array([upper_U[i] if bit else lower_U[i] for i, bit in enumerate(bits)]))
                for bits in itertools.product([0, 1], repeat=len(rows))]
    require(np.isclose(min(vertices), lower) and np.isclose(max(vertices), upper),
            'Interval fails independent corner enumeration')
    return dict(materials=[r['material_id'] for r in rows],
        d_band_coordinate=coordinate,occupied_d_centres_relative_EF_eV=x.tolist(),
        Ue_Yukawa_proxy_eV=y.tolist(),
        slope_eV_per_eV=value, conditional_lambda_allowance=allowance,
        conditional_slope_interval_eV_per_eV=[lower, upper],
        positive_sign_resolved_under_this_allowance=bool(lower > 0),
        nonpositive_sign_resolved_under_this_allowance=bool(upper <= 0),
        scope='OLS association conditional on the explicitly named occupied PAW d-centre coordinate and independent per-host +/-lambda allowance. Not statistical confidence, not uncertainty in the d-centre, not a universal or causal slope.')


def report(rows,require_complete=True):
    byid = {r['material_id']:r for r in rows}
    require(len(byid)==len(rows) and 0<len(rows)<=15,'Expected distinct commissioned materials')
    require(not require_complete or len(rows)==15,'Expected all fifteen Tier1 materials')
    rpa = [r for r in rows if 'screening' in r]
    require(len(rpa)<=9 and (not require_complete or len(rpa)==9),'Expected nine direct-response cases')
    require(all(r['screening']['eligible_for_Tier1'] for r in rpa), 'A response q-window is unresolved')
    strain = [byid['Pd_hydro_'+s+'pct'] for s in ['-2', '+0', '+2']]
    alloys = [byid['Pd3Rh_L12'], byid['Pd3Ag_L12']]
    joint = strain+alloys
    combined = slope(joint)
    leave_one = [slope([r for r in joint if r is not excluded]) for excluded in joint]
    return dict(schema='coherence-screen-tier1-scalars-v1',
        tier1_complete=(len(rows)==15 and len(rpa)==9),
        completed_DOS_cases=len(rows),completed_RPA_cases=len(rpa),
        cases=[scalar(r) for r in sorted(rows, key=lambda r:r['material_id'])],
        headline_static_Yukawa_proxy_slope=combined,
        all_metal_centroid_definition_sensitivity=slope(joint,coordinate='all_metal'),
        strain_only_static_proxy_slope=slope(strain),
        alloy_paths={r['material_id']:slope([byid['Pd_hydro_+0pct'], r]) for r in alloys},
        pooled_leave_one_out_central_slope_range_eV_per_eV=[
            min(s['slope_eV_per_eV'] for s in leave_one),max(s['slope_eV_per_eV'] for s in leave_one)],
        projection_scope='PAW-projector occupied[-20,0]eV d moments, Gaussian0.1eV. Centroids are not all-valence DOS centres. N(EF) tetrahedron, per-cell and per-metal both explicit.',
        response_scope='Valence-PAW PBE static low-q RPA, actual relaxed/fixed declared geometry; lambda from q^2(epsilon_M-1) intercept. Ue=C/lambda and force factors are explicitly Yukawa-model outputs, NOT near-origin full-RPA or all-electron screening potentials.',
        hypothesis_scope='Pooled and path-specific associations must not be interchanged. A sign unresolved by the declared allowance is neither support nor refutation of a positive physical d-band hypothesis.')
