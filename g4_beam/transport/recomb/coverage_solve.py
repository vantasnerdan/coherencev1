"""Solute/trap steady shooting with a source-bound frozen run-end vacancy field."""
import math
import numpy as np
from scipy.optimize import brentq
from coverage_boundary import coverage_from_mu, local_boundary, supplied_sheet
from coverage_inputs import make_model
from coverage_steady import loaded_state_evidence


def _traps(model, x, vacancy):
    q = (model.mu(np.asarray(x)) - model.iso.mu0() + model.mat['Eb']) / model.kT
    occupancy = 1 / (1 + np.exp(np.clip(-q, -200, 200)))
    # Match the unmodified Model's rhs projection; retain raw saved v in state.
    return model.mat['mtrap'] * np.maximum(vacancy, 0.0) * occupancy


def _shoot(model, xlast, sheet, capture_s0, omit_overflow, vacancy, minimum_x):
    """Anchor the stiff finite-current face, then sweep toward the sheet."""
    x = np.empty(model.N)
    x[-1] = xlast
    jd = model.D_of_x(x[-1:])[0]*model.nM*xlast*(model.mu_back-float(model.mu(x[-1:])[0]))/model.kT/(.5*model.h[-1])
    q = jd*model.Jmax_back/(abs(jd)+model.Jmax_back)
    for i in range(model.N-2, -1, -1):
        right = x[i+1]
        if not omit_overflow:
            trap = float(_traps(model, [right], vacancy[i+1])[0])
            q -= model.nM*model.h[i+1]*max(right+trap-model.mat['x_lmm'], 0)
        muright = float(model.mu(np.array([right]))[0])
        logright = math.log(right/(model.ths-right))
        def face(left):
            delta = logright-math.log(left/(model.ths-left))
            mean = (right-left)/delta if abs(delta) >= 1e-10 else .5*(left+right)*(1-.5*(left+right)/model.ths)
            diffusion = math.sqrt(float(model.D_of_x(np.array([left]))[0]) * float(model.D_of_x(np.array([right]))[0]))
            return diffusion*model.nM*mean*(muright-float(model.mu(np.array([left]))[0]))/model.kT/model.dzf[i]-q
        if face(minimum_x) < 0:
            return 1.0+abs(q)/model.Jmax_back, None
        if face(model.c['x_back']) > 0:
            return -1.0-abs(q)/model.Jmax_back, None
        x[i] = brentq(face, minimum_x, model.c['x_back'], xtol=5e-15)
    xt = _traps(model, x, vacancy)
    if not omit_overflow:
        q -= model.nM*model.h[0]*max(x[0]+xt[0]-model.mat['x_lmm'], 0)
    theta = coverage_from_mu(float(model.mu(x[:1])[0]), sheet, model.kT)
    boundary = local_boundary(model, theta, sheet=sheet, capture_s0=capture_s0)
    state = np.concatenate([x, xt, vacancy, [0.]])
    return (q-boundary['D_atom_net_outflux_m2_s'])/model.Jmax_back, (state, boundary)


def solve_loaded(source, row, *, thermodynamics, capture_s0, vacancies,
                 omit_overflow=False, grid=None, root_xtol=1e-18):
    if not row['loaded']:
        raise ValueError('Bulk-fed loaded row required')
    if (vacancies.provenance['case']['condition'] != 'PdD'
            or vacancies.provenance['case']['T_K'] != row['T_K']):
        raise ValueError('Loaded foil requires its matching run-end vacancy profile')
    sheet = supplied_sheet(thermodynamics, row['T_K'], capture_s0)
    model = make_model(source, row, s_effective=1., grid=grid)
    vacancy = vacancies.on_grid(model.e)
    xgas = model.iso.x_of_mu(.5*model.kT*math.log(model.p_gas))
    # Back-face depletion is tiny at low release flux. Anchoring it locally
    # avoids amplifying accumulated face-roundoff into the back-source flux.
    back = float(row['x_back'])
    f = lambda depletion: _shoot(model, back-depletion, sheet, capture_s0,
                                 omit_overflow, vacancy, xgas)[0]
    depletion, result = brentq(f, 0., back-xgas*(1+1e-12),
                              xtol=root_xtol, full_output=True)
    residual, payload = _shoot(model, back-depletion, sheet, capture_s0,
                               omit_overflow, vacancy, xgas)
    if payload is None:
        raise ArithmeticError('Steady shooting did not reach the complete foil')
    state, boundary = payload
    local = make_model(source, row, s_effective=boundary['s_effective'], grid=grid)
    evidence = loaded_state_evidence(local, state, sheet=sheet,
                                    capture_s0=capture_s0, frozen_v=vacancy)
    n = local.N
    steady_rhs = np.asarray(evidence['original_rhs']).copy()
    if omit_overflow:
        steady_rhs[:n] += np.maximum(state[:n]+state[n:2*n]-local.mat['x_lmm'], 0)
    scale = max(abs(evidence['front_atom_outflux_m2_s']), abs(evidence['back_atom_influx_m2_s']), 1)
    balance = abs(evidence['back_atom_influx_m2_s']-evidence['front_atom_outflux_m2_s'])/scale
    if not omit_overflow:
        balance = evidence['mass_balance_relative_residual']
    return dict(process_model='Bulk-fed solute/trap/sheet steady, saved run-end v(z) fixed',
        vacancy_provenance=vacancies.provenance,
        omit_overflow=omit_overflow, theta_candidate_ML=boundary['theta_ML'],
        s_effective=boundary['s_effective'], state=state.tolist(), grid_cells=n,
        grid_both=bool(local.c['grid_both']),
        first_cell_width_m=float(local.h[0]), last_cell_width_m=float(local.h[-1]),
        near_surface_depth_m=1e-7, first_cell_x_plus_xt=float(state[0]+state[n]),
        max_near_surface_x_plus_xt=float(np.max((state[:n]+state[n:2*n])[local.zc<=1e-7])),
        near_surface_reaches_095=bool(np.any((state[:n]+state[n:2*n])[local.zc<=1e-7]>=.95)),
        shooting_scaled_residual=float(residual), root_iterations=result.iterations,
        shooting_coordinate='back_face_reservoir_depletion',
        shooting_root_xtol=root_xtol, back_solute_depletion=depletion,
        front_solute_depletion=back-float(state[0]),
        steady_max_solute_rhs_per_s=float(np.max(np.abs(steady_rhs[:n]))),
        steady_max_cell_flux_relative_residual=float(np.max(np.abs(
            local.nM*local.h*(steady_rhs[:n]+steady_rhs[n:2*n]))))/scale,
        steady_boundary_mass_balance_relative_residual=balance,
        evidence=evidence, beam_stimulated_desorption=False)


def converge_loaded(source, row, *, thermodynamics, capture_s0, vacancies,
                    omit_overflow=False, coverage_target_ML=1e-4):
    if not 0 < coverage_target_ML <= 1e-4:
        raise ValueError('Positive coverage tolerance <=1e-4ML required')
    original = solve_loaded(source, row, thermodynamics=thermodynamics,
        capture_s0=capture_s0, vacancies=vacancies, omit_overflow=omit_overflow)
    previous, intermediate = original, None
    history = [dict(grid_cells=original['grid_cells'],
                    theta_ML=original['theta_candidate_ML'],
                    s_effective=original['s_effective'], refinement=1)]
    for factor in (2, 4):
        grid = dict(h0=.25e-9/factor, g1=1.05**(1/factor),
                    hmid=3e-9/factor, g2=1.13**(1/factor))
        refined = solve_loaded(source, row, thermodynamics=thermodynamics,
            capture_s0=capture_s0, vacancies=vacancies, omit_overflow=omit_overflow,
            grid=grid)
        delta = abs(previous['theta_candidate_ML']-refined['theta_candidate_ML'])
        s_delta = abs(previous['s_effective']-refined['s_effective']) / max(
            abs(previous['s_effective']), abs(refined['s_effective']), 1e-300)
        previous_rate = previous['evidence']['boundary']['D2_desorption_events_m2_s']
        refined_rate = refined['evidence']['boundary']['D2_desorption_events_m2_s']
        rate_delta = abs(previous_rate-refined_rate) / max(abs(previous_rate), abs(refined_rate), 1e-300)
        invariants = (refined['steady_max_cell_flux_relative_residual'] <= 1e-5
            and refined['steady_boundary_mass_balance_relative_residual'] <= 1e-5
            and refined['evidence']['max_trap_rhs_D_per_Pd_s'] <= 1e-10)
        converged = (delta <= coverage_target_ML and s_delta <= .05
                     and rate_delta <= .05 and invariants)
        history.append(dict(grid_cells=refined['grid_cells'],
            theta_ML=refined['theta_candidate_ML'], refinement=factor,
            s_effective=refined['s_effective'], effective_s_relative_change=s_delta,
            recombination_rate_relative_change=rate_delta,
            change_ML=delta, invariants_passed=invariants))
        if converged or not invariants:
            break
        if factor == 2:
            previous, intermediate = refined, refined
    return dict(status='converged_frozen_v_steady' if converged else 'not_converged',
        theta_ML=refined['theta_candidate_ML'] if converged else None,
        s_effective=refined['s_effective'] if converged else None,
        coverage_grid_change_ML=delta, coverage_target_ML=coverage_target_ML,
        effective_s_relative_grid_change=s_delta, effective_s_relative_target=.05,
        recombination_rate_relative_grid_change=rate_delta, recombination_rate_relative_target=.05,
        coarse=original, intermediate=intermediate, refined=refined,
        grid_refinement_history=history, maximum_refinement=4)
