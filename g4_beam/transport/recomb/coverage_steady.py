"""Audit original solute/trap equations while the saved vacancy profile is fixed."""
import numpy as np
from coverage_boundary import coverage_from_mu, local_boundary


def loaded_state_evidence(model, state, *, sheet, capture_s0, frozen_v):
    y = np.asarray(state, dtype=float)
    n = model.N
    if y.shape != (3 * n + 1,) or not np.all(np.isfinite(y)):
        raise ValueError('Finite original Model 3*N+1 state required')
    x, xt, v = y[:n], y[n:2*n], y[2*n:3*n]
    if (np.any(x <= 0) or np.any(x >= model.ths) or np.any(xt < 0)
            or not np.array_equal(v, frozen_v) or y[-1] != 0):
        raise ValueError('Unclipped solute/traps and the exact prescribed vacancy profile required')
    mu = model.mu(x)
    theta = coverage_from_mu(float(mu[0]), sheet, model.kT)
    boundary = local_boundary(model, theta, sheet=sheet, capture_s0=capture_s0)
    rhs = model.rhs(0.0, y, beam_on=False)
    jf = float(model.front_flux(float(mu[0]), 0.0)[0])
    jd = model.D_of_x(x[-1:])[0] * model.nM * x[-1] * (model.mu_back - mu[-1]) / model.kT / (0.5 * model.h[-1])
    jb = float(jd * model.Jmax_back / (abs(jd) + model.Jmax_back))
    lmm = np.maximum(x + xt - model.mat['x_lmm'], 0.0)
    sink = float(np.sum(model.nM * model.h * lmm))
    inventory_rate = float(np.sum(model.nM * model.h * (rhs[:n] + rhs[n:2*n])))
    scale = max(abs(jf), abs(jb), sink, 1.0)
    cell_rate = model.nM * model.h * (rhs[:n] + rhs[n:2*n])
    dmu_dx = (float(model.mu(np.array([x[0] * (1 + 1e-7)]))[0])
              - float(mu[0])) / (x[0] * 1e-7)
    theta_rate = sheet.dtheta_dmu(theta, model.kT) * dmu_dx * rhs[0]
    return dict(boundary=boundary, s_front_used=float(model.s_front),
        s_self_consistency_relative_error=abs(model.s_front-boundary['s_effective']) / max(boundary['s_effective'], 1e-300),
        front_atom_outflux_m2_s=jf, back_atom_influx_m2_s=jb,
        back_supply_limit_atoms_m2_s=float(model.Jmax_back),
        local_mixing_atom_loss_m2_s=sink, max_local_mixing_D_per_Pd_s=float(lmm.max()),
        inventory_atom_rate_m2_s=inventory_rate,
        telescoping_flux_identity_error_m2_s=inventory_rate-(jb-jf-sink),
        mass_balance_relative_residual=abs(jb-jf-sink)/scale,
        max_cell_inventory_flux_relative_residual=float(np.max(np.abs(cell_rate)))/scale,
        max_solute_rhs_D_per_Pd_s=float(np.max(np.abs(rhs[:n]))),
        max_trap_rhs_D_per_Pd_s=float(np.max(np.abs(rhs[n:2*n]))),
        original_vacancy_rhs_max_per_s=float(np.max(np.abs(rhs[2*n:3*n]))),
        vacancy_equation_enforced=False,
        vacancy_constraint='Saved run-end v(z) is held fixed, not thermally equilibrated',
        surface_coverage_rate_ML_s=float(theta_rate),
        vacancy_fraction_range=[float(v.min()), float(v.max())],
        theta_ML=theta, original_rhs=rhs.tolist())
