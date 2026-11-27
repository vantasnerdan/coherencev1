"""Analyze signed PdD displacements, fixed-volume excitation, and numerical sensitivity."""
import numpy as np


def signed_pair(zero, minus, plus, step, axis):
    axis = np.asarray(axis)
    force = np.array([minus['forces_eV_A'][1], plus['forces_eV_A'][1]])
    fm, fp = force @ axis
    energies = np.array([minus['energy_free_eV'], plus['energy_free_eV']])
    even_energy = (energies[0] + energies[1]) * 0.5 - zero['energy_free_eV']
    odd_restoring = (fm - fp) * 0.5
    return dict(step_A=step, raw_D_parallel_forces_eV_A=[float(fm), float(fp)],
                even_free_energy_change_eV=float(even_energy), odd_restoring_force_eV_A=float(odd_restoring),
                K_force_eV_A2=float(odd_restoring / step), K_free_energy_eV_A2=float(2 * even_energy / step**2),
                inversion_force_residual_eV_A=float(np.linalg.norm(force.sum(axis=0))),
                inversion_energy_residual_eV=float(abs(energies[1] - energies[0])),
                total_force_residual_eV_A=max((float(np.linalg.norm(np.asarray(q['forces_eV_A']).sum(axis=0)))
                                               for q in (zero, minus, plus))))


def quartic(pairs, k):
    """Jointly expose the two raw force-implied K4 values and energy holdout."""
    steps = np.array([p['step_A'] for p in pairs])
    force = np.array([p['odd_restoring_force_eV_A'] for p in pairs])
    energy = np.array([p['even_free_energy_change_eV'] for p in pairs])
    design = steps**3
    y = force - k * steps
    k4 = float(np.dot(design, y) / np.dot(design, design))
    fe = y - k4 * design
    ee = energy - (k * steps**2 / 2 + k4 * steps**4 / 4)
    return dict(K4_force_fit_eV_A4=k4, K4_per_amplitude_force_eV_A4=(y / design).tolist(),
                K4_per_amplitude_free_energy_eV_A4=(4 * (energy - k * steps**2 / 2) / steps**4).tolist(),
                force_residual_eV_A=fe.tolist(), free_energy_holdout_residual_eV=ee.tolist(),
                max_force_residual_eV_A=float(abs(fe).max()), max_free_energy_holdout_residual_eV=float(abs(ee).max()))


def thermal_delta(row, cold):
    du = row['zero']['energy_internal_eV'] - cold['zero']['energy_internal_eV']
    df = row['zero']['energy_free_eV'] - cold['zero']['energy_free_eV']
    dp = row['zero']['pressure_GPa'] - cold['zero']['pressure_GPa']
    k0 = cold['harmonic']['K_force_eV_A2']
    dk = row['harmonic']['K_force_eV_A2'] - k0
    dk_energy = row['harmonic']['K_free_energy_eV_A2'] - cold['harmonic']['K_free_energy_eV_A2']
    return dict(electronic_excitation_eV_per_D_site=du, electronic_free_energy_change_eV_per_D_site=df,
                pressure_change_GPa=dp, K_change_eV_A2=dk, K_free_energy_change_eV_A2=dk_energy,
                fractional_K_change=dk/k0,
                fractional_K_secant_per_eV=(dk/(k0 * du) if du > 0 else None),
                centered_D_force_change_eV_A=(np.array(row['zero']['forces_eV_A'][1]) - cold['zero']['forces_eV_A'][1]).tolist())


def mesh_check(a, b):
    ka = a['contrast']['K_change_eV_A2']
    kb = b['contrast']['K_change_eV_A2']
    err = abs(ka-kb)
    allowed = max(2e-4, 0.1 * abs(kb))
    fd = abs(kb - b['contrast']['K_free_energy_change_eV_A2'])
    return dict(delta_K_mesh_sensitivity_eV_A2=err, delta_K_force_energy_difference_eV_A2=fd,
                delta_K_tolerance_eV_A2=allowed, delta_K_mesh_agrees=bool(err <= allowed),
                delta_K_force_energy_agrees=bool(fd <= allowed),
                delta_K_sign_larger_than_sensitivities=bool(np.sign(ka)==np.sign(kb) and abs(kb)>max(err,fd)))
