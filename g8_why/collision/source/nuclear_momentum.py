"""Fourier-Hermite quadrature and canonical joint excitation/kinetic measure.

Finite-Hermite low-momentum results are numerical approximations to the declared
well. Their Gaussian ultraviolet behavior is never a physical Coulomb tail.
Only first/second kinetic moments are used for the reported numerical envelopes.
"""
import numpy as np

from quantum_dvr import (KINETIC_EV_A2, canonical_weights, kinetic_matrix,
                         quadrature_map, tensor_kinetic, tensor_transform)


def radial_shells(q, length_A):
    kinetic = KINETIC_EV_A2 / length_A**2 * (
        q[:, None, None]**2 + q[None, :, None]**2 + q[None, None, :]**2)
    rounded, inverse = np.unique(np.round(kinetic.ravel(), 14), return_inverse=True)
    shells = np.bincount(inverse, weights=kinetic.ravel()) / np.bincount(inverse)
    error = float(abs(kinetic.ravel() - shells[inverse]).max())
    return shells, inverse, error


def extended_moments(coefficients, length_A):
    """Continuum p^2 action with two extra Hermite levels, avoiding top-edge loss."""
    n = coefficients.shape[0]
    padded = np.pad(coefficients, ((0, 2), (0, 2), (0, 2), (0, 0)))
    action = tensor_kinetic(padded, kinetic_matrix(n + 2, length_A))
    mean = np.einsum('ijks,ijks->s', padded, action)
    second = np.einsum('ijks,ijks->s', action, action)
    return mean, second


def momentum_and_position(spectrum, weights, output):
    n = len(spectrum.nodes_A)
    order = n + 4
    q, gh_weights, transform = quadrature_map(n, order)
    fourier = transform * (-1j)**np.arange(n)[None, :]
    shells, inverse, grouping_error = radial_shells(q, spectrum.length_A)
    joint = np.empty((n**3, len(shells)))
    ground_shells = np.empty(len(shells))
    cubic = np.zeros((order, order, order))
    spatial = np.zeros_like(cubic)
    ground_cubic = np.empty_like(cubic)
    ground_spatial = np.empty_like(cubic)
    moment_error = np.zeros(2)
    max_state_norm_error = 0.0
    for start in range(0, n**3, 32):
        stop = min(start + 32, n**3)
        dvr = spectrum.vectors_dvr[:, start:stop].reshape(n, n, n, -1)
        coefficients = tensor_transform(dvr, spectrum.dvr_map.T)
        momentum = abs(tensor_transform(coefficients, fourier))**2
        position = tensor_transform(coefficients, transform)**2
        cubic += momentum @ weights[start:stop]
        spatial += position @ weights[start:stop]
        for j in range(stop - start):
            row = np.bincount(inverse, weights=momentum[..., j].ravel(), minlength=len(shells))
            joint[start + j] = weights[start + j] * row
            max_state_norm_error = max(max_state_norm_error, abs(float(row.sum()) - 1))
            if start + j == 0:
                ground_shells[:] = row
                ground_cubic[:] = momentum[..., j]
                ground_spatial[:] = position[..., j]
        first, second = extended_moments(coefficients, spectrum.length_A)
        rows = np.stack([np.bincount(inverse, weights=momentum[..., j].ravel(),
                                     minlength=len(shells)) for j in range(stop - start)])
        moment_error = np.maximum(moment_error,
            [abs(rows @ shells - first).max(), abs(rows @ shells**2 - second).max()])
    quadrature_measure = (gh_weights[:, None, None] * gh_weights[None, :, None]
                          * gh_weights[None, None, :] * np.exp(
                              q[:, None, None]**2 + q[None, :, None]**2 + q[None, None, :]**2))
    np.savez(output / 'joint-spectrum.npz', excitation_eV=spectrum.energies - spectrum.energies[0],
             kinetic_eV=shells, joint_probabilities=joint, canonical_weights=weights,
             ground_probabilities=ground_shells, omission_probability=0.0,
             shell_grouping_max_error_eV=grouping_error,
             shell_grouping_decimal_places=14, gauss_hermite_order=order)
    np.savez(output / 'momentum-position.npz', momentum_k_A_inverse=q / spectrum.length_A,
             momentum_cubic_probability=cubic, ground_momentum_cubic_probability=ground_cubic,
             momentum_density_per_A_inverse_cubed=cubic / quadrature_measure * spectrum.length_A**3,
             ground_momentum_density_per_A_inverse_cubed=ground_cubic / quadrature_measure * spectrum.length_A**3,
             momentum_volume_weights_A_inverse_cubed=quadrature_measure / spectrum.length_A**3,
             position_nodes_A=q * spectrum.length_A, position_cubic_probability=spatial,
             ground_position_cubic_probability=ground_spatial,
             radial_k_A_inverse=np.sqrt(shells / KINETIC_EV_A2),
             radial_probability=joint.sum(axis=0), ground_radial_probability=ground_shells,
             gh_nodes=q, gh_weights=gh_weights, weighted_hermite_map=transform,
             fourier_weighted_hermite_map=fourier)
    return dict(quadrature_order=order, max_state_normalization_error=max_state_norm_error,
                kinetic_first_identity_error_eV=float(moment_error[0]),
                kinetic_second_identity_error_eV2=float(moment_error[1]),
                shell_grouping_max_error_eV=grouping_error,
                shell_count=len(shells), cubic_normalization=float(cubic.sum()),
                spatial_normalization=float(spatial.sum()), omission_probability=0.0)


def consume_momentum(output, well):
    """Consume persisted joint and cubic measures; expose their actual scope."""
    with np.load(output / 'joint-spectrum.npz') as saved:
        joint = saved['joint_probabilities']
        kinetic = saved['kinetic_eV']
        radial = joint.sum(axis=0)
        ground = saved['ground_probabilities']
        state_weights_error = float(abs(joint.sum(axis=1) - saved['canonical_weights']).max())
        thermal = distribution_summary(radial, kinetic)
        ground_summary = distribution_summary(ground, kinetic)
        excitation = saved['excitation_eV']
        excitation_mean = float(joint.sum(axis=1) @ excitation)
        excitation_gt100 = float(joint[excitation > 100].sum())
    with np.load(output / 'momentum-position.npz') as saved:
        nodes = saved['position_nodes_A']
        r2 = (nodes[:, None, None]**2 + nodes[None, :, None]**2 + nodes[None, None, :]**2)
        position = saved['position_cubic_probability']
        ground_position = saved['ground_position_cubic_probability']
        localization = dict(mean_r2_A2=float((position * r2).sum()),
                            ground_mean_r2_A2=float((ground_position * r2).sum()))
        for radius in (0.65, 0.75, 0.90, 1.0):
            localization[f'probability_r_gt_{radius:g}A'] = float(position[r2 > radius**2].sum())
            localization[f'ground_probability_r_gt_{radius:g}A'] = float(ground_position[r2 > radius**2].sum())
        nearest = np.full_like(r2, np.inf)
        for nucleus in well.field.positions:
            delta = [(nodes + origin - atom + L / 2) % L - L / 2
                     for origin, atom, L in zip(well.origin, nucleus, well.field.lengths)]
            distance2 = delta[0][:, None, None]**2 + delta[1][None, :, None]**2 + delta[2][None, None, :]**2
            nearest = np.minimum(nearest, distance2)
        core = {f'within_{radius:g}A': dict(
            quadrature_node_count=int((nearest < radius**2).sum()),
            thermal_sampled_probability=float(position[nearest < radius**2].sum()),
            ground_sampled_probability=float(ground_position[nearest < radius**2].sum()))
            for radius in (0.10, 0.20, 0.30)}
        density_normalization = float((saved['momentum_density_per_A_inverse_cubed']
                                      * saved['momentum_volume_weights_A_inverse_cubed']).sum())
    return dict(joint_normalization=float(joint.sum()), state_weight_marginal_error=state_weights_error,
                cubic_density_normalization=density_normalization,
                thermal300K=thermal, ground0K=ground_summary, localization=localization,
                mean_excitation_eV=excitation_mean, probability_excitation_gt100eV=excitation_gt100,
                physical_point_core_diagnostic=core,
                core_scope='Indicator quadrature of the finite-Hermite state, not a certified core probability; zero nodes/weights are not physical absence.',
                tail_scope='Markov first/second-moment envelopes below bound ONLY the finite-Hermite distribution. Applying these numbers to the physical Coulomb state needs a continuum moment error bound not established here. Unconditional physical tail bound remains [0,1]. No T^3 or T^4 physical bound, Gaussian extrapolation, or core-removed claim is made.')


def distribution_summary(probability, kinetic):
    mean = float(probability @ kinetic)
    second = float(probability @ kinetic**2)
    return dict(normalization=float(probability.sum()), mean_kinetic_eV=mean,
                second_kinetic_moment_eV2=second,
                tails={str(threshold): dict(
                    finite_quadrature_mass=float(probability[kinetic > threshold].sum()),
                    finite_state_moment_envelope=min(1.0, mean / threshold, second / threshold**2),
                    physical_tail_lower=0.0, physical_tail_upper=1.0,
                    physical_point_Coulomb_contribution='unresolved; finite-Hermite ultraviolet cannot resolve the cusp')
                    for threshold in (100, 300, 1000)})
