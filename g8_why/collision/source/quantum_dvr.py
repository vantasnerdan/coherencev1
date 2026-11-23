"""Complete finite-Hermite DVR nuclear Hamiltonian, without a wavefunction box.

Coordinates are target-relative Angstrom; energies are eV. The potential uses
n-point Hermite collocation (DVR), not an exact potential Galerkin integral.
The kinetic operator is the exact projected continuum p^2, not (projected p)^2.
"""
from dataclasses import dataclass

import numpy as np
from scipy.constants import physical_constants, hbar, electron_volt, Boltzmann
from scipy.linalg import eigh
from scipy.special import roots_hermite

DEUTERON_MASS_KG = physical_constants['deuteron mass'][0]
KINETIC_EV_A2 = hbar**2 / (2 * DEUTERON_MASS_KG * electron_volt * 1e-20)
KB_EV_K = Boltzmann / electron_volt


def hermite_polynomials(nodes, n):
    """Normalized physicists' Hermite polynomials, weight exp(-q^2)."""
    values = np.empty((len(nodes), n))
    values[:, 0] = np.pi**(-0.25)
    if n > 1:
        values[:, 1] = np.sqrt(2) * nodes * values[:, 0]
    for j in range(1, n - 1):
        values[:, j + 1] = (np.sqrt(2 / (j + 1)) * nodes * values[:, j]
                             - np.sqrt(j / (j + 1)) * values[:, j - 1])
    return values


def quadrature_map(n, order):
    q, weights = roots_hermite(order)
    return q, weights, np.sqrt(weights)[:, None] * hermite_polynomials(q, n)


def kinetic_matrix(n, length_A):
    diagonal = np.arange(n, dtype=float) + 0.5
    matrix = np.diag(diagonal)
    off = -0.5 * np.sqrt(np.arange(1, n - 1) * np.arange(2, n))
    matrix += np.diag(off, 2) + np.diag(off, -2)
    return matrix * KINETIC_EV_A2 / length_A**2


def tensor_transform(vectors, matrix):
    """Apply the same one-axis map on all axes; last axis labels states."""
    return np.einsum('ai,bj,ck,ijks->abcs', matrix, matrix, matrix,
                     vectors, optimize=True)


def tensor_kinetic(vectors, kinetic_1d):
    return (np.einsum('ai,ijks->ajks', kinetic_1d, vectors, optimize=True)
            + np.einsum('bj,ijks->ibks', kinetic_1d, vectors, optimize=True)
            + np.einsum('ck,ijks->ijcs', kinetic_1d, vectors, optimize=True))


@dataclass
class NuclearSpectrum:
    energies: np.ndarray
    vectors_dvr: np.ndarray
    potential: np.ndarray
    kinetic_dvr: np.ndarray
    nodes_A: np.ndarray
    dvr_map: np.ndarray
    length_A: float


def solve(well, n, length_A=0.30):
    q, _, transform = quadrature_map(n, n)
    nodes = length_A * q
    potential = well.cartesian_grid(nodes, nodes, nodes)
    kinetic = transform @ kinetic_matrix(n, length_A) @ transform.T
    identity = np.eye(n)
    hamiltonian = (np.kron(np.kron(kinetic, identity), identity)
                   + np.kron(np.kron(identity, kinetic), identity)
                   + np.kron(np.kron(identity, identity), kinetic))
    hamiltonian.flat[::n**3 + 1] += potential.ravel()
    energies, vectors = eigh(hamiltonian, overwrite_a=True, check_finite=True,
                              driver='evd')
    return NuclearSpectrum(energies, vectors, potential, kinetic, nodes,
                           transform, length_A)


def canonical_weights(energies, temperature_K):
    excitation = energies - energies[0]
    if temperature_K == 0:
        weights = np.zeros_like(energies)
        weights[0] = 1
        return weights, 0.0
    exponent = -excitation / (KB_EV_K * temperature_K)
    shifted_partition = np.exp(exponent).sum()
    return (np.exp(exponent) / shifted_partition,
            float(np.log(shifted_partition) - energies[0] / (KB_EV_K * temperature_K)))


def save_spectrum(path, spectrum, weights, log_partition, temperature_K):
    np.savez(path, eigenenergies_eV=spectrum.energies,
             eigenvectors_dvr=spectrum.vectors_dvr,
             potential_grid_eV=spectrum.potential,
             kinetic_1d_dvr_eV=spectrum.kinetic_dvr,
             dvr_from_fbr=spectrum.dvr_map, coordinate_nodes_A=spectrum.nodes_A,
             basis_length_A=spectrum.length_A, deuteron_mass_kg=DEUTERON_MASS_KG,
             kinetic_eV_A2=KINETIC_EV_A2, canonical_weights=weights,
             log_partition=log_partition, temperature_K=temperature_K)


def consume_spectrum(path):
    """Replay the saved separable operator against every saved eigenvector."""
    with np.load(path) as saved:
        energies = saved['eigenenergies_eV']
        vectors = saved['eigenvectors_dvr']
        potential = saved['potential_grid_eV']
        kinetic = saved['kinetic_1d_dvr_eV']
        n = len(saved['coordinate_nodes_A'])
        residual = 0.0
        norm_error = 0.0
        for start in range(0, len(energies), 64):
            block = vectors[:, start:start + 64].reshape(n, n, n, -1)
            action = tensor_kinetic(block, kinetic) + potential[..., None] * block
            delta = action.reshape(n**3, -1) - block.reshape(n**3, -1) * energies[start:start + 64]
            residual = max(residual, float(np.linalg.norm(delta, axis=0).max()))
            norm_error = max(norm_error, float(abs((block**2).sum(axis=(0, 1, 2)) - 1).max()))
        weights, log_partition = canonical_weights(energies, float(saved['temperature_K']))
        return dict(max_eigenresidual_eV=residual, max_vector_norm_error=norm_error,
                    weight_sum=float(saved['canonical_weights'].sum()),
                    weight_replay_error=float(abs(weights - saved['canonical_weights']).max()),
                    log_partition_replay_error=abs(log_partition - float(saved['log_partition'])))
