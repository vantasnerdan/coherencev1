"""Saved-state convention sensitivity and floating-point omission diagnostics."""
from pathlib import Path

import numpy as np

from quantum_dvr import KINETIC_EV_A2, KB_EV_K


def saved_convention_sensitivity(output):
    """A missing Fourier phase maps position to momentum and must change K."""
    output = Path(output)
    with np.load(output / 'momentum-position.npz') as saved:
        q = saved['gh_nodes']
        probability = saved['ground_position_cubic_probability']
        false_kinetic = float((probability * (q[:, None, None]**2 + q[None, :, None]**2
                                             + q[None, None, :]**2)).sum() * KINETIC_EV_A2 / 0.30**2)
    with np.load(output / 'joint-spectrum.npz') as saved:
        true_kinetic = float(saved['ground_probabilities'] @ saved['kinetic_eV'])
    return dict(mutation='Omit every (-i)^Hermite-index Fourier phase on the actual saved ground state',
                true_ground_mean_kinetic_eV=true_kinetic, false_ground_mean_kinetic_eV=false_kinetic,
                difference_eV=abs(true_kinetic - false_kinetic),
                exposed=abs(true_kinetic - false_kinetic) > 1e-8)


def canonical_underflow_scope(output):
    with np.load(Path(output) / 'nuclear-spectrum.npz') as saved:
        energy = saved['eigenenergies_eV']
        shifted_log_partition = float(saved['log_partition']) + energy[0] / (KB_EV_K * float(saved['temperature_K']))
        underflow = saved['canonical_weights'] == 0
        # Each zero float has exp(exponent) below the least positive normal.
        # This conservative logarithmic bound remains representable even when
        # the actual missing sum is smaller than the least subnormal float.
        log_bound = (float(np.log(underflow.sum()) + np.log(np.finfo(float).tiny)
                           - shifted_log_partition) if underflow.any() else None)
        return dict(zero_weight_state_count=int(underflow.sum()),
                    log_omitted_canonical_probability_upper=log_bound,
                    scope='Numerical underflow only; all eigenvalues/states were solved and included. Joint omission_probability=0 means no deliberately omitted quadrature states/shells, not zero physical continuum tail.')
