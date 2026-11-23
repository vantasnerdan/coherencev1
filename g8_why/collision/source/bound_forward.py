"""Energy-conserving spectator/plane-wave impulse approximation at fixed profile.

A bound target's momentum component K is NOT additional initial kinetic energy.
For equal masses, available relative energy is
  Erel = (L-K)/2 - sqrt(L*K)*mu + excitation - removal.
The entrance rate uses the unchanged on-shell DD matrix at that available energy.
This is an explicit off-shell prescription, not an exact many-body fusion kernel.
The residual host is an infinite-mass spectator; nonnegative removal only lowers
available energy. Zero removal is a favorable model, not a computed binding energy.
"""
import numpy as np
from scipy.constants import electron_volt, proton_mass
from yield_kernel import FixedProfileShellKernel, DEFAULT_BASE


class BoundSpectatorKernel:
    def __init__(self, base_root=DEFAULT_BASE, angular_order=32):
        reference = FixedProfileShellKernel(base_root, angular_order)
        self.common = reference.common
        self.mu = reference.mu
        self.angular_weights = reference.angular_weights
        self.deuteron_mass_kg = reference.deuteron_mass_kg

    def local_samples(self, lab_keV, target_eV, channel='p', *, excitation_eV=0., removal_eV=0.):
        """Return rate including angular Jacobian, canonical mu, and available Erel.

        u=sqrt(Erel) smooths the screened Eq7 threshold integrable in angle;
        no negative-energy component is promoted to a positive on-shell state.
        Subthreshold virtual nuclear amplitudes are outside this approximation.
        """
        lab = np.atleast_1d(np.asarray(lab_keV, float))[:, None]
        target = float(target_eV)/1000
        if target < 0 or np.any(lab <= 0) or not np.isfinite(target+excitation_eV+removal_eV):
            raise ValueError('Finite nonnegative target and positive beam energies required')
        if channel not in ('p', 'n'):
            raise ValueError('Unknown D-D branch')
        available = (lab-target)/2 + (excitation_eV-removal_eV)/1000
        shape = (len(lab), len(self.mu))
        if target == 0:
            relative = np.broadcast_to(available, shape)
            cosine = np.broadcast_to(self.mu, shape)
            jacobian = np.ones(shape)
        else:
            coupling = np.sqrt(lab*target)
            lower = np.sqrt(np.maximum(available-coupling, 0))
            upper = np.sqrt(np.maximum(available+coupling, 0))
            u = (upper+lower)/2 + (upper-lower)*self.mu/2
            relative = u*u
            cosine = (available-relative)/coupling
            jacobian = u*(upper-lower)/coupling
        full_lab = np.broadcast_to(lab, shape)
        valid = (relative > 0) & (jacobian > 0)
        result = np.zeros(shape)
        result[valid] = (self.common.DD.sigma(relative[valid], self.common.U, channel)
                         * np.sqrt(2*relative[valid]/full_lab[valid]) * jacobian[valid])
        return result, cosine, relative

    def shell(self, material, nominal_cm_keV, target_eV, channel='p', **state):
        lab, path = self.common.response(material, float(nominal_cm_keV))
        rate, _, _ = self.local_samples(lab, target_eV, channel, **state)
        return float(path @ (rate @ self.angular_weights))

    def fold(self, material, energies, states, channel='p', *, removal_eV=0.):
        """States carry level excitation, kinetic nodes and joint probabilities."""
        total = sum(float(np.sum(s['probabilities'])) for s in states)
        if not np.isclose(total, 1., atol=1e-8, rtol=0):
            raise ValueError(f'Joint spectral probabilities sum to {total}, not one')
        result = np.zeros(len(energies))
        for state in states:
            kinetic = np.asarray(state['kinetic_eV'])
            probabilities = np.asarray(state['probabilities'])
            if kinetic.shape != probabilities.shape or np.any(probabilities < 0):
                raise ValueError('Matched nonnegative spectral probabilities required')
            for kinetic_eV, probability in zip(kinetic, probabilities):
                if probability:
                    result += probability*np.array([self.shell(material, e, kinetic_eV, channel,
                        excitation_eV=state['excitation_eV'], removal_eV=removal_eV) for e in energies])
        return result

    def proton_moments(self, material, nominal_cm_keV, target_eV, *, excitation_eV=0.,
                       removal_eV=0., detector_deg=90., proton_energy_keV=3020.):
        """Reaction-weighted spectator COM Doppler moments; no detector convolution."""
        lab, path = self.common.response(material, float(nominal_cm_keV))
        rate, cosine, _ = self.local_samples(lab, target_eV, 'p',
            excitation_eV=excitation_eV, removal_eV=removal_eV)
        weights = path[:, None]*rate*self.angular_weights
        norm = float(weights.sum())
        if norm <= 0:
            raise ValueError('No positive reaction weight for these moments')
        valid = weights > 0
        mass = self.deuteron_mass_kg
        vb = np.sqrt(2*lab[:, None]*1000*electron_volt/mass)
        vt = np.sqrt(2*target_eV*electron_volt/mass)
        vz = (vb+vt*cosine)/2
        transverse2 = vt*vt*(1-cosine*cosine)/4
        speed2 = vz*vz+transverse2
        average = lambda value: float(np.sum(weights[valid]*np.broadcast_to(value, weights.shape)[valid])/norm)
        theta = np.deg2rad(detector_deg)
        momentum = np.sqrt(2*proton_mass*proton_energy_keV*1000*electron_volt)
        shift = proton_mass*speed2/2+momentum*np.cos(theta)*vz
        mean = average(shift)
        variance = average((shift-mean)**2+momentum**2*np.sin(theta)**2*transverse2/2)
        return dict(physical_p_yield_per_true_C=norm, detector_deg=detector_deg,
                    fixed_proton_CM_energy_keV=proton_energy_keV,
                    mean_proton_shift_eV=mean/electron_volt,
                    proton_shift_std_eV=np.sqrt(variance)/electron_volt,
                    mean_COM_vz_m_s=average(vz), mean_COM_vx2_m2_s2=average(transverse2/2))
