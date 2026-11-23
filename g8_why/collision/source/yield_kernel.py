"""Fixed transported-profile isotropic free-target diagnostic, not a bound-state rate.

All nuclear cross sections, transport nodes, species weights, and paper readout
are consumed from the pinned g8_why/g4_beam implementation. A supplied bound-state
momentum density can be integrated here only as an explicitly on-shell diagnostic;
binding/excitation/recoil energy conservation needs a different collision operator.
"""
from functools import lru_cache
import hashlib
import importlib.util
from pathlib import Path
import sys

import numpy as np
from scipy.constants import atomic_mass, electron_volt, proton_mass

DEFAULT_BASE = Path('<PRIVATE_PATH>')
BASE_COMMIT = '231859350047c2ed9443394e512ea1697c987664'


@lru_cache(None)
def source_common(base_root=DEFAULT_BASE):
    base = Path(base_root).resolve()
    sys.path.insert(0, str(base))
    spec = importlib.util.spec_from_file_location('g8_fixed_profile_common', base/'g8_why/common.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FixedProfileShellKernel:
    """K(E_nominal,T_target,ch), in physical reactions per true charged C.

    T_target is LAB kinetic energy in eV. E_nominal is stationary-target CM keV;
    incident D+ LAB=2E_nominal and D2+ per-D LAB=E_nominal before transport.
    At each actual transport LAB energy L (keV), equal-mass relative CM energy
    is (L+T-2sqrt(LT)mu)/2, with T now in keV and mu=v_b dot v_t/(v_b v_t).
    The multiplier is sigma(E_relative,+40eV)*v_relative/v_beam, not sigma alone.
    Unlike transport's stationary 50eV CM cutoff, moving relative energies can
    lie below 50eV: the unchanged Eq7 law is evaluated there without a new floor.
    """

    def __init__(self, base_root=DEFAULT_BASE, angular_order=64):
        self.base_root = Path(base_root).resolve()
        self.common = source_common(self.base_root)
        self.mu, self.angular_weights = np.polynomial.legendre.leggauss(angular_order)
        self.angular_weights /= 2
        self.angular_order = angular_order
        self.deuteron_mass_kg = self.common.DD.M_D*atomic_mass

    def local_rate(self, lab_keV, target_eV, channel='p'):
        lab = np.atleast_1d(np.asarray(lab_keV, float))[:, None]
        target = float(target_eV)/1000
        if target < 0 or not np.isfinite(target) or np.any(lab <= 0) or np.any(~np.isfinite(lab)):
            raise ValueError('Positive finite beam and nonnegative finite target energies required')
        if channel not in ('p', 'n'):
            raise ValueError('Channel must be p or n')
        relative = (lab+target-2*np.sqrt(lab*target)*self.mu)/2
        # Eq7 sigma diverges at zero relative E but sigma*v has a finite limit.
        # Using the analytic product there preserves the original S and +40eV law.
        positive = np.maximum(relative, np.finfo(float).tiny)
        flux = np.sqrt(2*positive/lab)
        rate = self.common.DD.sigma(positive, self.common.U, channel)*flux
        return rate, relative

    def shell(self, material, nominal_cm_keV, target_eV, channel='p'):
        lab, path_weights = self.common.response(material, float(nominal_cm_keV))
        rate, _ = self.local_rate(lab, target_eV, channel)
        return float(path_weights @ (rate @ self.angular_weights))

    def matrix(self, material, nominal_cm_keV, target_eV, channel='p'):
        return np.array([[self.shell(material, e, t, channel) for t in target_eV]
                         for e in nominal_cm_keV])

    def mixture(self, material, nominal_cm_keV, target_eV, probabilities, channel='p'):
        probabilities = np.asarray(probabilities, float)
        target_eV = np.asarray(target_eV, float)
        if probabilities.shape != target_eV.shape or probabilities.ndim != 1:
            raise ValueError('One probability is required per target shell')
        if np.any(~np.isfinite(probabilities)) or np.any(probabilities < 0) or abs(probabilities.sum()-1) > 1e-9:
            raise ValueError('Nonnegative normalized target probabilities required')
        active = probabilities > 0
        return self.matrix(material, nominal_cm_keV, target_eV[active], channel) @ probabilities[active]

    def reported(self, material, energies, target_eV, probabilities, current_power=0):
        cache = {}
        def physical(e, channel):
            key = float(e), channel
            if key not in cache:
                cache[key] = float(self.mixture(material, [e], target_eV, probabilities, channel)[0])
            return cache[key]
        return self.common.readout(material, np.asarray(energies), physical,
                                   current_power=current_power)

    def shell_moments(self, material, nominal_cm_keV, target_eV, proton_energy_keV=3020., detector_deg=90.):
        """Reaction-weighted COM velocity and nominal proton Doppler moments.

        Beam defines z; target azimuth is integrated exactly. Proton direction is
        approximated by its CM direction (not a finite-aperture LAB detector).
        Doppler values use a declared FIXED 3020keV CM proton energy by default,
        not a nuclear Q-value calculation or an instrumental line convolution.
        """
        lab, path_weights = self.common.response(material, float(nominal_cm_keV))
        rate, _ = self.local_rate(lab, target_eV, 'p')
        weights = path_weights[:, None]*rate*self.angular_weights
        norm = weights.sum()
        mass = self.deuteron_mass_kg
        vb = np.sqrt(2*lab[:, None]*1000*electron_volt/mass)
        vt = np.sqrt(2*float(target_eV)*electron_volt/mass)
        vz = (vb+vt*self.mu)/2
        vperp2 = np.broadcast_to(vt**2*(1-self.mu**2)/4, vz.shape)
        v2 = vz**2+vperp2
        average = lambda value: float(np.sum(weights*value)/norm)
        mean_z = average(vz)
        var_z = average((vz-mean_z)**2)
        x2 = average(vperp2/2)
        theta = np.deg2rad(detector_deg)
        momentum = np.sqrt(2*proton_mass*proton_energy_keV*1000*electron_volt)
        mean_phi_shift_J = proton_mass*v2/2 + momentum*np.cos(theta)*vz
        mean_shift_J = average(mean_phi_shift_J)
        variance_J2 = average((mean_phi_shift_J-mean_shift_J)**2 + momentum**2*np.sin(theta)**2*vperp2/2)
        return dict(physical_p_yield_per_true_C=float(norm), mean_COM_vz_m_s=mean_z,
                    variance_COM_vz_m2_s2=var_z, mean_COM_vx2_m2_s2=x2,
                    mean_COM_speed2_m2_s2=average(v2), detector_deg=float(detector_deg),
                    assumed_fixed_CM_proton_energy_keV=float(proton_energy_keV),
                    mean_proton_shift_eV=mean_shift_J/electron_volt,
                    proton_shift_std_eV=np.sqrt(max(0., variance_J2))/electron_volt)

    def source_bindings(self):
        paths = ['g8_why/common.py', 'g8_why/reaction.py', 'g8_why/reaction_results.json',
                 'g8_why/reaction_retained.csv', 'g8_why/reaction_shape_bound.csv',
                 'g4_beam/merit/api.py', 'g4_beam/merit/profiles.json',
                 'g4_beam/forward/model.py', 'g4_beam/forward/consumer.py',
                 'g4_beam/forward/profile_thermal.py', 'g4_beam/forward/errbars.csv']
        return dict(base_commit=BASE_COMMIT, base_root=str(self.base_root), files=[
            dict(path=p, bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
            for p in paths for data in [(self.base_root/p).read_bytes()]])
