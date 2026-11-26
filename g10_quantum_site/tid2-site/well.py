"""Regular T_d interpolant of signed frozen-ion rays; see model.md."""
import json
from pathlib import Path

import numpy as np
from scipy.constants import atomic_mass, electron_volt, hbar
from scipy.interpolate import CubicHermiteSpline

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
MASSES = {'H': 1.00782503223, 'D': 2.01410177812, 'T': 3.01604928199}
KB = 8.617333262145e-5
HBAR_EV_FS = hbar / electron_volt * 1e15


def kinetic(mass_u):
    return hbar**2 / (2 * mass_u * atomic_mass) / electron_volt * 1e20


class Well:
    def __init__(self, variant='base'):
        self.variant = variant
        self.source = json.loads((HERE / 'inputs/cuts.json').read_text())
        rays, slopes = [], []
        for label in ('100', '110', '111'):
            energies, derivatives = [], []
            for row in self.source['directions'][label]:
                minus, plus = row['signed_points']
                energies.append([minus['delta_energy_free_eV'], plus['delta_energy_free_eV']])
                # Radial outward derivatives, not signed-coordinate derivatives.
                derivatives.append([minus['target_force_parallel_eV_A'], -plus['target_force_parallel_eV_A']])
            rays.append(np.asarray(energies)); slopes.append(np.asarray(derivatives))
        r = self.r = np.array([.12, .24, .48, .8])
        t = r*r
        e0, e1 = rays[0].mean(axis=1), rays[1].mean(axis=1)
        g0, g1 = slopes[0].mean(axis=1), slopes[1].mean(axis=1)
        even, odd = rays[2].mean(axis=1), (rays[2][:, 1]-rays[2][:, 0])/2
        ge, go = slopes[2].mean(axis=1), (slopes[2][:, 1]-slopes[2][:, 0])/2
        nums = [e0, e1-e0, even-e0-4*(e1-e0)/3, odd]
        ders = [g0, g1-g0, ge-g0-4*(g1-g0)/3, go]
        self.coefficients = []
        for num, der, power in zip(nums, ders, (2, 4, 6, 3)):
            value = num/r**power
            derivative = (der/r**power-power*num/r**(power+1))/(2*r)
            self.coefficients.append((CubicHermiteSpline(t, value, derivative), value[0], derivative[0]))
        self.closure_scale = .05/self.null_shape(np.array([1., 2., 3.])/np.sqrt(14)*.4)

    @staticmethod
    def null_shape(x):
        q = np.asarray(x)**2
        t = q.sum(axis=-1)
        return ((q[..., 0]-q[..., 1])*(q[..., 1]-q[..., 2])*(q[..., 2]-q[..., 0]))**2/(t+.16)**4

    def interior(self, x, derivative=False):
        x = np.asarray(x)
        q = x*x; t = q.sum(axis=-1)
        s2 = q[..., 0]*q[..., 1]+q[..., 1]*q[..., 2]+q[..., 2]*q[..., 0]
        s3 = q.prod(axis=-1)
        multipliers = [t, 4*s2, 27*s3, 3*np.sqrt(3)*x.prod(axis=-1)]
        powers = [2, 4, 6, 3]
        out = np.zeros_like(t)
        for (spl, val, der), mul, power in zip(self.coefficients, multipliers, powers):
            f = np.where(t < .0144, val+der*(t-.0144), spl(t))
            if derivative:
                df = np.where(t < .0144, der, spl(t, 1))
                out += np.divide(mul*(power*f+2*t*df), np.sqrt(t), out=np.zeros_like(t), where=t > 0)
            else:
                out += mul*f
        return out

    def __call__(self, x):
        x = np.asarray(x)
        radius = np.linalg.norm(x, axis=-1)
        scale = np.minimum(1., np.divide(.8, radius, out=np.ones_like(radius), where=radius > 0))
        edge = x*scale[..., None]
        extra = np.maximum(radius-.8, 0)
        curvature = 4. if self.variant == 'base' else 8.
        value = self.interior(edge)+extra*self.interior(edge, True)+curvature/2*extra**2
        if self.variant != 'base':
            value += self.closure_scale*self.null_shape(x)
        return value

    def anchors(self):
        rows = []
        for label, direction in [('100', [1, 0, 0]), ('110', [1, 1, 0]), ('111', [1, 1, 1])]:
            direction = np.asarray(direction)/np.linalg.norm(direction)
            for group in self.source['directions'][label]:
                for point in group['signed_points']:
                    u = point['displacement_A']; x = u*direction
                    numerical_slope = float((self(x+1e-6*direction)-self(x-1e-6*direction))/2e-6)
                    rows.append(dict(direction=label, signed_u_A=u, source_eV=point['delta_energy_free_eV'],
                                     model_eV=float(self(x)), source_slope_eV_A=-point['target_force_parallel_eV_A'],
                                     model_slope_eV_A=numerical_slope))
        return rows
