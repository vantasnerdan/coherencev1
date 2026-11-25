"""NiH [111] signed-anchor model; Cartesian Å -> eV, eV/Å."""
from hashlib import sha256
import json

import numpy as np
from scipy.interpolate import CubicHermiteSpline

SIGNED_SHA256 = '18fa6232525243391b84f2fb495cfe9947f7d41eef270931fb3dc8b99857e28d'
DIRECTIONS = np.array([[1., 0., 0.], [1., 1., 0.], [1., 1., 1.]])
DIRECTIONS /= np.linalg.norm(DIRECTIONS, axis=1)[:, None]
VARIANTS = ('old', 'w', 'w2', 'w_stiff_exterior')
RMAX = 0.80


def load_anchors(path, moving_index):
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != SIGNED_SHA256:
        raise ValueError('Signed wide-well export differs from the accepted input')
    source = json.loads(raw)
    rows = [row for row in source['rows'] if row['material'] == 'NiH']
    rows.sort(key=lambda row: row['displacement_A'])
    if [row['displacement_A'] for row in rows] != [-.8, -.48, -.24, .24, .48, .8]:
        raise ValueError('Expected precisely the six signed NiH anchors')
    anchors = []
    for row in rows:
        s = row['displacement_A']
        u = np.asarray(row['ray_unit'])
        f = np.asarray(row['forces_eV_A'][moving_index]) @ u
        kind = 'native_scf' if s > 0 else 'symmetry_derived'
        if (row['moving_index'] != moving_index or not np.allclose(u, DIRECTIONS[2], atol=1e-14, rtol=0)
                or row['sample_kind'] != kind or abs(f - row['native_force_parallel_eV_A']) > 1e-12):
            raise ValueError(f"Inconsistent signed NiH row: {row['case']}")
        anchors.append(dict(signed_displacement_A=s, delta_energy_free_eV=row['delta_energy_free_eV'],
                            force_parallel_eV_A=row['native_force_parallel_eV_A'], case=row['case'],
                            sample_kind=kind, parent_native_case=row.get('parent_native_case'),
                            native_row_sha256=row['native_row_sha256']))
    for neg, pos in zip(reversed(anchors[:3]), anchors[3:]):
        if (neg['parent_native_case'] != pos['case']
                or abs(neg['delta_energy_free_eV'] - pos['delta_energy_free_eV']) > 1e-12
                or abs(neg['force_parallel_eV_A'] + pos['force_parallel_eV_A']) > 1e-12):
            raise ValueError('Negative NiH rows are not the declared inversion derivatives')
    return anchors


class NiHWell:
    def __init__(self, k, a, b, anchors):
        self.k, self.a, self.b = k, a, b
        self.q111 = (a + 2 * b) / 3
        self.q110 = (a + b) / 2
        self.anchors = anchors
        positive = anchors[3:]
        r = np.array([row['signed_displacement_A'] for row in positive])
        e = np.array([row['delta_energy_free_eV'] for row in positive])
        f = np.array([row['force_parallel_eV_A'] for row in positive])
        self.spline = CubicHermiteSpline(np.r_[0., r * r], np.r_[0., e],
                                        np.r_[.5 * k, -f / (2 * r)], extrapolate=False)
        self.outer_e, self.outer_slope = e[-1], -f[-1]
        # On each interval p'(u) is a quadratic: exact location of its minimum.
        minima = []
        for j, width in enumerate(np.diff(self.spline.x)):
            c3, c2, c1, _ = self.spline.c[:, j]
            z = [0., width]
            if c3 > 0 and 0 < -c2 / (3 * c3) < width:
                z.append(-c2 / (3 * c3))
            minima.append(min((3 * c3 * x + 2 * c2) * x + c1 for x in z))
        self.min_dE_du = min(minima)
        if min(k, a, self.q111, self.min_dE_du, self.outer_slope) <= 0:
            raise ValueError('NiH closure lost positivity or confinement')

    def energy_force(self, x, variant, gradient=True):
        xx = x * x
        u = xx.sum(axis=-1)
        old = .5 * self.k * u + .25 * ((self.a - self.b) * (xx * xx).sum(axis=-1) + self.b * u * u)
        if gradient:
            g = (self.k + (self.a - self.b) * xx + self.b * u[..., None]) * x
        if variant == 'old':
            return old, -g if gradient else None
        r = np.sqrt(u)
        dr = np.maximum(r - RMAX, 0.)
        # The stress changes only d²E/dr² outside the last measured anchor.
        curvature = self.k * (8 if variant == 'w_stiff_exterior' else 1)
        clipped = np.minimum(u, RMAX * RMAX)
        radial = self.spline(clipped) + dr * (self.outer_slope + .5 * curvature * dr)
        inc = radial - (.5 * self.k + .25 * self.q111 * u) * u
        inv_u = np.divide(1., u, out=np.zeros_like(u), where=u > 0)
        l = xx * inv_u[..., None]
        w = 27 * np.prod(l, axis=-1)
        w_eff = w * w if variant == 'w2' else w
        energy = old + w_eff * inc
        if not gradient:
            return energy, None
        radial_slope = np.where(r <= RMAX, 2 * r * self.spline(clipped, 1), self.outer_slope + curvature * dr)
        inv_r = np.divide(1., r, out=np.zeros_like(r), where=r > 0)
        dinc_over_r = radial_slope * inv_r - self.k - self.q111 * u
        # ∂w/∂x_i = (54 x_i l_j l_k - 6 w x_i)/u. No division by x_i.
        dw = ((54 * np.roll(l, 1, axis=-1) * np.roll(l, 2, axis=-1) - 6 * w[..., None])
              * x * inv_u[..., None])
        if variant == 'w2':
            dw *= 2 * w[..., None]
        return energy, -(g + (w_eff * dinc_over_r)[..., None] * x + inc[..., None] * dw)

    def probe(self):
        """Differentiate the real signed model independently of its force code."""
        s = np.array([row['signed_displacement_A'] for row in self.anchors])
        energy = np.array([row['delta_energy_free_eV'] for row in self.anchors])
        force = np.array([row['force_parallel_eV_A'] for row in self.anchors])
        n = DIRECTIONS[2]
        x = s[:, None] * n
        hc = 1e-6
        found = {}
        for variant in VARIANTS[1:]:
            e, f = self.energy_force(x, variant)
            eh = [self.energy_force(x + side * hc * n, variant, False)[0] for side in (-1, 1)]
            left, right = (e - eh[0]) / hc, (eh[1] - e) / hc
            fd = -.5 * (left + right)
            keep_x = s[:, None, None] * DIRECTIONS[None, :2, :]
            old_e, old_f = self.energy_force(keep_x, 'old')
            new_e, new_f = self.energy_force(keep_x, variant)
            point = np.array([.37, -.21, .18])
            _, point_f = self.energy_force(point, variant)
            finite_f = np.array([-(self.energy_force(point + hc * d, variant, False)[0]
                                   - self.energy_force(point - hc * d, variant, False)[0]) / (2 * hc)
                                 for d in np.eye(3)])
            # An outward-force convention accepted as fixed-axis would fail at the negative branch.
            wrong_sign_error = np.max(abs(fd + force))
            info = dict(max_energy_anchor_error_eV=np.max(abs(e - energy)),
                        max_force_anchor_error_eV_A=np.max(abs(f @ n - force)),
                        max_anchor_force_finite_difference_error_eV_A=np.max(abs(fd - force)),
                        max_one_sided_dEds_jump_eV_A=np.max(abs(left - right)),
                        wrong_signed_force_error_eV_A=wrong_sign_error,
                        old_100_110_energy_max_abs_change_eV=np.max(abs(old_e - new_e)),
                        old_100_110_force_max_abs_change_eV_A=np.max(abs(old_f - new_f)),
                        off_line_force_fd_max_abs_error_eV_A=np.max(abs(point_f - finite_f)))
            if (max(info['max_energy_anchor_error_eV'], info['max_force_anchor_error_eV_A']) > 1e-12
                    or max(info['max_anchor_force_finite_difference_error_eV_A'],
                           info['max_one_sided_dEds_jump_eV_A']) > 2e-5
                    or wrong_sign_error < 1.
                    or max(info['old_100_110_energy_max_abs_change_eV'],
                           info['old_100_110_force_max_abs_change_eV_A']) > 1e-12
                    or info['off_line_force_fd_max_abs_error_eV_A'] > 1e-8):
                raise ValueError(f'{variant}: signed fit consumer failed: {info}')
            found[variant] = info
        return found
