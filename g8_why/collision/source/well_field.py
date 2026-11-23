"""Frozen44-electron PdD-vacancy field with point nuclei and resolved Pd cores.

The auxiliary cell has charge-1. Its G=0/background convention only fixes an
additive constant, removed at the occupied target's reference position. The
regularized PAW ionic potential is NOT used as a nuclear short-range wall.
"""
from itertools import product
from pathlib import Path
import hashlib

import numpy as np
from ase.units import Bohr, Hartree
from scipy.interpolate import CubicSpline
from scipy.special import erfc

COULOMB_EV_A = Hartree * Bohr


def file_binding(path):
    path = Path(path)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


class FrozenHostField:
    """Periodic valence-density field plus physical point-nucleus/core field."""
    def __init__(self, path):
        self.binding = file_binding(path)
        with np.load(path) as saved:
            density = saved['density_valence_e_A3']
            self.cell = saved['cell_A']
            self.positions = saved['positions_A']
            self.Z = saved['nuclear_Z']
            r = saved['pd_r_g_bohr'] * Bohr
            dr = saved['pd_dr_g_bohr'] * Bohr
            nc = saved['pd_nc_g'] / Bohr**3
        diagonal = np.diag(self.cell)
        if not np.allclose(self.cell, np.diag(diagonal), atol=1e-12, rtol=0):
            raise ValueError('This frozen PdD producer is an orthorhombic cell')
        self.lengths = diagonal
        self.volume = float(np.prod(diagonal))
        self.electrons = float(density.sum()*self.volume/density.size)
        if abs(self.electrons-44) > 1e-5:
            raise ValueError(f'Expected the same44-electron valence field, got{self.electrons}')
        # nc is the radial l=0 coefficient, not the physical spherical density.
        shell_charge = np.sqrt(4*np.pi)*nc*r*r*dr
        self.core_electrons = float(shell_charge.sum())
        if abs(self.core_electrons-36) > 1e-4:
            raise ValueError(f'Pd core normalization is not36:{self.core_electrons}')
        outer_charge = np.cumsum(shell_charge[::-1])[::-1]
        shell_over_r = np.divide(shell_charge, r, out=np.zeros_like(r), where=r > 0)
        outer_inverse_radius = np.cumsum(shell_over_r[::-1])[::-1]
        # A(r)/r = Qcore/r + Phi_core(r)/C. A(0)=Qcore; no forced normalization.
        A = outer_charge-r*outer_inverse_radius
        self.core_A = CubicSpline(r, A,
            bc_type=((1, -outer_inverse_radius[0]), (1, 0.)), extrapolate=False)
        self.core_radius_max = float(r[-1])
        self.ionic_charge = self.Z.astype(float)-np.where(self.Z == 46, self.core_electrons, 0.)
        self.g_axes = [2*np.pi*np.fft.fftfreq(n, d=L/n)
                       for n,L in zip(density.shape, diagonal)]
        gx,gy,gz = np.meshgrid(*self.g_axes, indexing='ij', sparse=True)
        g2 = gx*gx+gy*gy+gz*gz
        self.eta = 3/min(diagonal)
        charge_q = np.zeros(density.shape, complex)
        for charge, position in zip(self.ionic_charge, self.positions):
            charge_q += charge*np.exp(-1j*(gx*position[0]+gy*position[1]+gz*position[2]))
        numerator = (charge_q*np.exp(-g2/(4*self.eta**2))/self.volume
                     - np.fft.fftn(density)/density.size)
        self.phi_q = np.divide(4*np.pi*COULOMB_EV_A*numerator, g2,
                               out=np.zeros_like(numerator), where=g2 > 0)
        self.constant = -np.pi*COULOMB_EV_A*self.ionic_charge.sum()/(self.eta**2*self.volume)

    def core_correction(self, radius):
        radius = np.asarray(radius)
        if np.any(radius <= 0):
            raise ValueError('Point-nucleus coincidences have infinite potential')
        A = self.core_A(np.minimum(radius, self.core_radius_max))
        return COULOMB_EV_A*np.where(radius < self.core_radius_max, A/radius, 0.)

    def cartesian_grid(self, x, y, z):
        """Evaluate on arbitrary Cartesian quadrature nodes without FFT box edges."""
        axes = [np.atleast_1d(v) for v in (x,y,z)]
        phase = [np.exp(1j*np.outer(v,g)) for v,g in zip(axes,self.g_axes)]
        potential = np.einsum('ia,jb,kc,abc->ijk', *phase, self.phi_q,
                              optimize=True).real+self.constant
        for atom,charge,position in zip(self.Z,self.ionic_charge,self.positions):
            delta = [(v-p+L/2)%L-L/2 for v,p,L in zip(axes,position,self.lengths)]
            for image in product((-1,0,1), repeat=3):
                a,b,c = [d+n*L for d,n,L in zip(delta,image,self.lengths)]
                radius = np.sqrt(a[:,None,None]**2+b[None,:,None]**2+c[None,None,:]**2)
                if np.any(radius == 0):
                    raise ValueError('Quadrature node coincides with a point nucleus')
                potential += charge*COULOMB_EV_A*erfc(self.eta*radius)/radius
            if atom == 46:
                a,b,c = delta
                radius = np.sqrt(a[:,None,None]**2+b[None,:,None]**2+c[None,None,:]**2)
                potential += self.core_correction(radius)
        return potential

    def points(self, points):
        return np.asarray([self.cartesian_grid([x],[y],[z]).item() for x,y,z in points])

    def diagnostics(self):
        return dict(density_binding=self.binding, valence_electrons=self.electrons,
                    Pd_core_electrons=self.core_electrons,
                    ionic_charge_sum=float(self.ionic_charge.sum()),
                    net_auxiliary_cell_charge=float(self.ionic_charge.sum()-self.electrons),
                    Coulomb_eV_A=COULOMB_EV_A, Ewald_eta_A_inverse=self.eta,
                    real_image_shell=1,
                    scope='Frozen auxiliary electrons, physical point nuclei and Pd core; not the relaxed BO surface')
