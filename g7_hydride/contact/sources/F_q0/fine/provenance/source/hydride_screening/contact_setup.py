"""Built, projector-free Z=2 operator in the actual Hae absolute convention.

This object deliberately has no ``build`` method: create_setup's public
setupdata branch returns a prebuilt BaseSetup directly, without LeanSetup.
The analytic 1s basis is an initial guess, not a change to the Hamiltonian.
"""
import math

import numpy as np
from gpaw.ae import HydrogenAllElectronSetup
from gpaw.atom.radialgd import EquidistantRadialGridDescriptor
from gpaw.basis_data import Basis, BasisFunction
from gpaw.setup import BaseSetup
from gpaw.spline import Spline


class ContactHeliumSetup(BaseSetup):
    name = 'HaeZ2'

    def __init__(self):
        hydrogen = HydrogenAllElectronSetup(alpha1=10.0, alpha2=300.0)
        # Copy the actual operator's zero-sized correction data and normalized
        # compensation charge, not an independently reconstructed PAW He setup.
        self.__dict__.update(hydrogen.__dict__)
        self.symbol = 'He'
        self.Z = self.Nv = 2
        self.Nc = 0
        self.E = 0.0
        self.Delta0 = 2 * hydrogen.Delta0
        r = np.linspace(0.0, 2.0, 100)
        self.vbar = Spline.from_data(0, 2.0, 2 * hydrogen.vbar.map(r))
        self.M = 4 * hydrogen.M
        # A normalized real hydrogenic Z=2 radial 1s, smoothly truncated at
        # 8 Bohr; normalization uses the actual initial-guess radial mesh.
        rgd = EquidistantRadialGridDescriptor(0.01, 801)
        r = rgd.r_g
        radial = np.exp(-2 * r) * (1 - (r / r[-1]) ** 2) ** 2
        radial /= math.sqrt(float(np.sum(radial ** 2 * r ** 2) * (r[1] - r[0])))
        self.basis = Basis('He', 'contact-real-1s', rgd,
                           bf_j=[BasisFunction(n=1, l=0, rc=float(r[-1]),
                                               phit_g=radial, type='real-hydrogenic-1s')])
        self.basis_functions_J = self.basis.tosplines()
        self.pseudo_partial_waves_j = None
        self.f_j = [2.0]
        self.nao = self.basis.nao

    def get_default_nbands(self):
        return 1

    def print_info(self, text):
        text('Hae-derived Z=2 all-electron contact operator; Eref=0, M=4*M_H')

    def descriptor(self, sources):
        return dict(name=self.name, type=self.type, Z=self.Z, Nv=self.Nv,
                    Nc=self.Nc, ni=self.ni, E=self.E, alpha1=self.alpha1,
                    alpha2=self.alpha2, Delta0=float(self.Delta0), M=float(self.M),
                    occupation=list(self.f_j), source=sources,
                    initial_guess='normalized real hydrogenic Z=2 1s; smooth 8-Bohr cutoff',
                    ghat='unchanged normalized alpha1 Gaussian',
                    vbar_scale=2, M_scale=4,
                    public_path='prebuilt BaseSetup -> setups[He] -> create_setup(setupdata); no build method')
