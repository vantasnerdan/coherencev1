"""Zero-G periodic point/Gaussian nuclear Green functions (P255).

This qualifies the *bare* block. Replacing its regularizer within a common
periodic gauge is not by itself an isolated material-pair image correction.
The total neutral target-conditioned interaction needs its own cell control.
"""
import numpy as np
from scipy.special import erfc
from .collision_observables import COULOMB_EV_A


def _triplets(shell):
    axis=np.arange(-shell,shell+1)
    return np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1).reshape(-1,3)


class PeriodicNuclearGreen:
    def __init__(self,cell_A,real_shell=4,reciprocal_shell=6,eta_Ainv=None):
        self.cell=np.asarray(cell_A,dtype=float)
        self.volume=float(np.linalg.det(self.cell))
        if self.cell.shape!=(3,3) or not np.isfinite(self.cell).all() or self.volume<=0:
            raise ValueError('A finite positive-volume cell is required')
        self.eta=2/self.volume**(1/3) if eta_Ainv is None else float(eta_Ainv)
        if self.eta<=0 or not np.isfinite(self.eta):
            raise ValueError('A positive Ewald splitting scale is required')
        self.shifts=_triplets(real_shell)@self.cell
        reciprocal=2*np.pi*np.linalg.inv(self.cell).T
        modes=_triplets(reciprocal_shell)
        modes=modes[np.any(modes!=0,axis=1)]
        self.G=modes@reciprocal
        g2=np.einsum('ij,ij->i',self.G,self.G)
        self.weights=4*np.pi/self.volume*np.exp(-g2/(4*self.eta**2))/g2

    def _distances(self,displacement_A):
        point=np.asarray(displacement_A,dtype=float)
        if point.shape!=(3,) or not np.isfinite(point).all():
            raise ValueError('A finite target-relative displacement is required')
        distances=np.linalg.norm(self.shifts+point,axis=1)
        if np.min(distances)<1e-12:
            raise ValueError('Point-nuclear coalescence is singular')
        return point,distances

    def point_eV(self,displacement_A):
        point,distances=self._distances(displacement_A)
        real=np.sum(erfc(self.eta*distances)/distances)
        reciprocal=np.dot(self.weights,np.cos(self.G@point))
        return float(COULOMB_EV_A*(real+reciprocal-np.pi/(self.eta**2*self.volume)))

    def point_minus_gaussian_eV(self,displacement_A,pair_width_Ainv):
        """Exact Ewald split: point minus a normalized pair Gaussian.

        For nuclear densities exp(-alpha_i r^2), the pair Gaussian width is
        sqrt(alpha_1*alpha_2/(alpha_1+alpha_2)), in inverse Angstrom.
        """
        _,distances=self._distances(displacement_A)
        width=float(pair_width_Ainv)
        if width<=0 or not np.isfinite(width):
            raise ValueError('A finite positive Gaussian pair width is required')
        return float(COULOMB_EV_A*(np.sum(erfc(width*distances)/distances)
                                  -np.pi/(width**2*self.volume)))

    def gaussian_eV(self,displacement_A,pair_width_Ainv):
        return self.point_eV(displacement_A)-self.point_minus_gaussian_eV(displacement_A,pair_width_Ainv)
