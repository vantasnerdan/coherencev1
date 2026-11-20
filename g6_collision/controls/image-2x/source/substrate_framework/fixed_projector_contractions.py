"""Exact dense contractions of Native GPAW's fixed serial FD PAW projectors.

Native ``add`` constructs every projector at every actual Bloch q. This only
changes its contraction implementation; it does not approximate a projector,
Hamiltonian, overlap, occupation or propagation scheme. Moving a nonempty
projector centre is refused before its Native position update.
"""
from time import perf_counter
import numpy as np


class _FixedProjectors:
    def __init__(self, calc):
        self.pt = calc.wfs.pt
        if calc.wfs.world.size != 1 or self.pt.gd.comm.size != 1:
            raise ValueError('Frozen dense projectors require the serial case')
        if self.pt.xp is not np or np.dtype(self.pt.dtype) != np.dtype(complex):
            raise ValueError('Frozen dense projectors require NumPy complex FD arrays')
        self.native_add = self.pt.add
        self.native_integrate = self.pt.integrate
        self.native_set_positions = self.pt.set_positions
        self.active = [a for a in self.pt.my_atom_indices if self.pt.get_function_count(a)]
        self.slices = {a: slice(self.pt.M_a[a],self.pt.M_a[a]+self.pt.get_function_count(a))
                       for a in self.active}
        self.ni = self.pt.Mmax
        if not self.ni or sum(self.pt.get_function_count(a) for a in self.active) != self.ni:
            raise ValueError('Incomplete serial projector ownership')
        self.fixed_spos = calc.atoms.get_scaled_positions()[self.active].copy()
        self.shape = tuple(calc.wfs.kpt_u[0].psit_nG.shape[-3:])
        self.ngrid = int(np.prod(self.shape))
        self.projectors = {}
        self.adjoints = {}
        self.scratch = {}
        coefficients = self.pt.dict(self.ni,zero=True)
        identity = np.eye(self.ni,dtype=complex)
        for atom,section in self.slices.items():
            coefficients[atom][:] = identity[:,section]
        for q in sorted({k.q for k in calc.wfs.kpt_u}):
            basis = np.zeros((self.ni,)+self.shape,dtype=complex)
            self.native_add(basis,coefficients,q)
            matrix = basis.reshape(self.ni,self.ngrid)
            self.projectors[q] = matrix
            self.adjoints[q] = matrix.conj().T

    def _arrays(self, array, q):
        if (array.dtype != np.dtype(complex) or not array.flags.c_contiguous or
                array.ndim not in (3,4) or tuple(array.shape[-3:]) != self.shape):
            raise ValueError('Unsupported frozen-projector array layout')
        if q not in self.projectors:
            raise ValueError('Bloch q differs from the qualified projector set')
        count = 1 if array.ndim==3 else array.shape[0]
        if count not in self.scratch:
            self.scratch[count] = (np.empty((count,self.ni),dtype=complex),
                                   np.empty((count,self.ngrid),dtype=complex))
        return array.reshape(count,self.ngrid),self.scratch[count]

    def add(self, a_xG, c_axi=1.0, q=-1):
        if isinstance(c_axi,float):
            raise ValueError('Complex Bloch projectors require atomic coefficients')
        flat,(packed,image) = self._arrays(a_xG,q)
        for atom,section in self.slices.items():
            packed[:,section] = c_axi[atom].reshape(len(flat),-1)
        np.dot(packed,self.projectors[q],out=image)
        flat += image

    def integrate(self, a_xG, c_axi, q=-1, add_to=False):
        if add_to:
            raise ValueError('Native projector integration does not support add_to')
        flat,(packed,_) = self._arrays(a_xG,q)
        np.dot(flat,self.adjoints[q],out=packed)
        packed *= self.pt.gd.dv
        for atom,section in self.slices.items():
            c_axi[atom][...] = packed[:,section].reshape(c_axi[atom].shape)

    def set_positions(self, spos_ac, atom_partition=None):
        if not np.allclose(np.asarray(spos_ac)[self.active],self.fixed_spos,
                           rtol=0,atol=1e-12):
            raise ValueError('A nonempty projector centre moved; the fixed cache is invalid')
        return self.native_set_positions(spos_ac,atom_partition)

    def qualify(self, calc):
        checks=[]
        for point in calc.wfs.kpt_u:
            psi=point.psit_nG
            native=self.pt.dict(len(psi))
            dense=self.pt.dict(len(psi))
            t0=perf_counter(); self.native_integrate(psi,native,point.q)
            integrate_native_s=perf_counter()-t0
            t0=perf_counter(); self.integrate(psi,dense,point.q)
            integrate_dense_s=perf_counter()-t0
            integral_error=max(np.max(np.abs(native[a]-dense[a])) for a in self.active)
            expected=np.zeros_like(psi); actual=np.zeros_like(psi)
            t0=perf_counter(); self.native_add(expected,native,point.q)
            add_native_s=perf_counter()-t0
            t0=perf_counter(); self.add(actual,native,point.q)
            add_dense_s=perf_counter()-t0
            add_error=float(np.max(np.abs(expected-actual)))
            if max(integral_error,add_error)>1e-11:
                raise ValueError(f'Native projector contraction mismatch at q={point.q}: '
                                 f'{integral_error}, {add_error}')
            checks.append(dict(q=int(point.q),integration_max_abs=float(integral_error),
                               addition_max_abs=add_error,native_integrate_s=integrate_native_s,
                               dense_integrate_s=integrate_dense_s,native_add_s=add_native_s,
                               dense_add_s=add_dense_s))
        wrong=calc.atoms.get_scaled_positions()
        wrong[self.active[0],0] += .01
        try:
            self.set_positions(wrong)
        except ValueError as error:
            if 'nonempty projector centre moved' not in str(error):
                raise
            refused=True
        else:
            raise AssertionError('A changed projector centre was accepted')
        return dict(nonlocal_projectors=int(self.ni),grid_shape=list(self.shape),
                    cached_bytes=sum(p.nbytes+self.adjoints[q].nbytes
                                     for q,p in self.projectors.items()),
                    actual_native_checks=checks,moved_projector_refused=refused,
                    scope='Exact contractions on actual saved-state orbitals; not material convergence')


def install_fixed_projectors(calc):
    """Qualify against actual Native operators before installing case-local methods."""
    cache=_FixedProjectors(calc)
    report=cache.qualify(calc)
    cache.pt.add=cache.add
    cache.pt.integrate=cache.integrate
    cache.pt.set_positions=cache.set_positions
    return report
