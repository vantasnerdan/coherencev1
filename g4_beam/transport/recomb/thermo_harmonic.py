"""Harmonic free energies and pair-normalized, no-recrossing TST (eV, K, s).

Mode energies are h*nu, not angular frequencies or Hessian eigenvalues. They
must come from H_cart v = omega**2 M v in the SAME constrained mobile space,
with isotope masses in M. A TS has exactly one unstable mode, excluded from
its harmonic partition function; no positive/soft modes are silently dropped.
"""
from dataclasses import dataclass
import math
from ase import units

KB = units.kB
H_EV_S = units._hplanck / units._e
D_MASS_U = 2.014101778


def positive(value, name):
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be finite and positive')


@dataclass(frozen=True)
class ModeSet:
    """A settled spectrum and its kinetic metric; signed unstable magnitude is separate.

    active_atoms and masses_u refer to identical ordered Cartesian DOFs in IS/TS.
    constraint_id binds the shared fixed atoms/cell/projection convention.
    source_id binds the settled geometry/Hessian evidence (not merely a label).
    stable_eV has 3*N values for IS, 3*N-1 for a first-order TS; unstable_eV is
    None in IS, positive |h*Im(nu)| in TS. Classification/connection evidence
    remains the producing DFT calculation's responsibility.
    """
    active_atoms: tuple[int, ...]
    masses_u: tuple[float, ...]
    constraint_id: str
    source_id: str
    stable_eV: tuple[float, ...]
    unstable_eV: float | None = None

    def __post_init__(self):
        if not self.source_id or not self.constraint_id or not self.active_atoms:
            raise ValueError('Spectrum requires geometry/Hessian and constraint bindings')
        if len(set(self.active_atoms)) != len(self.active_atoms):
            raise ValueError('Repeated active atom')
        if len(self.active_atoms) != len(self.masses_u):
            raise ValueError('One kinetic mass per active atom required')
        for mass in self.masses_u:
            positive(mass, 'mass_u')
        for energy in self.stable_eV:
            positive(energy, 'stable h*nu in eV')
        if self.unstable_eV is not None:
            positive(self.unstable_eV, 'unstable magnitude in eV')
        expected = 3 * len(self.active_atoms) - (self.unstable_eV is not None)
        if len(self.stable_eV) != expected:
            raise ValueError('Full constrained Cartesian spectrum required')


def harmonic_parts(temperature_K, stable_eV):
    """Return distinct ZPE, thermal Helmholtz, total F in eV; no electronic E."""
    positive(temperature_K, 'temperature_K')
    energies = tuple(stable_eV)
    for energy in energies:
        positive(energy, 'stable h*nu in eV')
    kt = KB * temperature_K
    zpe = math.fsum(energies) / 2
    thermal = kt * math.fsum(math.log(-math.expm1(-e / kt)) for e in energies)
    return dict(zpe_eV=zpe, thermal_F_eV=thermal, F_eV=zpe + thermal)


def pair_tst(temperature_K, electronic_barrier_eV, initial, transition,
             *, initial_basin_degeneracy, transition_path_degeneracy):
    """Return D2 events/(specified OCCUPIED unordered pair*s), kappa=1.

    Barrier is E_TS-E_IS at matched electronic settings, WITHOUT ZPE. Explicit
    degeneracies count equivalent occupied basins/paths not already in spectra;
    these are NOT the lattice coordination or the pair occupancy probability.
    """
    if not math.isfinite(electronic_barrier_eV):
        raise ValueError('Finite matched electronic barrier required')
    positive(initial_basin_degeneracy, 'initial_basin_degeneracy')
    positive(transition_path_degeneracy, 'transition_path_degeneracy')
    if initial.unstable_eV is not None or transition.unstable_eV is None:
        raise ValueError('Stable IS and first-order TS required')
    if (initial.active_atoms, initial.masses_u, initial.constraint_id) != (
            transition.active_atoms, transition.masses_u, transition.constraint_id):
        raise ValueError('IS and TS kinetic metric/constraint mismatch')
    fi = harmonic_parts(temperature_K, initial.stable_eV)
    ft = harmonic_parts(temperature_K, transition.stable_eV)
    kt = KB * temperature_K
    dg = electronic_barrier_eV + ft['F_eV'] - fi['F_eV']
    log_rate = (math.log(kt / H_EV_S) + math.log(transition_path_degeneracy /
                initial_basin_degeneracy) - dg / kt)
    return dict(electronic_barrier_eV=electronic_barrier_eV,
                delta_ZPE_eV=ft['zpe_eV'] - fi['zpe_eV'],
                delta_thermal_F_eV=ft['thermal_F_eV'] - fi['thermal_F_eV'],
                activation_F_eV=dg, log_rate_per_occupied_pair_s=log_rate,
                rate_per_occupied_pair_s=math.exp(log_rate))
