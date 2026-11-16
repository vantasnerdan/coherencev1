"""Saved-state DOS, spin/occupation normalization, and explicit PAW projections."""
import numpy as np
from gpaw import GPAW
from gpaw.dos import DOSCalculator
from gpaw.mpi import serial_comm
from scipy.special import expit

from screen_io import require, write_npz


def moments(energies, projection):
    occupied = float(np.trapezoid(projection, energies))
    require(occupied > 0 and np.min(projection) >= -1e-12,
            'Invalid occupied metal d projection')
    centre = float(np.trapezoid(energies*projection, energies)/occupied)
    width = float(np.sqrt(np.trapezoid((energies-centre)**2*projection,
                                     energies)/occupied))
    return dict(occupied_d_states_per_cell=occupied,
                occupied_d_centre_relative_EF_eV=centre,
                occupied_d_rms_width_eV=width)


def analyse_dos(gpw, atoms, config, output, electrons):
    # Isolate the DOS reader: it shifts its eigenvalue array relative to EF.
    calc = GPAW(str(gpw), legacy_gpaw=False, communicator=serial_comm, txt=None)
    dos = DOSCalculator.from_calculator(calc)
    require(dos.collinear, 'Only collinear charge response is in this model')
    settings = config['ground_state']
    count = float(dos.degeneracy*np.einsum(
        'k,skn->', dos.weight_k,
        expit(-dos.eig_skn/settings['occupation_width_eV'])))
    require(abs(count-electrons) < 1e-5,
            'Saved-band occupations do not reproduce the electron count')
    highest_band_min = float(np.min(dos.eig_skn[:, :, -1]))
    if config['kind'] == 'dos':
        require(highest_band_min > 1.,
                'DOS-only band ceiling is too close to EF; omitted crossings are not excluded')
    n_ef = {}
    for width in config['dos']['widths_eV']:
        value = float(dos.raw_dos([0.0], width=width)[0])
        independent = float(dos.degeneracy*np.einsum(
            'k,skn->', dos.weight_k, np.exp(-(dos.eig_skn/width)**2))
            / (np.sqrt(np.pi)*width))
        require(np.isclose(value, independent, rtol=1e-10, atol=1e-12),
                'DOS spin, energy, or k-weight normalization mismatch')
        n_ef[str(width)] = value
    energies = np.arange(*config['dos']['occupied_window_eV'])
    total_d = np.zeros_like(energies)
    by_element = {}
    for a in config['dos']['metal_indices']:
        symbol = atoms[a].symbol
        if symbol not in by_element:
            by_element[symbol] = np.zeros_like(energies)
        atomic_d = dos.raw_pdos(energies, a=a, l=2,
                               width=config['dos']['moment_width_eV'])
        total_d += atomic_d
        by_element[symbol] += atomic_d
    hydrogen = [a for a, atom in enumerate(atoms) if atom.symbol == 'H']
    channels = {str(a): sorted(set(int(l) for l in dos.setups[a].l_j))
                for a in hydrogen}
    d_ef = None
    if hydrogen:
        d_ef = sum(float(dos.raw_pdos([0.0], a=a, l=l, width=0.0)[0])
                   for a in hydrogen for l in channels[str(a)])
    return dict(
        settings=config['dos'], fermi_energy_eV=float(dos.fermi_level),
        nspins=dos.nspins, spin_degeneracy=dos.degeneracy,
        highest_saved_band_min_relative_EF_eV=highest_band_min,
        n_EF_states_per_eV_per_cell=n_ef,
        n_EF_tetrahedron_states_per_eV_per_cell=float(dos.raw_dos([0.0], width=0.0)[0]),
        **moments(energies, total_d),
        by_metal_element={s: moments(energies, projection)
                          for s, projection in by_element.items()},
        D_projection_present=bool(hydrogen), D_projector_channels=channels,
        D_projected_n_EF_tetrahedron_states_per_eV_per_cell=d_ef,
        D_projected_n_EF_tetrahedron_states_per_eV_per_D=(
            d_ef/len(hydrogen) if hydrogen else None),
        numerical_electron_count=count, eligible=False,
        reason='Requires independent observable-level convergence qualification',
        spectrum=write_npz(output/'dos.npz', energies_eV=energies,
                           occupied_metal_d_states_per_eV=total_d,
                           **{'occupied_'+s+'_d_states_per_eV': p
                              for s, p in by_element.items()}))
