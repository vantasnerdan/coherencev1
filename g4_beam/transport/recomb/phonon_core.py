"""Exact clamped (1,2,1) ASE-Phonons reduction; no force or spectrum repair.

C_N is UNWEIGHTED, displacement-row/force-column, before reciprocity averaging.
H=[[C0,C1],[C1,C0]] in cell-first active order. At q_y=0,1/2 the
orthogonal cell Fourier transform gives C0+C1,C0-C1, respectively.
The fixed bottom layers exert reactions: no ASR/Frederiksen/drift removal.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
from ase.phonons import Phonons
from hessian_modes import mass_spectrum

from ase.utils.filecache import MultiFileJSONCache
SUPERCELL = (1, 2, 1)
ROW_TOLERANCE = 0.05


def manifest(path):
    path = Path(path)
    data = path.read_bytes()
    return {'path': str(path.resolve()), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest()}


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def layout(primitive):
    symbols = primitive.get_chemical_symbols()
    if symbols not in (['Pd'] * 3, ['Pd'] * 3 + ['H']):
        raise ValueError('Primitive must be cell-first Pd3 or Pd3+D(H isotope)')
    if list(primitive.pbc) != [True, True, False] or abs(primitive.cell.volume) <= 0:
        raise ValueError('Nondegenerate xy-periodic, z-open primitive required')
    fixed = set()
    for constraint in primitive.constraints:
        if constraint.todict()['name'] != 'FixAtoms':
            raise ValueError('Only bottom-layer FixAtoms supported')
        fixed.update(int(i) for i in constraint.get_indices())
    if fixed != {0, 1}:
        raise ValueError('Exactly primitive bottom Pd0/Pd1 must be fixed')
    masses = primitive.get_masses()
    if np.any(~np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError('Positive finite actual isotope masses required')
    if len(primitive) == 4 and abs(masses[3] - 2.014101778) > 1e-6:
        raise ValueError('Adsorbate must carry explicit deuterium isotope mass')
    active = list(range(2, len(primitive)))
    full = active + [i + len(primitive) for i in active]
    return active, full


def reconstruct(blocks):
    blocks = np.asarray(blocks, dtype=float)
    if (blocks.ndim != 3 or blocks.shape[0] != 2
            or blocks.shape[1] != blocks.shape[2] or not blocks.shape[1]
            or np.any(~np.isfinite(blocks))):
        raise ValueError('Two finite square raw C_N blocks required')
    c0, c1 = blocks
    return np.block([[c0, c1], [c1, c0]])


def classify(spectrum):
    quanta = np.asarray(spectrum['signed_quanta_eV'])
    return {'modes': len(quanta), 'negative': int(np.sum(quanta < 0)),
            'exact_zero': int(np.sum(quanta == 0)),
            'positive': int(np.sum(quanta > 0)),
            'rule': 'Signed raw-eigenvalue classification; no tolerance or clipping'}


class RecordedPhonons(Phonons):
    """ASE cache lifecycle; record all raw forces, energy and optional result bank."""
    def __init__(self, *args, recorder=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.recorder = recorder

    def calculate(self, atoms_N, disp):
        # Constraint reactions must survive in the saved force bank.
        forces = atoms_N.get_forces(apply_constraint=False)
        energy = float(atoms_N.get_potential_energy())
        result = {'forces': forces, 'energy_eV': energy,
                  'positions_A': atoms_N.positions.copy()}
        if self.recorder is not None:
            result.update(self.recorder(atoms_N, disp.name))
        return result


def withheld_rows(phonons, raw, cache, full_indices):
    """Read independent image-atom +/-x calculations, NOT translated cache rows."""
    rows = []
    for atom in full_indices[len(phonons.indices):]:
        minus = np.asarray(cache[f'{atom}x-']['forces'])
        plus = np.asarray(cache[f'{atom}x+']['forces'])
        expected_shape = (2 * len(phonons.atoms), 3)
        if minus.shape != expected_shape or plus.shape != expected_shape:
            raise ValueError('Withheld forces do not match repeated geometry')
        observed = ((minus - plus)[full_indices] / (2 * phonons.delta)).ravel()
        row = 3 * full_indices.index(atom)
        difference = observed - raw[row]
        rows.append({'image_atom': atom, 'axis': 'x', 'row': row,
                     'observed_eV_A2': observed.tolist(),
                     'predicted_eV_A2': raw[row].tolist(),
                     'max_abs_discrepancy_eV_A2': float(np.max(np.abs(difference)))})
    maximum = max(r['max_abs_discrepancy_eV_A2'] for r in rows)
    return {'rows': rows, 'max_abs_discrepancy_eV_A2': maximum,
            'threshold_eV_A2': ROW_TOLERANCE,
            'accepted': bool(np.isfinite(maximum) and maximum <= ROW_TOLERANCE)}


def run_withholding(phonons, directory):
    """Only independent translated-cell x pairs; ordinary ASE cache locks."""
    _, full = layout(phonons.atoms)
    cache = MultiFileJSONCache(directory)
    removed = cache.strip_empties()
    atoms = phonons.atoms * SUPERCELL
    atoms.calc = phonons.calc
    for atom in full[len(phonons.indices):]:
        original = float(atoms.positions[atom, 0])
        for sign in (-1, 1):
            disp = phonons._disp(atom, 0, sign)
            with cache.lock(disp.name) as handle:
                if handle is None:
                    continue
                try:
                    atoms.positions[atom, 0] = original + sign * phonons.delta
                    handle.save(phonons.calculate(atoms, disp))
                finally:
                    atoms.positions[atom, 0] = original
    return cache, removed


def consume(phonons, withholding_cache):
    active, full = layout(phonons.atoms)
    if (list(phonons.indices) != active or phonons.supercell != SUPERCELL
            or phonons.offset != 0):
        raise ValueError('Exact corner-reference active primitive two-cell route required')
    phonons.read(method='standard', symmetrize=0, acoustic=False)
    raw = reconstruct(phonons.C_N)
    withholding = withheld_rows(phonons, raw, withholding_cache, full)
    if not withholding['accepted']:
        # Caller persists this evidence before refusing a consumable spectrum.
        return {'accepted': False, 'withholding': withholding,
                'raw_C_N_eV_A2': phonons.C_N.tolist(), 'raw_H_eV_A2': raw.tolist()}
    masses = (phonons.atoms * SUPERCELL).get_masses()[full]
    spectrum = mass_spectrum(raw, masses)
    sectors = []
    for qy, sign in ((0.0, 1), (0.5, -1)):
        sector = mass_spectrum(phonons.C_N[0] + sign * phonons.C_N[1],
                               phonons.atoms.get_masses()[active])
        sectors.append({'q_scaled': [0, qy, 0], 'spectrum': sector,
                        'classification': classify(sector)})
    return {'accepted': True, 'active_primitive': active, 'active_supercell': full,
            'raw_C_N_eV_A2': phonons.C_N.tolist(), 'withholding': withholding,
            'spectrum': spectrum, 'classification': classify(spectrum),
            'fourier_sectors': sectors}
