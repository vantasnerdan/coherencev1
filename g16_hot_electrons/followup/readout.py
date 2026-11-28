"""Rebuild the G16 native material response without rerunning GPAW."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.units import GPa, kB

from native import AXES, BANDS, CONVERGENCE, PAW, TEMPERATURES_K, reference, require, write_once
from thermal_response import mesh_check, quartic, signed_pair, thermal_delta


class NativeTable:
    def __init__(self, root, producer):
        self.root = root
        self.producer_sha256 = hashlib.sha256(producer.read_bytes()).hexdigest()
        self.atoms = reference()
        self.sources = {}

    def get(self, mesh, temp, folder, direction, u):
        rel = Path(f'k{mesh}') / str(temp) / folder / f'u{u:+.4f}.json'
        data = (self.root / rel).read_bytes()
        row = json.loads(data)
        self.verify(row, mesh, temp, direction, u)
        self.sources[rel.as_posix()] = dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
        return row

    def verify(self, row, mesh, temp, direction, u):
        spec = row['specification']
        expected_paw = {name: dict(bytes=n, sha256=d) for name, (n, d) in PAW.items()}
        require(row['schema'] == 'g16-pdd-thermal-point-v1' and spec['producer_sha256'] == self.producer_sha256
                and (row['direction'], row['displacement_A']) == (direction, u)
                and (spec['temperature_K'], spec['kmesh']) == (temp, [mesh] * 3)
                and abs(spec['smearing_eV'] - kB * temp) < 1e-14 and spec['PAW'] == expected_paw
                and spec['cutoff_eV'] == 500 and spec['nbands'] == BANDS and spec['gamma'] is True and spec['charge'] == 0
                and spec['spinpol'] is False and spec['xc'] == 'PBE' and spec['convergence'] == CONVERGENCE,
                'Native lattice point identity or electronic occupation differs')
        xyz = self.atoms.positions.copy()
        if direction != 'zero':
            xyz[1] += u * np.asarray(AXES[direction])
        require(np.array_equal(self.atoms.cell.array, row['cell_A']) and np.allclose(xyz, row['positions_A'], rtol=0, atol=1e-12)
                and self.atoms.get_chemical_symbols() == row['symbols'] and self.atoms.get_masses().tolist() == row['masses_u']
                and abs(row['volume_A3'] - self.atoms.get_volume()) < 1e-12, 'Native coordinates or volume differ')
        f, ns, e, ex = (row[k] for k in ('energy_free_eV', 'minus_T_S_eV', 'energy_internal_eV', 'energy_extrapolated_eV'))
        force = np.asarray(row['forces_eV_A'])
        require(np.isfinite([f, ns, e, ex]).all() and force.shape == (2, 3) and np.isfinite(force).all()
                and ns <= 0 and abs(e - (f - ns)) < 1e-10 and abs(ex - (f - ns/2)) < 1e-9
                and abs(row['electron_count'] - 11) < 1e-8 and 0 <= row['maximum_last_band_occupation'] <= 1e-6
                and 0 < row['n_ibz_kpoints'] <= mesh**3 and 0 < row['scf_iterations'] <= 120,
                'Native energy, force or finite-band electronic closure differs')
        if direction == 'zero':
            stress = np.asarray(row['stress_eV_A3'])
            require(stress.shape == (3, 3) and np.isfinite(stress).all() and
                    abs(row['pressure_GPa'] + np.trace(stress) / (3 * GPa)) < 1e-10, 'Native pressure sign or units differ')

    def pair(self, zero, mesh, temp, folder, direction, u):
        minus, plus = (self.get(mesh, temp, folder, direction, s * u) for s in (-1, 1))
        return signed_pair(zero, minus, plus, u, AXES[direction])

    def one_temperature(self, mesh, temp, reference_wave):
        zero = self.get(mesh, temp, 'zero', 'zero', 0.)
        harm = self.pair(zero, mesh, temp, '100-fd', '100', 0.01)
        extra = self.pair(zero, mesh, temp, '100-fd', '100', 0.02) if temp in (300, 350, 30000) else None
        centered = {k: zero[k] for k in ('volume_A3', 'energy_internal_eV', 'energy_free_eV', 'minus_T_S_eV',
                                         'pressure_GPa', 'stress_eV_A3', 'forces_eV_A')}
        quartics = {}
        if not reference_wave and mesh == 12:
            for direction in ('100', '111'):
                samples = [self.pair(zero, mesh, temp, f'{direction}-quartic', direction, u) for u in (0.12, 0.24)]
                quartics[direction] = dict(samples=samples, fit=quartic(samples, harm['K_force_eV_A2']))
        return dict(zero=centered, harmonic=harm, step_control=extra, quartic=quartics)


def build(root, producer, reference_wave):
    reader = NativeTable(root, producer)
    temperatures = (300, 350) if reference_wave else TEMPERATURES_K
    rows = {t: {k: reader.one_temperature(k, t, reference_wave) for k in (12, 16)} for t in temperatures}
    for t, values in rows.items():
        for k, row in values.items():
            row['contrast'] = thermal_delta(row, rows[300][k])
            step = row['step_control']
            if step:
                alt = step['K_force_eV_A2'] - rows[300][k]['step_control']['K_force_eV_A2']
                row['contrast']['alternate_step_K_change_eV_A2'] = alt
                row['contrast']['delta_K_step_sensitivity_eV_A2'] = abs(alt - row['contrast']['K_change_eV_A2'])
    return dict(schema='g16-pdd-thermal-material-v1', reference_temperature_K=300,
                rows=[dict(temperature_K=t, mesh12=v[12], mesh16=v[16], mesh_check=mesh_check(v[12], v[16])) for t,v in rows.items()],
                native_sources=reader.sources)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--producer', type=Path, default=Path(__file__).with_name('native.py'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--reference-wave', action='store_true')
    args = p.parse_args()
    write_once(args.output, build(args.root, args.producer, args.reference_wave))


if __name__ == '__main__':
    main()
