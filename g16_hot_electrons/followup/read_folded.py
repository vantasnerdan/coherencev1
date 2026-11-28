"""Rebuild the folded PdD force projections from qualified Gamma and X points."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.units import kB

from fold import SUPERCELL, geometry, mesh
from native import BANDS, CONVERGENCE, CUTOFF_EV, PAW, require, write_once
from read_crossover import ITEM1_SHA

T = (300, 3000, 10000)


def load(root, rel):
    data = (root / rel).read_bytes()
    return json.loads(data), dict(path=rel.as_posix(), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def validate(row, producer_sha, n, t, mode, axis, u):
    a, sites, phases, unit = geometry(mode, axis)
    a.positions[sites] += (u*phases)[:, None]*unit
    spec = row['specification']
    require(row['schema'] == 'g16-pdd-folded-point-v1' and row['displacement_A'] == u
            and spec['producer_sha256'] == producer_sha and spec['temperature_K'] == t
            and abs(spec['smearing_eV']-kB*t) < 1e-14 and spec['mode'] == mode and spec['polarization'] == axis
            and spec['primitive_gamma_mesh'] == [n]*3 and spec['primitive_to_four_atom_integer_rows'] == SUPERCELL.tolist()
            and np.array_equal(spec['folded_electronic_kpoints_fractional'], mesh(n)), 'Native fold, point or electronic mesh differs')
    require(spec['D_indices'] == sites.tolist() and spec['D_displacement_signs'] == phases.tolist()
            and spec['polarization_cartesian'] == unit.tolist() and spec['cutoff_eV'] == CUTOFF_EV and spec['nbands'] == 2*BANDS
            and spec['xc'] == 'PBE' and spec['PAW'] == {k: dict(bytes=b,sha256=h) for k,(b,h) in PAW.items()}
            and spec['net_charge'] == 0 and spec['spinpol'] is False and spec['convergence'] == CONVERGENCE,
            'Native mode or calculator differs')
    require(row['symbols'] == ['Pd','Pd','H','H'] and np.array_equal(row['cell_A'], a.cell.array)
            and np.array_equal(row['masses_u'], a.get_masses()) and abs(row['volume_A3']-a.get_volume()) < 1e-12
            and np.allclose(row['positions_A'], a.positions, atol=1e-12, rtol=0), 'Native ionic coordinates or masses differ')
    f, ns, e, ex = (row[k] for k in ('energy_free_eV','minus_T_S_eV','energy_internal_eV','energy_extrapolated_eV'))
    forces = np.asarray(row['forces_eV_A'])
    count = n**3//2
    require(np.isfinite([f,ns,e,ex]).all() and ns <= 0 and abs(e-f+ns)<1e-10 and abs(ex-f+ns/2)<1e-9
            and forces.shape == (4,3) and np.isfinite(forces).all() and abs(row['electron_count']-22)<1e-8
            and 0 <= row['maximum_last_band_occupation'] <= 1e-6 and row['bz_kpoints_count'] == count
            and 0 < row['ibz_kpoints_count'] <= count and 0 < row['scf_iterations'] <= 120, 'Native finite-T state or forces differ')
    return forces, sites, phases, unit


def projected(zero, minus, plus, forces_minus, forces_plus, sites, phases, unit, u):
    m = forces_minus[sites] @ unit
    p = forces_plus[sites] @ unit
    per_site = (m-p)*phases/(2*u)
    energy_k = (minus['energy_free_eV']+plus['energy_free_eV']-2*zero['energy_free_eV'])/(2*u*u)
    return dict(K_force_eV_A2=float(np.mean(per_site)), individual_D_projected_K_eV_A2=per_site.tolist(),
                K_energy_eV_A2=energy_k, folded_force_energy_difference_eV_A2=abs(float(np.mean(per_site))-energy_k),
                D_site_projection_difference_eV_A2=float(np.ptp(per_site)),
                inversion_force_residual_eV_A=float(np.linalg.norm(forces_minus+forces_plus,axis=1).max()),
                total_force_residual_eV_A=float(max(np.linalg.norm(forces_minus.sum(axis=0)),np.linalg.norm(forces_plus.sum(axis=0)))))


class FoldedNative:
    def __init__(self, root, producer):
        self.root = root
        self.producer_sha = hashlib.sha256(producer.read_bytes()).hexdigest()
        self.sources = []

    def point(self, n, t, mode, polarization, u):
        label = 'zero' if u == 0 else f'{mode}-{polarization}'
        path = Path(f'k{n}') / str(t) / label / f'u{u:+.4f}.json'
        row, file = load(self.root, path)
        f, sites, phases, unit = validate(row, self.producer_sha, n, t, mode, polarization, u)
        self.sources.append(file)
        return row, f, sites, phases, unit

    def pair(self, zero, n, t, mode, axis):
        m, fm, sites, phases, unit = self.point(n,t,mode,axis,-0.01)
        p, fp, _, _, _ = self.point(n,t,mode,axis,0.01)
        return projected(zero,m,p,fm,fp,sites,phases,unit,0.01)


def folded_combo(gamma, boundary):
    return dict(K_diagonal_folded_eV_A2=(boundary+gamma)/2, J_cross_folded_eV_A2=(boundary-gamma)/2,
                J_shell1_plus_shell3_scalar_projection_eV_A2=0.75*(boundary-gamma))


def build(root, producer, item1, check_only):
    raw = item1.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == ITEM1_SHA, 'Wrong accepted primitive counterpart')
    old = {r['temperature_K']:r for r in json.loads(raw)['rows']}
    native = FoldedNative(root, producer)
    outputs = []
    bases = {}
    for t,n in ((300,12),) if check_only else ((t,n) for t in T for n in (12,16)):
        z = native.point(n,t,'gamma','100',0.)[0]
        row = dict(temperature_K=t, primitive_mesh_n=n, U_eV_per_D_site=z['energy_internal_eV']/2)
        if n not in bases:
            bases[n] = row['U_eV_per_D_site']
        row['electronic_excitation_eV_per_D_site'] = row['U_eV_per_D_site']-bases[n]
        g = native.pair(z,n,t,'gamma','100')
        reference = old[t][f'mesh{n}']['harmonic']['K_force_eV_A2']
        limit = max(.002,.01*abs(reference))
        diff = abs(g['K_force_eV_A2']-reference)
        row.update(gamma_100=g, accepted_primitive_K_eV_A2=reference, primitive_fold_difference_eV_A2=diff,
                   gamma_matches_primitive=diff<=limit and g['folded_force_energy_difference_eV_A2']<=limit,
                   gamma_tolerance_eV_A2=limit)
        if not check_only:
            for axis in ('100','010'):
                x = native.pair(z,n,t,'x',axis)
                row[f'x_{axis}'] = x
                if row['gamma_matches_primitive']:
                    row[f'folded_split_{axis}'] = folded_combo(g['K_force_eV_A2'],x['K_force_eV_A2'])
        outputs.append(row)
    return dict(schema='g16-pdd-folded-mode-v1', item1_response_sha256=ITEM1_SHA, rows=outputs, sources=native.sources)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native-root', type=Path, required=True)
    p.add_argument('--item1', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--check-only', action='store_true')
    a = p.parse_args()
    write_once(a.output, build(a.native_root, Path(__file__).with_name('folded_native.py'), a.item1, a.check_only))


if __name__ == '__main__':
    main()
