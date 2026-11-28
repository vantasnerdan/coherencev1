"""Bind saved GPAW arrays to the exact Phonopy displacement and IFC reconstruction."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from native import PAW, require
from phonon_geometry import electronic_mesh, inputs, inversion, phonon
from phonon_native import SCHEMA, profile

_SOURCE_DIR = Path(__file__).resolve().parent


def point(root, geo, case, cell, mesh, temperature, grid, source_sha, *, compact_bands=False, tight_forces=False):
    rel = Path(case) / 'point.json'
    raw = (root / rel).read_bytes()
    row = json.loads(raw)
    atoms, atom, vector = geo[case]
    spec = row['specification']
    cfg = profile(len(atoms),temperature,compact_bands=compact_bands,tight_forces=tight_forces)
    require(row['schema']==SCHEMA
            and (spec['cell'],spec['primitive_mesh_n'],spec['temperature_K'],spec['case'])==(cell,mesh,temperature,case)
            and spec['calculator']==cfg and spec['source_sha256']==source_sha
            and spec['PAW']=={name:dict(bytes=size,sha256=h) for name,(size,h) in PAW.items()}
            and spec['displaced_atom']==atom and spec['displacement_A']==(None if vector is None else vector.tolist())
            and np.array_equal(spec['electronic_grid'],grid), f'{case}: not the matched calculation')
    require(np.array_equal(row['positions_A'],atoms.positions) and np.array_equal(row['cell_A'],atoms.cell.array)
            and row['symbols']==atoms.get_chemical_symbols() and np.array_equal(row['masses_u'],atoms.get_masses()),
            f'{case}: order, geometry or masses differ')
    force = np.asarray(row['forces_eV_A'])
    require(force.shape==atoms.positions.shape and np.isfinite(force).all()
            and 0<=row['occupation_tail']<=cfg.get('occupation_tail_max',1e-6) and np.isfinite(row['energy_free_eV'])
            and abs(row['electron_count']-cfg['valence_electrons'])<1e-8, f'{case}: native force or occupation differs')
    return force, row['energy_free_eV'], dict(path=rel.as_posix(), sha256=sha256(raw).hexdigest(), bytes=len(raw))


def assemble(root, cell, mesh, temperature, producer_sources=_SOURCE_DIR, *, compact_bands=False, tight_forces=False):
    ph = phonon(cell)
    geo = inputs(ph)
    grid = electronic_mesh(mesh,cell)
    pins = {f:sha256((producer_sources/f).read_bytes()).hexdigest() for f in ('native.py','phonon_geometry.py','phonon_native.py')}

    def read_case(case):
        return point(root,geo,case,cell,mesh,temperature,grid,pins,compact_bands=compact_bands,tight_forces=tight_forces)

    zero = read_case('zero') if (root/'zero'/'point.json').is_file() else None
    ideal = geo['zero'][0]
    records = [zero[2]] if zero else []
    metadata = {}
    polarity = {}
    for species in ('pd','d'):
        fp, ep, identity = read_case(f'{species}-p')
        records.append(identity)
        _, mobile, step = geo[f'{species}-p']
        ids = inversion(ideal,mobile)
        inferred = -fp[ids]
        path = root / f'{species}-m' / 'point.json'
        ferr = energy_k = energy_diff = None
        if path.is_file():
            fm, em, identity = read_case(f'{species}-m')
            records.append(identity)
            ferr = float(np.linalg.norm(fm-inferred,axis=1).max())
            energy_diff = abs(ep-em)
            if zero:
                energy_k = (ep+em-2*zero[1]) / float(step@step)
            status = 'measured'
        else:
            fm = inferred
            status = 'inferred_by_inversion'
        odd = (fp-fm)/2
        polarity[mobile] = odd
        k = -float(odd[mobile] @ step) / float(step@step)
        metadata[species] = dict(negative_sign=status, inversion_force_residual_eV_A=ferr, sign_energy_residual_eV=energy_diff,
                                 max_raw_net_force_eV_A=float(max(np.linalg.norm(v.sum(axis=0)) for v in (fp,fm))),
                                 site_direction_force_K_eV_A2=k, site_direction_energy_K_eV_A2=energy_k)
    forces = []
    for mobile, *u in ph.displacements:
        element = {'Pd':'pd','H':'d'}[ideal.symbols[mobile]]
        s = np.sign(np.asarray(u) @ geo[f'{element}-p'][2])
        forces.append(s * polarity[mobile])
    ph.produce_force_constants(forces=forces,calculate_full_force_constants=True,fc_calculator='traditional',show_drift=False)
    return ph, dict(sources=records, displacement_diagnostics=metadata,
                    zero_max_force_eV_A=float(np.linalg.norm(zero[0],axis=1).max()) if zero else None)
