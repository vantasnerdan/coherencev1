"""Byte-bound primitive records; no scientific defaults in the consumer."""
import hashlib
import json
import os
from pathlib import Path


def manifest(path):
    p = Path(path)
    raw = p.read_bytes()
    return {'path': str(p), 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def save(path, value):
    p = Path(path)
    tmp = p.with_name(p.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(tmp, p)


def geometry(atoms):
    return {'symbols': atoms.get_chemical_symbols(), 'positions_A': atoms.positions.tolist(),
            'cell_A': atoms.cell.tolist(), 'pbc': atoms.pbc.tolist(),
            'constraints': [c.todict() for c in atoms.constraints]}
