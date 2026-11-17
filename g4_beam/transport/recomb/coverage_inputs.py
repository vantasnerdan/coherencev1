"""Authenticate accepted operating points and load the original, unmodified Model."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

SOURCE_SHA256 = '383aad3068c5d5b62fc5f031ccfe8c35c1c6ab426f41daeca6c18aeba3d4a47e'


def binding(path):
    path = Path(path)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=sha256(data).hexdigest())


def load_source(source_path):
    source = Path(source_path).resolve()
    if binding(source)['sha256'] != SOURCE_SHA256:
        raise ValueError('Pinned transport source digest mismatch')
    # Construction unconditionally reads BCA tables even for zero beam. Never
    # allow the source to generate new tables (and launch multiprocessing).
    if not (source.parent / 'bca_tables.npz').is_file():
        raise FileNotFoundError('Existing pinned BCA tables required; generation forbidden')
    spec = importlib.util.spec_from_file_location('coverage_original_transport', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def operating_points(points_path):
    data = json.loads(Path(points_path).read_bytes())
    if data['transport_source_sha256'] != SOURCE_SHA256:
        raise ValueError('Operating points bind a different Model')
    rows = data['rows']
    if len(rows) != 2 or {r['loaded'] for r in rows} != {False, True}:
        raise ValueError('Exactly one accepted gas and one accepted loaded case required')
    return rows


def make_model(source, row, *, s_effective, grid=None):
    """Use public case parameters only; no patches, alloys or process overrides."""
    case = source.default_case('Pd', loaded=row['loaded'], E_cm=row['E_cm_keV'])
    case.update(T=row['T_K'], cooling=row['cooling'], beam=False, sigma_des=0.0,
                t_pre=0.0, s_front=s_effective, grid=grid,
                grid_both=bool(row['loaded']))
    if row['loaded']:
        case['x_back'] = row['x_back']
    model = source.Model(case)
    if model.p_gas != row['p_D2_Pa']:
        raise ValueError('Original Model pressure differs from accepted pressure')
    if row['loaded'] and abs(model.mu_back - row['mu_back_eV']) > 1e-12:
        raise ValueError('Original Model back fugacity differs from accepted value')
    return model
