"""Pass 5 selected-input ownership and units; never invoke baseline producers."""
from pathlib import Path
import csv, hashlib, importlib.util, json, sys
import numpy as np
HERE = Path(__file__).resolve().parent
BASE = 'f6010129f11ba8010815b78a67b3b6bed1db7194'
QE = 1.602176634e-19
MD = 2.0141 * 1.66053906660e-27
KB = 8.617333262145e-5
L = np.sqrt(np.pi) * .012
H = .00025
RHO = 12007.
ENERGIES = (.25, .5, 1., 1.5, 2., 2.5)
BEAMS = ('even_steady', 'hot_steady', 'hot_bursts')

def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def emit(path, value):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False, default=convert) + '\n')

def table(path, rows):
    if not rows: raise ValueError('Empty scientific table: ' + str(path))
    names = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader(); writer.writerows(rows)

def inherited(path4):
    """One exact baseline import epoch; all baseline output functions stay unused."""
    path4 = Path(path4).resolve()
    if str(path4) not in sys.path: sys.path.insert(0, str(path4))
    import sources, native, fields, launch, response, acoustic, yielding
    sources.HERE = path4; native.HERE = path4
    return sources, native, fields, launch, response, acoustic, yielding

def calendar(path4):
    """Consume the frozen Cox edge/phase bytes, never regenerate the rain."""
    fields=inherited(path4)[2]
    result=fields.Calendar.__new__(fields.Calendar)
    with np.load(Path(path4)/'calendar.npz') as saved:
        result.edges=saved['edges_s'];result.state=saved['burst_on']
        result.phase=float(saved['phase']);result.norm=float(saved['normalization'])
        result.seconds=float(saved['seconds'])
    result.omega=2*np.pi*1e5
    result.burst_integral=np.r_[0.,np.cumsum(np.diff(result.edges)*result.state)]
    return result

def verify_selected(path4, names):
    path4 = Path(path4)
    expected = {}
    for name in ('manifest_core.json', 'manifest_evidence.json'):
        for row in json.loads((path4 / name).read_text())['files']:
            expected[row['path'].removeprefix('g15_track/pass4/')] = row
    rows = []
    for name in sorted(set(names)):
        row = expected.get(name)
        if row is None: raise ValueError('Input absent from accepted manifests: ' + name)
        path = path4 / name
        if path.stat().st_size != row['bytes'] or digest(path) != row['sha256']:
            raise ValueError('Accepted input changed: ' + name)
        rows.append(dict(path=name, bytes=row['bytes'], sha256=row['sha256']))
    return rows
