"""One unpublished pulse/size capsule, exact safe extraction and saved consumer."""
from pathlib import Path, PurePosixPath
import gzip
import hashlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREFIX = HERE.relative_to(ROOT)
OWN = ('README.md', 'method.md', 'inward.py', 'inward_consumer.py',
       'inward_probe.json', 'pulse_size_consumer.py', 'pulse_size_package.py',
       'inward_N131072.json', 'inward_N131072.npz',
       'inward_N16384.json', 'inward_N16384.npz')
CONTACT = ('contact.py', 'g12_readout.py', 'inputs.json', 'results.json',
           'readout_contract.md', 'method.md')

def package():
    payload, roles = {}, {}
    def add(path, role):
        relative = str(path.relative_to(ROOT))
        payload[relative], roles[relative] = path.read_bytes(), role
    baseline = json.loads((HERE/'nonlinear_manifest.json').read_text())
    for item in baseline['files']:
        path = ROOT/item['path']
        content = path.read_bytes()
        if path not in (HERE/'README.md', HERE/'method.md'):
            assert len(content) == item['bytes'] and hashlib.sha256(content).hexdigest() == item['sha256']
        add(path, 'unchanged exact nonlinear replay binding')
    for name in OWN:
        add(HERE/name, 'pulse/size delivery')
    for name in CONTACT:
        add(HERE.parent/'contact'/name, 'unchanged conditional radial model binding')
    add(HERE/'continuation_traces.npz', 'unchanged exact accepted seed-restoration binding')
    add(HERE/'nonlinear_manifest.json', 'unchanged original nonlinear replay manifest')
    for path in sorted((HERE/'pulse').iterdir()):
        if path.is_file() and path.suffix in ('.py', '.json', '.npz', '.md'):
            add(path, 'integrated conditional pulse result/source')
    manifest = dict(schema='pdd-pulse-size-capsule-v1', publication='unpublished local only',
        accepted_nonlinear_revision='49b8831d5f5e5ed3b7f2e72bc0e7fa536765e25a',
        plan='coherencev1-g11-pulse-size-20261125',
        files=[dict(path=p, bytes=len(b), sha256=hashlib.sha256(b).hexdigest(), role=roles[p])
               for p, b in sorted(payload.items())],
        semantics='strictly inward classical selection and explicitly conditional radial pulse probabilities; no physical foil rate',
        archive='pulse_size.tar.gz', archive_hash='reported separately; no circular identity')
    encoded = (json.dumps(manifest, indent=2)+'\n').encode()
    (HERE/'pulse_size_manifest.json').write_bytes(encoded)
    payload[str(PREFIX/'pulse_size_manifest.json')] = encoded
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode='w') as archive:
        for path, content in sorted(payload.items()):
            entry = tarfile.TarInfo(path)
            entry.size, entry.mode, entry.mtime = len(content), 0o644, 0
            archive.addfile(entry, io.BytesIO(content))
    compressed = gzip.compress(raw.getvalue(), compresslevel=9, mtime=0)
    (HERE/'pulse_size.tar.gz').write_bytes(compressed)
    with tempfile.TemporaryDirectory(prefix='pdd-pulse-size-receiver-') as temp:
        receiver = Path(temp)
        with tarfile.open(fileobj=io.BytesIO(compressed), mode='r:gz') as archive:
            for item in archive.getmembers():
                path = PurePosixPath(item.name)
                assert item.isfile() and not path.is_absolute() and '..' not in path.parts
            archive.extractall(receiver, filter='data')
        assert all((receiver/path).read_bytes() == content for path, content in payload.items())
        subprocess.run([sys.executable, 'pulse_size_consumer.py'], cwd=receiver/PREFIX, check=True)
    identity = dict(path=str(PREFIX/'pulse_size.tar.gz'), bytes=len(compressed),
        sha256=hashlib.sha256(compressed).hexdigest(), members=len(payload),
        verdict='safe fresh extraction; exact archived bytes; settled integrated saved consumer passed',
        publication='none; unpublished local capsule')
    (HERE/'pulse_size_archive_identity.json').write_text(json.dumps(identity, indent=2)+'\n')
    print(json.dumps(identity, indent=2), flush=True)


if __name__ == '__main__':
    package()
