"""Freeze exactly tid2-site/; do not rerun scientific producers or mutate releases."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BASE = '463910dcf13ce69d74b138a21d53183129a779de'
EXCLUDE = {'capsule.tar.gz', 'manifest.json', 'SHA256SUMS', 'delivery.json', 'plan-receipt.json'}


def main():
    payloads = {}
    for path in sorted(HERE.rglob('*')):
        if path.is_symlink():
            raise ValueError('Symlinks are not capsule members')
        if not path.is_file() or path.name in EXCLUDE or '__pycache__' in path.parts:
            continue
        if path.name in ('SUMMARY.md', 'LEAD.md'):
            raise ValueError('Lead pages are outside the capsule contract')
        data = path.read_bytes()
        if path.suffix in ('.py', '.md', '.csv', '.json'):
            text = data.decode('utf-8')
            if '/'+'home'+'/' in text or '<PRIVATE'+'_'+'PATH>' in text:
                raise ValueError('Private path in '+str(path.relative_to(ROOT)))
        payloads[str(path.relative_to(ROOT))] = data
    rows = [dict(path=name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest()) for name, data in payloads.items()]
    manifest = dict(schema='tid2-site-capsule-v1', accepted_base=BASE,
                    scope='New directory overlay exactly g10_quantum_site/tid2-site/; no other campaign modified',
                    source_semantics='frozen constrained fully occupied cubic fluorite; modeled 3D angular/exterior closure; no foil phase or surface assignment',
                    scientific_outputs='H/D/T at 77/300K; signed marginals and Fourier conditional G8 measures; [100] transitions; actual G11 frozen nearest-pair comparison with reused PdD evidence',
                    files=rows, self_excluded=['g10_quantum_site/tid2-site/manifest.json', 'g10_quantum_site/tid2-site/SHA256SUMS'],
                    external_receipts='delivery.json and plan-receipt.json bind later transport and shared-plan closure without recursive archive hashes')
    manifest_bytes = (json.dumps(manifest, indent=2)+'\n').encode()
    sums = ''.join(f"{r['sha256']}  {r['path']}\n" for r in rows).encode()
    (HERE/'manifest.json').write_bytes(manifest_bytes)
    (HERE/'SHA256SUMS').write_bytes(sums)
    payloads['g10_quantum_site/tid2-site/manifest.json'] = manifest_bytes
    payloads['g10_quantum_site/tid2-site/SHA256SUMS'] = sums
    archive = HERE/'capsule.tar.gz'
    with archive.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream, mode='wb', mtime=0, filename='') as compressed:
            with tarfile.open(fileobj=compressed, mode='w') as tar:
                for name, data in sorted(payloads.items()):
                    member = tarfile.TarInfo(name)
                    member.size = len(data); member.mode = 0o644; member.mtime = 0
                    tar.addfile(member, io.BytesIO(data))
    data = archive.read_bytes()
    print(json.dumps(dict(path=str(archive.relative_to(ROOT)), bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
                          members=len(payloads), accepted_base=BASE)))


if __name__ == '__main__':
    main()
