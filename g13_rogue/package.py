"""Package exactly g13_rogue/, without rerunning scientific producers."""
from pathlib import Path
import gzip
import hashlib
import io
import json
import tarfile

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'g13_rogue'
EXCLUDE = {'capsule.tar.gz','manifest.json','SHA256SUMS','delivery.json'}


def main():
    files = []
    for path in sorted(HERE.rglob('*')):
        if path.is_symlink():
            raise ValueError('Capsule does not accept symlinks')
        if path.is_file() and path.name not in EXCLUDE and '__pycache__' not in path.parts:
            files.append(path)
    rows=[]
    payloads={}
    for path in files:
        data=path.read_bytes()
        name=str(path.relative_to(ROOT))
        payloads[name]=data
        rows.append(dict(path=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    manifest=dict(schema='g13-rogue-pass1-capsule-v1',
        patch_base='02fe7cc19a121706b70ca201195b36fc779e0fb8',
        patch_kind='new directory overlay; no source outside g13_rogue/',
        files=rows,self_excluded=['g13_rogue/manifest.json','g13_rogue/SHA256SUMS'],
        external_delivery_record='g13_rogue/delivery.json is not an archive member; authenticated release readback binds archive bytes')
    manifest_bytes=(json.dumps(manifest,indent=2)+'\n').encode()
    sums=(''.join(f"{r['sha256']}  {r['path']}\n" for r in rows)).encode()
    (HERE/'manifest.json').write_bytes(manifest_bytes)
    (HERE/'SHA256SUMS').write_bytes(sums)
    payloads['g13_rogue/manifest.json']=manifest_bytes
    payloads['g13_rogue/SHA256SUMS']=sums
    packed=io.BytesIO()
    with gzip.GzipFile(fileobj=packed,mode='wb',mtime=0,filename='') as gz:
        with tarfile.open(fileobj=gz,mode='w') as tar:
            for name,data in sorted(payloads.items()):
                info=tarfile.TarInfo(name)
                info.size=len(data)
                info.mode=0o644
                info.mtime=0
                tar.addfile(info,io.BytesIO(data))
    data=packed.getvalue()
    path=HERE/'capsule.tar.gz'
    path.write_bytes(data)
    result=dict(path='g13_rogue/capsule.tar.gz',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),
                patch_base=manifest['patch_base'],members=len(payloads))
    print(json.dumps(result))


if __name__=='__main__':
    main()
