"""Build one unpublished, deterministic, whitelist-only focusing capsule."""
from pathlib import Path, PurePosixPath
import gzip
import hashlib
import io
import json
import tarfile
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FILES=(
    'README.md','method.md','keep-drop.md','patterns.md','rate_mapping.md',
    'engine.py','patterns.py','rate_mapping.py','run.py','analyze.py','verify.py','package.py',
    'results.json','comparison.csv','trajectories.npz','launches.npz','consumer_proof.json','focusing.png',
    'input_manifest.json')
INPUTS=(
    'g10_quantum_site/model.py',
    'g3_loading/x1/dense_Pd/k3/final-x1-dense-readout-20261114T0213Z.json',
    'g10_quantum_site/bulk-refit-080/well.py',
    'g10_quantum_site/bulk-refit-048/well.py',
    'g10_quantum_site/bulk-refit-048/inputs/old-readout.json',
    'g10_quantum_site/bulk-refit-048/inputs/bulk-wells.json',
    'g10_quantum_site/.input-080/g10_quantum_site/bulk-wide-wells/bulk-wells.json',
    'g11_cmns/collective/.pair-input/g10_quantum_site/pdd-nearest-pair/pair-cut.csv',
    'g11_cmns/collective/contact/results.json',
    'g11_cmns/collective/contact/g12_readout.py')


def record(path,data):
    return dict(path=path,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def package():
    inputs=[record(p,(ROOT/p).read_bytes()) for p in INPUTS]
    # Exact byte equality, rather than matching filenames, supplies these joins.
    equivalences=[]
    for a,b in (
        (INPUTS[1],INPUTS[4]),
        (INPUTS[7],'g11_cmns/collective/contact/sources/pair-cut.csv')):
        aa,bb=(ROOT/a).read_bytes(),(ROOT/b).read_bytes()
        equivalences.append(dict(first=a,second=b,byte_equal=aa==bb,sha256_first=hashlib.sha256(aa).hexdigest(),sha256_second=hashlib.sha256(bb).hexdigest()))
    main_equivalences=[
        dict(local=INPUTS[6],canonical='g10_quantum_site/bulk-wide-wells/bulk-wells.json',bytes=319722,sha256='18fa6232525243391b84f2fb495cfe9947f7d41eef270931fb3dc8b99857e28d'),
        dict(local=INPUTS[7],canonical='g10_quantum_site/pdd-nearest-pair/pair-cut.csv',bytes=844,sha256='8b610c379ab585b01a2a17b140de97835d0733fe29a09d5bfe2e350f3b30f8b6'),
        dict(local=INPUTS[2],canonical=INPUTS[2],bytes=2881,sha256='a791fe4a4b13e5cf60ab0fc1e20c8c4c8be2e49a67ad2b2e0b373465a314bc9e')]
    for item in main_equivalences:
        item.update(revision='c6f079ce30333158160aed50a391f0cae6ccc0f0',
                    verification='GitHub contents API bytes decoded and compared directly to this local input; byte-equal')
        assert record(item['local'],(ROOT/item['local']).read_bytes())['sha256']==item['sha256']
    manifest=dict(schema='focusing-assumed-inputs-v1',inputs=inputs,checked_equivalences=equivalences,verified_main_equivalences=main_equivalences,
        source_revisions={'workspace_input_HEAD':'7349c26e626d41262f9ad66ab94e9e266373425c','accepted_080':'e07c49e','pair_source':'e8f7286','accepted_contact':'c6f079ce30333158160aed50a391f0cae6ccc0f0'},
        pair_asset={'id':790324076,'archive_sha256':'716978d9797e12683a937a5f1fe79c715b4d3d000e8ac4da5614cc2f25ef7777'},
        assumptions='Extract at g11_cmns/collective/focusing/ in the adopted source tree with the exact hashed prior inputs. Prior arrays are not duplicated. No contact rerun. Revision labels are authoritative supplied provenance; exact current input bytes are recorded above.')
    (HERE/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    payload={p:(HERE/p).read_bytes() for p in FILES}
    records=[record(p,b) for p,b in payload.items()]
    whitelist=''.join(p+'\n' for p in FILES)+ 'ARCHIVE_FILES.txt\nSHA256SUMS\n'
    sums=''.join(f"{r['sha256']}  {r['path']}\n" for r in records)
    (HERE/'ARCHIVE_FILES.txt').write_text(whitelist)
    (HERE/'SHA256SUMS').write_text(sums)
    payload['ARCHIVE_FILES.txt']=whitelist.encode();payload['SHA256SUMS']=sums.encode()
    raw=io.BytesIO()
    with tarfile.open(fileobj=raw,mode='w') as tf:
        for name,data in payload.items():
            info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644;info.mtime=0
            tf.addfile(info,io.BytesIO(data))
    archive=gzip.compress(raw.getvalue(),compresslevel=9,mtime=0)
    (HERE/'focusing.tar.gz').write_bytes(archive)
    with tempfile.TemporaryDirectory(prefix='focusing-receiver-') as directory:
        destination=Path(directory)
        with tarfile.open(fileobj=io.BytesIO(archive),mode='r:gz') as tf:
            for member in tf.getmembers():
                name=PurePosixPath(member.name)
                if name.is_absolute() or '..' in name.parts or not member.isfile() or member.name not in payload:
                    raise ValueError('unsafe or unlisted archive member')
                data=tf.extractfile(member).read()
                (destination/member.name).write_bytes(data)
        for name,data in payload.items():
            assert (destination/name).read_bytes()==data
    capsule=dict(schema='g11-focusing-capsule-v1',archive='focusing.tar.gz',bytes=len(archive),sha256=hashlib.sha256(archive).hexdigest(),
        whitelist='ARCHIVE_FILES.txt',files=records,
        verification='Fresh-destination safe-member extraction and every-member byte equality; physics consumer_proof.json is separate.',
        plan_id='coherencev1-g11-pdd-focusing-20261125',publication='Unpublished local staged bytes only; parent owns handoff and any later publication.')
    (HERE/'capsule.json').write_text(json.dumps(capsule,indent=2)+'\n')
    print(json.dumps({k:v for k,v in capsule.items() if k!='files'},indent=2))

if __name__=='__main__':package()
