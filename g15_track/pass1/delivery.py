"""Package and read back the authorized unpublished G15 Pass1 draft.
No Git commands, pushes, branch changes or publication.
"""
from pathlib import Path, PurePosixPath
import argparse, csv, gzip, hashlib, io, json, subprocess, tarfile
import numpy as np

HERE=Path(__file__).resolve().parent
PREFIX='g15_track/pass1/'
BASE='d8b935df173ae6cb02ffb0d689a64c55039c9a04'
REPO='vantasnerdan/coherencev1'
TAG='g15-pass1-all-track-energy-20261127'
EXCLUDE={'manifest.json','delivery.json','remote_delivery.json','capsule.tar.gz'}


def sha(data):return hashlib.sha256(data).hexdigest()


def verify(data,manifest):
    if len(data)>=50_000_000:raise ValueError('Capsule exceeds50MB')
    expected={r['path']:r for r in manifest['files']};opened=[]
    with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as archive:
        for member in archive.getmembers():
            p=PurePosixPath(member.name)
            if not member.isfile() or p.is_absolute() or '..' in p.parts:raise ValueError('Unsafe archive member')
            raw=archive.extractfile(member).read()
            if member.name==PREFIX+'manifest.json':
                if json.loads(raw)!=manifest:raise ValueError('Manifest content changed')
                continue
            row=expected.pop(member.name,None)
            if row is None or len(raw)!=row['bytes'] or sha(raw)!=row['sha256']:raise ValueError('Member identity mismatch')
            if member.name.endswith('.npz'):
                with np.load(io.BytesIO(raw),allow_pickle=False) as arrays:
                    shapes={k:list(arrays[k].shape) for k in arrays.files}
                opened.append(dict(path=member.name,numeric_members=len(shapes)))
            elif member.name.endswith('.png'):
                from PIL import Image
                image=Image.open(io.BytesIO(raw));image.verify();opened.append(dict(path=member.name,image_verified=True))
            elif member.name.endswith('.json'):json.loads(raw)
            elif member.name.endswith('.csv'):
                rows=list(csv.reader(io.StringIO(raw.decode())))
                if not rows or any(len(r)!=len(rows[0]) for r in rows):raise ValueError('CSV not rectangular')
            else:raw.decode('utf-8')
    if expected:raise ValueError('Archive incomplete')
    return opened


def package():
    files={}
    for p in sorted(HERE.iterdir()):
        if p.is_symlink():raise ValueError('Symlink in payload')
        if not p.is_file() or p.name in EXCLUDE:continue
        if p.suffix not in ('.py','.json','.csv','.npz','.png','.md'):raise ValueError('Unexpected payload type')
        files[PREFIX+p.name]=p.read_bytes()
    required={'README.md','inputs_ledger.md','run.py','inputs.json','results.json','conservation.csv','single_impact_channels.png','single_response.npz','whole_foil.csv','yield_consumers.csv','pair_consumers.csv'}
    if not required<={Path(n).name for n in files}:raise ValueError('Missing required result')
    manifest=dict(base=BASE,unpublished=True,files=[dict(path=n,bytes=len(b),sha256=sha(b)) for n,b in files.items()])
    mn=(json.dumps(manifest,indent=2)+'\n').encode();(HERE/'manifest.json').write_bytes(mn)
    stream=io.BytesIO()
    with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0) as compressed:
        with tarfile.open(fileobj=compressed,mode='w') as archive:
            for n,b in sorted({**files,PREFIX+'manifest.json':mn}.items()):
                info=tarfile.TarInfo(n);info.size=len(b);info.mtime=0;info.mode=0o644
                archive.addfile(info,io.BytesIO(b))
    data=stream.getvalue();opened=verify(data,manifest);(HERE/'capsule.tar.gz').write_bytes(data)
    receipt=dict(filename='capsule.tar.gz',sha256=sha(data),bytes=len(data),files=len(files),all_member_hashes_verified=True,opened=opened)
    (HERE/'delivery.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


def gh(*args,body=None):
    return subprocess.run(['gh','api',*args],input=body,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout


def upload():
    local=json.loads((HERE/'delivery.json').read_text());manifest=json.loads((HERE/'manifest.json').read_text())
    data=(HERE/local['filename']).read_bytes()
    if sha(data)!=local['sha256']:raise ValueError('Local archive changed')
    # Save the release id immediately: an interrupted upload is resumed by id,
    # never by creating a second draft or changing a published release.
    receipt_path=HERE/'remote_delivery.json'
    if receipt_path.exists():
        prior=json.loads(receipt_path.read_text());release=json.loads(gh(f'repos/{REPO}/releases/{prior["draft_id"]}'))
    else:
        body=json.dumps(dict(tag_name=TAG,target_commitish=BASE,name='G15 Pass1: all-track energy and optical/foil consumers',draft=True,prerelease=False,body='Unpublished conditional coarse carrier/whole-foil model, frozen-profile kinematics and exact capsule bytes. See README and inputs ledger. No main-branch push or publication.')).encode()
        release=json.loads(gh(f'repos/{REPO}/releases','--method','POST','--input','-',body=body))
        receipt_path.write_text(json.dumps(dict(draft_id=release['id'],state='created_not_verified',draft=True),indent=2)+'\n')
    if not release['draft'] or release['target_commitish']!=BASE:raise ValueError('Draft/base boundary changed')
    existing=next((a for a in release['assets'] if a['name']==local['filename']),None)
    if existing is None:
        endpoint=release['upload_url'].split('{')[0]+'?name='+local['filename']
        gh(endpoint,'--method','POST','-H','Content-Type: application/gzip','--input',str(HERE/local['filename']))
    remote=json.loads(gh(f'repos/{REPO}/releases/{release["id"]}'))
    if not remote['draft'] or remote['target_commitish']!=BASE:raise ValueError('Readback draft/base mismatch')
    if {a['name'] for a in remote['assets']}!={local['filename']}:raise ValueError('Unexpected draft asset inventory')
    asset=remote['assets'][0]
    raw=gh(f'repos/{REPO}/releases/assets/{asset["id"]}','-H','Accept: application/octet-stream')
    if len(raw)!=local['bytes'] or sha(raw)!=local['sha256']:raise ValueError('Remote archive byte mismatch')
    opened=verify(raw,manifest)
    receipt=dict(draft_id=remote['id'],tag=TAG,draft=True,base=BASE,asset_id=asset['id'],filename=asset['name'],sha256=sha(raw),bytes=len(raw),remote_byte_readback=True,all_member_hashes_verified=True,opened=opened)
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['package','upload']);a=p.parse_args();globals()[a.action]()
