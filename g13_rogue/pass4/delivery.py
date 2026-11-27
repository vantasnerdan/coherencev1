"""Whitelist/package and authenticate the authorized unpublished draft bytes.
No Git commit, tag push, publication, or merge. Uses the operator's existing gh.
"""
from pathlib import Path, PurePosixPath
import argparse
import gzip
import hashlib
import io
import json
import subprocess
import tarfile

HERE=Path(__file__).resolve().parent
PREFIX='g13_rogue/pass4/'
REPO='vantasnerdan/coherencev1'
BASE='650a2ebd57c35f371ecd23c4db975ebd6e3cc55a'
TAG='g13-rogue-pass4-transport-20261127'


def digest(data):return hashlib.sha256(data).hexdigest()


def package():
    groups={'capsule':{},'evidence':{}}
    excluded={'manifest.json','evidence_manifest.json','delivery.json','remote_delivery.json'}
    for p in sorted(HERE.iterdir()):
        if p.is_symlink():raise ValueError('Payload contains a symlink')
        if not p.is_file() or p.name in excluded or p.suffix not in ('.py','.md','.csv','.json','.npz'):continue
        if p.name.endswith('_r0.npz') or p.name.endswith('_r1.npz') or p.name.endswith('_r2.npz') or p.name.endswith('_r3.npz') or p.name in ('joint_distributions.npz','slow_background.npz'):
            group='evidence'
        else:group='capsule'
        data=p.read_bytes()
        if p.suffix!='.npz' and (b'/' + b'home/') in data:raise ValueError('Host-local path in portable payload: '+p.name)
        groups[group][PREFIX+p.name]=data
    receipts={}
    for group,files in groups.items():
        manifest={'base':BASE,'unpublished':True,'files':[dict(path=n,bytes=len(b),sha256=digest(b)) for n,b in files.items()]}
        mn='manifest.json' if group=='capsule' else 'evidence_manifest.json'
        mb=(json.dumps(manifest,indent=2)+'\n').encode();(HERE/mn).write_bytes(mb)
        files[PREFIX+mn]=mb
        out=io.BytesIO()
        with gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as gz:
            with tarfile.open(mode='w|',fileobj=gz) as tf:
                for name,data in files.items():
                    ti=tarfile.TarInfo(name);ti.size=len(data);ti.mode=0o644;ti.mtime=0
                    tf.addfile(ti,io.BytesIO(data))
        data=out.getvalue();name=group+'.tar.gz';(HERE/name).write_bytes(data)
        verify(data,manifest,mn)
        if group=='capsule' and len(data)>50_000_000:raise ValueError('Main capsule exceeds authorized 50 MB')
        receipts[group]=dict(path=PREFIX+name,bytes=len(data),sha256=digest(data),members=len(files),whitelist_and_member_hashes_verified=True)
    (HERE/'delivery.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print(json.dumps(receipts))


def verify(data,manifest,mn):
    expected={r['path']:r for r in manifest['files']};seen=False
    with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as tf:
        for member in tf:
            p=PurePosixPath(member.name)
            if not member.isfile() or p.is_absolute() or '..' in p.parts or not member.name.startswith(PREFIX):raise ValueError('Unsafe archive member')
            b=tf.extractfile(member).read()
            if member.name==PREFIX+mn:
                if seen or json.loads(b)!=manifest:raise ValueError('Manifest identity')
                seen=True;continue
            r=expected.pop(member.name)
            if member.size!=r['bytes'] or digest(b)!=r['sha256']:raise ValueError('Member identity '+member.name)
    if expected or not seen:raise ValueError('Incomplete archive')


def gh(*args,body=None):
    result=subprocess.run(['gh','api',*args],input=body,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
    return result.stdout


def upload():
    local=json.loads((HERE/'delivery.json').read_text())
    body=json.dumps(dict(tag_name=TAG,target_commitish=BASE,name='G13 Pass 4 — conditional whole-foil transport',draft=True,prerelease=False,body='Unpublished scoped Pass-4 capsule and separate evidence. No publication or merge. Source/limitations and root-relative replay are in the capsule.')).encode()
    release=json.loads(gh(f'repos/{REPO}/releases','--method','POST','--input','-',body=body))
    if not release['draft']:raise ValueError('Draft boundary violated')
    (HERE/'remote_delivery.json').write_text(json.dumps(dict(draft_id=release['id'],tag=TAG,draft=True,state='created_not_yet_verified'),indent=2)+'\n')
    assets=[]
    for group,row in local.items():
        path=HERE/Path(row['path']).name
        if digest(path.read_bytes())!=row['sha256']:raise ValueError('Changed archive')
        endpoint=release['upload_url'].split('{')[0]+'?name='+path.name
        gh(endpoint,'--method','POST','-H','Content-Type: application/gzip','--input',str(path))
    remote=json.loads(gh(f'repos/{REPO}/releases/{release["id"]}'))
    if not remote['draft']:raise ValueError('Release ceased being draft')
    if remote['target_commitish']!=BASE:raise ValueError('Draft target base changed')
    if {x['name'] for x in remote['assets']}!={Path(x['path']).name for x in local.values()}:raise ValueError('Draft asset inventory differs from whitelist')
    for asset in remote['assets']:
        group=asset['name'].removesuffix('.tar.gz')
        row=local[group]
        raw=gh(f'repos/{REPO}/releases/assets/{asset["id"]}','-H','Accept: application/octet-stream')
        if len(raw)!=row['bytes'] or digest(raw)!=row['sha256']:raise ValueError('Remote archive identity')
        mn='manifest.json' if group=='capsule' else 'evidence_manifest.json'
        verify(raw,json.loads((HERE/mn).read_text()),mn)
        assets.append(dict(group=group,id=asset['id'],name=asset['name'],bytes=len(raw),sha256=digest(raw),remote_byte_readback=True,all_members_read_back=True))
    receipt=dict(draft_id=release['id'],tag=TAG,draft=True,assets=assets)
    (HERE/'remote_delivery.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['package','upload']);a=p.parse_args();globals()[a.action]()
