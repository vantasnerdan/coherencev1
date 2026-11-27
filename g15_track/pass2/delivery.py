"""Package only g15_track/pass2 and verify the unpublished draft's exact bytes.
Run from the accepted repository root. No Git mutation or publication.
"""
from pathlib import Path,PurePosixPath
import argparse,gzip,hashlib,io,json,subprocess,tarfile
import numpy as np
HERE=Path(__file__).resolve().parent
BASE='c037c8c479bc08f6ae3ed01916cf8582cafdcf38'
PREFIX='g15_track/pass2/'
REPO='vantasnerdan/coherencev1'
TAG='g15-pass2-causal-optical-phase-20261127'
CORE = {
    'README.md', 'inputs_ledger.md', 'run.py', 'whole.py', 'consumers.py',
    'mechanics.py', 'motion.py', 'linear_join.py', 'observer.py', 'observables.py',
    'checks.py', 'delivery.py', 'results.json', 'source_index.json',
    'source_run.json', 'producer_contract.json', 'controls.json', 'force_checks.json',
    'full_window_checks.json', 'patch_summary.csv', 'phase_map.csv', 'phase_evolution.png',
    'phase_map.png', 'whole_patch_modes.csv', 'critical_phase_drive.csv',
    'phase_transport.csv', 'whole_foil_rates.csv', 'yield_consumers.csv',
    'yield_mixture.csv', 'pair_consumers.csv', 'state_difference.csv',
    'patch_level_rates.csv',
    'velocity_support.csv', 'continued_moving.csv',
}
EVIDENCE = {
    'all_kick_marks.csv', 'phase_fields.npz', 'joint_states.npz',
    'patch_times.csv', 'parent_trials.csv', 'local_modes.csv', 'motif_mode_probes.csv',
}

def sha(b):return hashlib.sha256(b).hexdigest()

def verify(data,manifest):
    expected={r['path']:r for r in manifest['files']};opened=[]
    with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as archive:
        seen=set()
        for member in archive.getmembers():
            name=PurePosixPath(member.name)
            if not member.isfile() or name.is_absolute() or '..' in name.parts or not str(name).startswith(PREFIX) or member.name in seen:
                raise ValueError('Archive member outside additive scope')
            seen.add(member.name);raw=archive.extractfile(member).read()
            if member.name==PREFIX+manifest['manifest']:
                if json.loads(raw)!=manifest:raise ValueError('Embedded manifest mismatch')
                continue
            row=expected.pop(member.name,None)
            if row is None or row['bytes']!=len(raw) or row['sha256']!=sha(raw):raise ValueError('Member byte mismatch')
            if name.suffix=='.npz':
                with np.load(io.BytesIO(raw),allow_pickle=False) as arrays:
                    for k in arrays.files:
                        a=arrays[k]
                        if a.dtype.kind in 'fc' and not np.isfinite(a).all():raise ValueError('Nonfinite saved state')
            elif name.suffix=='.json':json.loads(raw)
            elif name.suffix=='.png':
                from matplotlib.image import imread
                imread(io.BytesIO(raw),format='png')
            else:raw.decode('utf-8')
            opened.append(member.name)
    if expected or PREFIX+manifest['manifest'] not in seen:raise ValueError('Missing archive members')
    return opened

def pack(group,names):
    files={}
    for name in sorted(names):
        p=HERE/name
        if p.is_symlink() or not p.is_file():raise ValueError(f'Expected a regular payload file: {name}')
        b=p.read_bytes()
        if p.suffix in ('.md','.csv','.json'):
            for bad in (b'/home/',b'/tmp/',b'BEGIN PRIVATE KEY',b'github_pat_',b'ghp_'):
                if bad in b:raise ValueError('Private path or credential in public payload')
        files[PREFIX+p.name]=b
    mn='manifest.json' if group=='capsule' else 'evidence_manifest.json'
    manifest=dict(base=BASE,unpublished=True,scope=PREFIX,manifest=mn,
                  files=[dict(path=p,bytes=len(b),sha256=sha(b)) for p,b in files.items()])
    mb=(json.dumps(manifest,indent=2)+'\n').encode();(HERE/mn).write_bytes(mb)
    files[PREFIX+mn]=mb
    stream=io.BytesIO()
    with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0) as compressed:
        with tarfile.open(fileobj=compressed,mode='w') as archive:
            for name,b in files.items():
                info=tarfile.TarInfo(name);info.size=len(b);info.mode=0o644;info.mtime=0
                archive.addfile(info,io.BytesIO(b))
    data=stream.getvalue()
    cap=50_000_000 if group=='capsule' else 2_000_000_000
    if len(data)>=cap:raise ValueError(f'{group} exceeds the byte limit')
    opened=verify(data,manifest);filename=group+'.tar.gz';(HERE/filename).write_bytes(data)
    return dict(filename=filename,manifest=mn,sha256=sha(data),bytes=len(data),
                files=len(opened),all_member_hashes_verified=True)

def package():
    peaks={p.name for p in HERE.glob('peak_*.npz')}
    if len(peaks)!=72:raise ValueError('Expected the 72 accepted source-case arrays')
    receipt=[pack('capsule',CORE),pack('evidence',EVIDENCE|peaks)]
    (HERE/'delivery.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

def gh(*args,body=None):
    return subprocess.run(['gh','api',*args],input=body,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout

def upload():
    local={a['filename']:a for a in json.loads((HERE/'delivery.json').read_text())}
    receipt_path=HERE/'remote_delivery.json'
    if receipt_path.exists():
        old=json.loads(receipt_path.read_text());release=json.loads(gh(f'repos/{REPO}/releases/{old["draft_id"]}'))
    else:
        body=json.dumps(dict(tag_name=TAG,target_commitish=BASE,name='G15 Pass2: causal optical phase and same-foil consumers',draft=True,prerelease=False,
            body='Unpublished conditional signed-force and conserved-carrier refinement. Exact source epoch, finite native pair support and all six fixed-profile settings. No branch push or publication.')).encode()
        release=json.loads(gh(f'repos/{REPO}/releases','--method','POST','--input','-',body=body))
        receipt_path.write_text(json.dumps(dict(draft_id=release['id'],draft=True,state='created_not_verified'),indent=2)+'\n')
    if not release['draft'] or release['target_commitish']!=BASE:raise ValueError('Draft/base boundary mismatch')
    existing={a['name'] for a in release['assets']}
    for filename, row in local.items():
        path=HERE/filename
        with path.open('rb') as stream:
            if hashlib.file_digest(stream,'sha256').hexdigest()!=row['sha256']:
                raise ValueError('Local archive changed')
        if filename not in existing:
            endpoint=release['upload_url'].split('{')[0]+'?name='+filename
            gh(endpoint,'--method','POST','-H','Content-Type: application/gzip','--input',str(path))
    remote=json.loads(gh(f'repos/{REPO}/releases/{release["id"]}'))
    if not remote['draft'] or remote['target_commitish']!=BASE:raise ValueError('Readback draft/base mismatch')
    if {a['name'] for a in remote['assets']}!=set(local):raise ValueError('Unexpected asset inventory')
    verified=[]
    for asset in remote['assets']:
        row=local[asset['name']]
        raw=gh(f'repos/{REPO}/releases/assets/{asset["id"]}','-H','Accept: application/octet-stream')
        if len(raw)!=row['bytes'] or sha(raw)!=row['sha256']:raise ValueError('Remote archive mismatch')
        manifest=json.loads((HERE/row['manifest']).read_text())
        opened=verify(raw,manifest)
        verified.append(dict(asset_id=asset['id'],filename=asset['name'],sha256=row['sha256'],
                              bytes=len(raw),remote_byte_readback=True,
                              all_member_hashes_verified=True,files=len(opened)))
    receipt=dict(draft_id=remote['id'],tag=TAG,draft=True,base=BASE,url=remote['html_url'],assets=verified)
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['package','upload']);args=parser.parse_args();globals()[args.action]()
