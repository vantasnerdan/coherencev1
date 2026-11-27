"""Package and byte-read-back the authorized UNPUBLISHED Pass 3 draft only."""
from pathlib import Path,PurePosixPath
import argparse,gzip,hashlib,io,json,subprocess,tarfile
HERE=Path(__file__).resolve().parent
BASE='2d85230f2580efbee4d7de5e8aece2c9540b3c1c'
PREFIX='g15_track/pass3/'
REPO='vantasnerdan/coherencev1'
TAG='g15-pass3-direct-stiffness-squeeze-20261127'
CORE={'README.md','inputs_ledger.md','run.py','consumers.py','joined.py','receivers.py','four_mode.py','delivery.py',
      'results.json','controls.json','source_index.json','consumer_checks.json','whole_mode_summary.csv','whole_mode_pairs.csv',
      'phase_summary.csv','phase_evolution.csv','pair_instantaneous.csv','pair_trajectory_minima.csv',
      'native_quiet_tails.csv','yield_mixture.csv','acoustic_join.csv','acoustic_transmission.csv',
      'carrier_delta.csv','whole_foil.csv','sibling_consumers.csv','neighbor_consumers.csv',
      'yield_consumers.csv','struck_consumers.csv','governing_plan.json','cage_axial_transfer.csv','strength_and_support.csv'}


def sha(data):return hashlib.sha256(data).hexdigest()


def verify(path,manifest):
    expected={r['path']:r for r in manifest['files']};found=set();manifest_seen=False
    with tarfile.open(path,mode='r:gz') as archive:
        for member in archive.getmembers():
            name=PurePosixPath(member.name)
            if not member.isfile() or name.is_absolute() or '..' in name.parts:raise ValueError('Unsafe archive member')
            if member.name==PREFIX+f'manifest_{manifest["group"]}.json':
                if manifest_seen or json.load(archive.extractfile(member))!=manifest:raise ValueError('Embedded manifest mismatch')
                manifest_seen=True;continue
            if member.name in found or member.name not in expected:raise ValueError('Unexpected or repeated member')
            item=expected[member.name]
            with archive.extractfile(member) as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
            if member.size!=item['bytes'] or digest!=item['sha256']:raise ValueError('Member bytes mismatch')
            found.add(member.name)
    if found!=set(expected) or not manifest_seen:raise ValueError('Incomplete archive')
    return len(found)


def pack(group,names,source_epoch):
    manifest={'base':BASE,'scope':PREFIX,'group':group,'source_epoch_sha256':source_epoch,'files':[]}
    filename=f'coherencev1-g15-pass3-{group}-20261127.tar.gz';target=HERE/filename
    with target.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0,compresslevel=6) as compressed:
            with tarfile.open(fileobj=compressed,mode='w|') as archive:
                for name in sorted(names):
                    path=HERE/name
                    if not path.is_file():raise FileNotFoundError(name)
                    data=path.read_bytes()
                    manifest['files'].append(dict(path=PREFIX+name,bytes=len(data),sha256=sha(data)))
                    info=tarfile.TarInfo(PREFIX+name);info.size=len(data);info.mode=0o644;info.mtime=0
                    archive.addfile(info,io.BytesIO(data))
                encoded=(json.dumps(manifest,indent=2)+'\n').encode()
                info=tarfile.TarInfo(PREFIX+f'manifest_{group}.json');info.size=len(encoded);info.mode=0o644;info.mtime=0
                archive.addfile(info,io.BytesIO(encoded))
    (HERE/f'manifest_{group}.json').write_bytes(encoded)
    size=target.stat().st_size
    if group=='capsule' and size>=50_000_000:raise ValueError('Main capsule exceeds 50 MB')
    if size>=1_900_000_000:raise ValueError('Evidence asset exceeds the bounded transport size')
    with target.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(filename=filename,bytes=size,sha256=digest,manifest=f'manifest_{group}.json',
        files=verify(target,manifest),all_file_hashes_verified=True)


def package():
    scientific={p.name for pattern in ('cov_*.npz','field_pairs_*.npz','band_*.npz') for p in HERE.glob(pattern)}|{'parent_windows.csv'}
    epoch_inputs={'source_index_sha256':sha((HERE/'source_index.json').read_bytes()),'covariances':[]}
    for path in sorted(HERE.glob('cov_*.npz')):
        with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
        epoch_inputs['covariances'].append({'filename':path.name,'sha256':digest})
    source_epoch=sha(json.dumps(epoch_inputs,sort_keys=True,separators=(',',':')).encode())
    receipt=[pack('capsule',CORE,source_epoch)];assigned=set()
    for group,energies in (('evidence-low',(.25,.5)),('evidence-mid',(1.,1.5)),('evidence-high',(2.,2.5))):
        names={name for name in scientific if any(name.startswith(f'{kind}_{ec:g}_') for kind in ('cov','field_pairs','band') for ec in energies)}
        if group=='evidence-low':names.add('parent_windows.csv')
        if names&assigned:raise ValueError('Repeated scientific evidence member')
        assigned|=names;receipt.append(pack(group,names,source_epoch))
    if assigned!=scientific:raise ValueError('Unassigned scientific evidence')
    (HERE/'delivery.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


def gh(*args,body=None):
    return subprocess.run(['gh','api',*args],input=body,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout


def upload():
    local={a['filename']:a for a in json.loads((HERE/'delivery.json').read_text())}
    receipt_path=HERE/'remote_delivery.json'
    if receipt_path.exists():
        old=json.loads(receipt_path.read_text());release=json.loads(gh(f'repos/{REPO}/releases/{old["draft_id"]}'))
    else:
        body=json.dumps(dict(tag_name=TAG,target_commitish=BASE,name='G15 Pass3: direct stiffness and covariance breathing',draft=True,prerelease=False,
            body='Unpublished same-carrier signed-curvature refinement: six fixed-profile settings, three closures and beams, native widths and finite pair support. No branch push or publication.')).encode()
        release=json.loads(gh(f'repos/{REPO}/releases','--method','POST','--input','-',body=body))
        receipt_path.write_text(json.dumps(dict(draft_id=release['id'],draft=True,state='created_not_verified'),indent=2)+'\n')
    if not release['draft'] or release['target_commitish']!=BASE:raise ValueError('Draft/base authority mismatch')
    existing={a['name'] for a in release['assets']}
    if existing-set(local):raise ValueError('Unexpected existing draft assets')
    for filename,row in local.items():
        path=HERE/filename
        with path.open('rb') as stream:
            if hashlib.file_digest(stream,'sha256').hexdigest()!=row['sha256']:raise ValueError('Local archive changed')
        if filename not in existing:
            gh(release['upload_url'].split('{')[0]+'?name='+filename,'--method','POST','-H','Content-Type: application/gzip','--input',str(path))
    remote=json.loads(gh(f'repos/{REPO}/releases/{release["id"]}'))
    if not remote['draft'] or remote['target_commitish']!=BASE:raise ValueError('Readback draft/base mismatch')
    if {a['name'] for a in remote['assets']}!=set(local):raise ValueError('Readback inventory mismatch')
    verified=[]
    for asset in remote['assets']:
        row=local[asset['name']]
        readback=HERE/'remote_readback.tar.gz'
        try:
            with readback.open('wb') as stream:
                subprocess.run(['gh','api',f'repos/{REPO}/releases/assets/{asset["id"]}','-H','Accept: application/octet-stream'],
                    stdout=stream,stderr=subprocess.PIPE,check=True)
            size=readback.stat().st_size
            with readback.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
            if size!=row['bytes'] or digest!=row['sha256']:raise ValueError('Remote archive bytes mismatch')
            manifest=json.loads((HERE/row['manifest']).read_text())
            verified.append(dict(asset_id=asset['id'],filename=asset['name'],sha256=row['sha256'],bytes=size,
                remote_byte_readback=True,all_member_hashes_verified=True,files=verify(readback,manifest)))
        finally:
            readback.unlink(missing_ok=True)
    receipt=dict(draft_id=remote['id'],tag=TAG,draft=True,base=BASE,url=remote['html_url'],assets=verified)
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['package','upload'])
    args=parser.parse_args();globals()[args.action]()
