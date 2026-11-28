"""Private draft packaging, cold consumption and authenticated byte readback only."""
from pathlib import Path,PurePosixPath
import argparse,gzip,hashlib,io,json,os,subprocess,tarfile,tempfile
HERE=Path(__file__).resolve().parent;PREFIX='g15_track/pass4/'
BASE='fd76b0cc3a26ac5e77f9da0577d46c9c583b7eb1';REPO='vantasnerdan/coherencev1';TAG='g15-pass4-whole-foil-waves-20261127'
EXCLUDE={'delivery.json','remote_delivery.json','cold_bundle_consume.json','cold_consume.json','native_first_closure.json'}


def sha(data):return hashlib.sha256(data).hexdigest()


def verify(path,manifest,extract=None):
    expected={r['path']:r for r in manifest['files']};seen=set();embedded=False
    with tarfile.open(path,'r:gz') as archive:
        for member in archive.getmembers():
            name=PurePosixPath(member.name)
            if not member.isfile() or name.is_absolute() or '..' in name.parts:raise ValueError('Unsafe archive member')
            data=archive.extractfile(member).read()
            if member.name==PREFIX+f'manifest_{manifest["group"]}.json':
                if embedded or json.loads(data)!=manifest:raise ValueError('Embedded manifest differs')
                embedded=True;continue
            row=expected.get(member.name)
            if row is None or member.name in seen:raise ValueError('Unexpected archive member')
            if len(data)!=row['bytes'] or sha(data)!=row['sha256']:raise ValueError('Member byte mismatch')
            seen.add(member.name)
            if extract is not None:
                target=Path(extract)/member.name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    if not embedded or seen!=set(expected):raise ValueError('Incomplete archive inventory')
    return len(seen)


def package():
    groups={'core':[],'evidence':[]}
    import sources as source
    obsolete={f'{prefix}_{beam}_{arm["name"]}{suffix}.npz'
        for prefix in ('foil','phase','acoustic_phase') for beam in source.BEAMS
        for arm in source.PHYSICAL_ARMS for suffix in ('','_temperature_off','_strain_off','_frozen')}
    obsolete|={f'paths_{arm["name"]}.npz' for arm in source.PHYSICAL_ARMS}
    phase_kernels={f'acoustic_phase_E2.5_{beam}_{arm["name"]}.npz' for beam in source.BEAMS
        for arm in source.ARMS if arm['closure']=='nominal_retained'}
    for path in HERE.rglob('*'):
        if not path.is_file():continue
        rel=path.relative_to(HERE).as_posix()
        if '__pycache__' in path.parts or path.suffix in ('.gz','.log') or path.name.startswith(('manifest_','.case_')) or path.name in EXCLUDE or path.name in obsolete:continue
        if path.name.startswith('acoustic_phase_E') and path.name not in phase_kernels:continue
        if path.suffix not in ('.py','.md','.csv','.json','.npz','.png'):raise ValueError('Unclassified pass4 file '+rel)
        group='evidence' if path.suffix=='.npz' or rel.startswith('frozen/') or path.name=='conditioned_optical.csv' else 'core'
        groups[group].append(rel)
    receipts=[]
    for group,names in groups.items():
        manifest=dict(base=BASE,scope=PREFIX,group=group,files=[]);filename=f'coherencev1-g15-pass4-{group}-20261127.tar.gz';target=HERE/filename
        raw=sum((HERE/n).stat().st_size for n in names);cap=50_000_000 if group=='core' else 300_000_000
        if raw>=cap:raise ValueError('Uncompressed '+group+' exceeds requested cap')
        with target.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0,compresslevel=6) as compressed:
                with tarfile.open(fileobj=compressed,mode='w|') as archive:
                    for name in sorted(names):
                        data=(HERE/name).read_bytes();manifest['files'].append(dict(path=PREFIX+name,bytes=len(data),sha256=sha(data)))
                        info=tarfile.TarInfo(PREFIX+name);info.size=len(data);info.mode=0o644;info.mtime=0;archive.addfile(info,io.BytesIO(data))
                    encoded=(json.dumps(manifest,indent=2)+'\n').encode();info=tarfile.TarInfo(PREFIX+f'manifest_{group}.json');info.size=len(encoded);info.mode=0o644;info.mtime=0;archive.addfile(info,io.BytesIO(encoded))
        (HERE/f'manifest_{group}.json').write_bytes(encoded)
        if target.stat().st_size>=cap:raise ValueError('Compressed '+group+' exceeds cap')
        receipts.append(dict(group=group,filename=filename,bytes=target.stat().st_size,uncompressed_member_bytes=raw,
            sha256=sha(target.read_bytes()),manifest=f'manifest_{group}.json',files=verify(target,manifest)))
    (HERE/'delivery.json').write_text(json.dumps(receipts,indent=2)+'\n');print(json.dumps(receipts),flush=True)


def cold(python):
    receipts=json.loads((HERE/'delivery.json').read_text())
    with tempfile.TemporaryDirectory(prefix='g15-wave-cold-') as directory:
        files=0
        for row in receipts:
            target=HERE/row['filename']
            if sha(target.read_bytes())!=row['sha256']:raise ValueError('Archive changed before cold extraction')
            files+=verify(target,json.loads((HERE/row['manifest']).read_text()),directory)
        command=[python,str(Path(directory)/PREFIX/'run.py'),'--consume'];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        result=subprocess.run(command,cwd=directory,env=env,capture_output=True,text=True,check=True)
        saved=json.loads((Path(directory)/PREFIX/'cold_consume.json').read_text())
        # Compare final consumer outputs from the fresh extraction, not only exit status.
        compared=[]
        for name in ('yield_consumers.csv','yield_mixture.csv','optical_phase.csv','consumer_checks.json'):
            original=(HERE/name).read_bytes();replayed=(Path(directory)/PREFIX/name).read_bytes()
            if sha(original)!=sha(replayed):raise ValueError('Cold consumer output differs: '+name)
            compared.append(dict(path=PREFIX+name,bytes=len(replayed),sha256=sha(replayed),byte_identical=True))
        receipt=dict(success=True,archives=receipts,cold_regular_files_verified=files,stdout=result.stdout.strip(),consumer=saved,outputs=compared,
            evidence_bytes=sum(x['bytes'] for x in receipts if x['group']=='evidence'),no_original_worktree_or_modal_archive_needed=True)
    (HERE/'cold_bundle_consume.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


def gh(*args,body=None):
    return subprocess.run(['gh','api',*args],input=body,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout


def upload():
    local={x['filename']:x for x in json.loads((HERE/'delivery.json').read_text())};receipt_path=HERE/'remote_delivery.json'
    proof=HERE/'cold_bundle_consume.json';cold_result=json.loads(proof.read_text())
    if not cold_result['success'] or cold_result['archives']!=list(local.values()):raise ValueError('Cold evidence does not bind these archives')
    local[proof.name]=dict(filename=proof.name,bytes=proof.stat().st_size,sha256=sha(proof.read_bytes()),manifest=None)
    if receipt_path.exists():
        old=json.loads(receipt_path.read_text());release=json.loads(gh(f'repos/{REPO}/releases/{old["draft_id"]}'))
    else:
        body=json.dumps(dict(tag_name=TAG,target_commitish=BASE,name='G15 Pass4: whole-foil coupled traveling fields',draft=True,prerelease=False,
            body='UNPUBLISHED replacement wave-model draft. Same accepted foil, three beams, causal carrier ports, reciprocal strain/optical/electronic work, field histories, controls and six-setting final velocity consumer. No source branch push or publication.')).encode()
        release=json.loads(gh(f'repos/{REPO}/releases','--method','POST','--input','-',body=body))
        receipt_path.write_text(json.dumps(dict(draft_id=release['id'],draft=True,state='created_not_verified'),indent=2)+'\n')
    if not release['draft'] or release['target_commitish']!=BASE:raise ValueError('Draft/base authority mismatch')
    existing={x['name'] for x in release['assets']}
    if existing-set(local):raise ValueError('Unexpected existing asset; preserve owner state')
    for name,row in local.items():
        path=HERE/name
        if sha(path.read_bytes())!=row['sha256']:raise ValueError('Local bytes changed')
        media='application/json' if row['manifest'] is None else 'application/gzip'
        if name not in existing:gh(release['upload_url'].split('{')[0]+'?name='+name,'--method','POST','-H','Content-Type: '+media,'--input',str(path))
    remote=json.loads(gh(f'repos/{REPO}/releases/{release["id"]}'))
    if not remote['draft'] or remote['target_commitish']!=BASE or {x['name'] for x in remote['assets']}!=set(local):raise ValueError('Remote draft inventory/base mismatch')
    verified=[]
    for asset in remote['assets']:
        row=local[asset['name']]
        with tempfile.NamedTemporaryFile(suffix='.tar.gz') as tmp:
            subprocess.run(['gh','api',f'repos/{REPO}/releases/assets/{asset["id"]}','-H','Accept: application/octet-stream'],stdout=tmp,stderr=subprocess.PIPE,check=True)
            tmp.flush();data=Path(tmp.name).read_bytes()
            if len(data)!=row['bytes'] or sha(data)!=row['sha256']:raise ValueError('Remote archive byte mismatch')
            if row['manifest'] is None:
                if json.loads(data)!=cold_result:raise ValueError('Remote cold evidence differs')
                count=1
            else:count=verify(tmp.name,json.loads((HERE/row['manifest']).read_text()))
        verified.append(dict(asset_id=asset['id'],filename=asset['name'],bytes=row['bytes'],sha256=row['sha256'],regular_member_hashes_verified=count,remote_byte_readback=True))
    receipt=dict(draft_id=remote['id'],tag=TAG,base=BASE,draft=True,url=remote['html_url'],published=False,assets=verified)
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['package','cold','upload']);parser.add_argument('--python',default='python3');args=parser.parse_args()
    if args.action=='cold':cold(args.python)
    else:globals()[args.action]()
