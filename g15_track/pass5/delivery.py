"""Selected SHA-bound inputs, lossless thermal codec and private draft readback."""
from pathlib import Path,PurePosixPath
import csv,gzip,hashlib,io,json,lzma,os,subprocess,sys,tarfile,tempfile,zipfile
import numpy as np
from common import HERE,BASE,ENERGIES,verify_selected,digest,emit,table

PREFIX='g15_track/pass5/'
REPO='vantasnerdan/coherencev1'
TAG='g15-pass5-whole-foil-coherence-20261128'


def prepare(path4,out=HERE):
    from medium import PREDECESSOR_FILES
    names=set(PREDECESSOR_FILES)-{'manifest_core.json','manifest_evidence.json'}
    names.update(('yielding.py','yield_mixture.csv','frozen/bca.py','frozen/native_transfer.npz','frozen/native_optical_basis.npz','decision.json','frozen/prior_inflight_yields.csv'))
    for ec in ENERGIES:
        for closure in ('fast_Cu','nominal_retained','long_retained'):
            for sp in (0,1):
                names.add(f'source_E{ec:g}_{closure}_sp{sp}.npz')
                for allocation in ('onsite','coupling'):names.add(f'native_E{ec:g}_{closure}_{allocation}_sp{sp}_m16.npz')
    source=Path(path4);out=Path(out);rows=verify_selected(source,names)
    for row in rows:
        target=out/'inherited'/row['path'];target.parent.mkdir(parents=True,exist_ok=True)
        data=(source/row['path']).read_bytes()
        if target.exists() and target.read_bytes()!=data:raise ValueError('Selected input already differs: '+row['path'])
        if not target.exists():target.write_bytes(data)
    manifests=[]
    for name in ('manifest_core.json','manifest_evidence.json'):
        target=out/'inherited'/name;target.write_bytes((source/name).read_bytes());manifests.append(dict(path=name,sha256=digest(target)))
    emit(out/'selected_inputs.json',dict(base=BASE,source='exact accepted Pass4 selected members',files=rows,accepted_manifests=manifests))
    from low import build
    build(out,source)
    return rows


def compact_medium(out=HERE):
    out=Path(out);path=out/'medium.npz';rows=[]
    with np.load(path) as source:
        if any(k.endswith('_lzma') for k in source.files):return json.loads((out/'medium_codec.json').read_text())
        target=out/'medium.compact.npz'
        with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
            def save(name,array,compressed=True):
                buffer=io.BytesIO();np.save(buffer,array,allow_pickle=False)
                archive.writestr(name+'.npy',buffer.getvalue(),compress_type=zipfile.ZIP_DEFLATED if compressed else zipfile.ZIP_STORED)
            for name in source.files:
                array=source[name]
                if name.startswith('T_') and name.endswith('_depth'):
                    if array.dtype!=np.dtype('float32'):raise ValueError('Thermal codec requires original binary32 history')
                    delta=np.diff(array.view(np.uint32),axis=-1,prepend=np.uint32(0))
                    shuffled=delta.reshape(-1).view(np.uint8).reshape(-1,4).T.copy().tobytes()
                    payload=lzma.compress(shuffled,preset=3)
                    raw=np.frombuffer(lzma.decompress(payload),dtype=np.uint8).reshape(4,-1).T.copy().reshape(-1).view(np.uint32).reshape(array.shape)
                    decoded=np.cumsum(raw,axis=-1,dtype=np.uint32).view(np.float32)
                    if decoded.tobytes()!=array.tobytes():raise ValueError('Thermal codec changes physical bits')
                    save(name+'_lzma',np.frombuffer(payload,dtype=np.uint8),False);save(name+'_shape',np.array(array.shape,dtype=np.int64))
                    rows.append(dict(array=name,shape=list(array.shape),decoded_bytes=array.nbytes,sha256=hashlib.sha256(array.tobytes()).hexdigest(),packed_bytes=len(payload),bit_identical=True))
                else:save(name,array)
    receipt=dict(codec='uint32-depth-delta-byte-shuffle-lzma-v1',original_bytes=path.stat().st_size,compact_bytes=target.stat().st_size,arrays=rows,physical_recomputation=False)
    target.replace(path);emit(out/'medium_codec.json',receipt);return receipt


def verify(archive,manifest,extract=None):
    expected={r['path']:r for r in manifest['files']};seen=set()
    with tarfile.open(archive,'r:gz') as stream:
        for member in stream.getmembers():
            if not member.isfile() or member.name not in expected:raise ValueError('Unexpected archive member')
            rel=PurePosixPath(member.name)
            if rel.is_absolute() or '..' in rel.parts or not member.name.startswith(PREFIX):raise ValueError('Unsafe archive path')
            data=stream.extractfile(member).read();row=expected[member.name]
            if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:raise ValueError('Archive member changed')
            if member.name in seen:raise ValueError('Duplicate archive member')
            seen.add(member.name)
            if extract:
                target=Path(extract)/member.name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    if seen!=set(expected):raise ValueError('Incomplete selected archive')
    return len(seen)


def package(out=HERE):
    out=Path(out);groups={'core':[],'evidence':[]}
    result=json.loads((out/'results.json').read_text())
    if (result['receiver_rows'],result['yield_rows'],result['meeting_cases'],result['matched_frozen_launch_control_cases'])!=(648,648,162,54):raise ValueError('Incomplete final case grid')
    if result['consumer_identity']!={p.name:digest(p) for p in sorted(out.glob('*.py'))}:raise ValueError('Final report uses a different consumer source')
    exclude={'delivery.json','remote_delivery.json','cold_bundle_consume.json','cold_consume.json'}
    evidence={'medium.npz','low_modes.npz','meetings.npz','acoustic_probe.npz','nativeGX_meetings.npz','nativeGX_acoustic_probe.npz','native_Gamma_X_receiver.npz'}
    control_rows=[r for r in csv.DictReader((out/'nativeGX_yield_outcomes.csv').open()) if r['allocation']=='coupling_frozen_launch']
    table(out/'frozen_launch_outcomes.csv',control_rows)
    for path in out.rglob('*'):
        relative=path.relative_to(out)
        if not path.is_file() or any(part.startswith('.') or part=='__pycache__' for part in relative.parts):continue
        if path.name in exclude or path.suffix=='.gz' or path.name.startswith('manifest_') and path.parent==out:continue
        if path.parent==out:
            if path.name.startswith('first_'):continue
            if path.name.startswith('nativeGX_E'):continue
            if path.name.startswith('nativeGX_') and path.suffix=='.csv' or path.name=='nativeGX_consumer_checks.json':continue
            if path.name.startswith('E') and not path.name.endswith('_receiver.npz'):continue
            if path.suffix=='.log' or path.name.endswith('.compact.npz'):continue
            if path.suffix=='.npz' and path.name not in evidence and not path.name.startswith('spectrum_') and not path.name.endswith('_receiver.npz'):continue
        group='evidence' if path.suffix=='.npz' else 'core'
        groups[group].append(path)
    receipts=[]
    for group,paths in groups.items():
        files=[dict(path=PREFIX+p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(paths)]
        manifest=dict(base=BASE,scope=PREFIX,group=group,files=files)
        emit(out/f'manifest_{group}.json',manifest)
        filename=f'coherencev1-g15-pass5-{group}-20261128.tar.gz';target=out/filename
        with target.open('wb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w') as archive:
            for row,path in zip(files,sorted(paths)):
                info=tarfile.TarInfo(row['path']);info.size=row['bytes'];info.mtime=0;info.mode=0o644;archive.addfile(info,io.BytesIO(path.read_bytes()))
        cap=50_000_000 if group=='core' else 300_000_000
        if target.stat().st_size>=cap:raise ValueError('Approved archive cap exceeded: '+group)
        receipts.append(dict(filename=filename,bytes=target.stat().st_size,sha256=digest(target),manifest=f'manifest_{group}.json',files=verify(target,manifest)))
    emit(out/'delivery.json',receipts);return receipts


def cold(python,out=HERE):
    out=Path(out);receipts=json.loads((out/'delivery.json').read_text())
    with tempfile.TemporaryDirectory(prefix='g15-pass5-cold-') as root:
        count=0
        for row in receipts:
            path=out/row['filename']
            if digest(path)!=row['sha256']:raise ValueError('Archive changed before consumer')
            count+=verify(path,json.loads((out/row['manifest']).read_text()),root)
        command=[python,str(Path(root)/PREFIX/'run.py'),'--consume']
        result=subprocess.run(command,cwd=root,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,check=True)
        compared=[]
        for name in ('yield_outcomes.csv','consumer_checks.json','probe_projection_holes.csv','probe_projection_samples.csv','frozen_launch_controls.csv'):
            original=(out/name).read_bytes();actual=(Path(root)/PREFIX/name).read_bytes()
            if original!=actual:raise ValueError('Cold consumer changed result: '+name)
            compared.append(dict(path=name,sha256=hashlib.sha256(actual).hexdigest(),byte_identical=True))
        receipt=dict(success=True,archives=receipts,regular_files_verified=count,consumer_command=['python','g15_track/pass5/run.py','--consume'],
            consumer=json.loads((Path(root)/PREFIX/'cold_consume.json').read_text()),outputs=compared,stdout=result.stdout.strip(),selected_archives_only=True,
            member_manifests={row['manifest']:json.loads((out/row['manifest']).read_text()) for row in receipts})
    emit(out/'cold_bundle_consume.json',receipt);return receipt


def gh(*args,body=None):
    return subprocess.run(['gh','api',*args],input=body,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout


def upload(out=HERE):
    out=Path(out);local={r['filename']:r for r in json.loads((out/'delivery.json').read_text())}
    repository=json.loads(gh(f'repos/{REPO}'));actor=json.loads(gh('user'))['login']
    if repository['full_name'].lower()!=REPO.lower():raise ValueError('Wrong authenticated repository')
    proof=out/'cold_bundle_consume.json';cold_result=json.loads(proof.read_text())
    if not cold_result['success'] or cold_result['archives']!=list(local.values()):raise ValueError('Cold proof does not bind selected bytes')
    local[proof.name]=dict(filename=proof.name,bytes=proof.stat().st_size,sha256=digest(proof),manifest=None)
    receipt_path=out/'remote_delivery.json'
    if receipt_path.exists():release=json.loads(gh(f'repos/{REPO}/releases/{json.loads(receipt_path.read_text())["draft_id"]}'))
    else:
        body=json.dumps(dict(tag_name=TAG,target_commitish=BASE,name='G15 Pass 5: conditional whole-foil signed receiver and velocity consumer',draft=True,prerelease=False,
            body='Private unpublished scientific draft on the exact accepted Pass4 base. Selected-input replay, signed receiver laws, depth landscape, conditioned independent-parent optical meetings and six-energy replacement velocity consumer. Conditional model limitations remain explicit. No Git branch or public release.')).encode()
        release=json.loads(gh(f'repos/{REPO}/releases','--method','POST','--input','-',body=body))
        emit(receipt_path,dict(draft_id=release['id'],draft=True,state='created_not_verified'))
    if not release['draft'] or release['target_commitish']!=BASE:raise ValueError('Private draft/base boundary violated')
    existing={a['name'] for a in release['assets']}
    if existing-set(local):raise ValueError('Unexpected existing asset; preserve owner state')
    for name,row in local.items():
        if digest(out/name)!=row['sha256']:raise ValueError('Selected local bytes changed')
        if name not in existing:gh(release['upload_url'].split('{')[0]+'?name='+name,'--method','POST','-H','Content-Type: '+('application/json' if row['manifest'] is None else 'application/gzip'),'--input',str(out/name))
    remote=json.loads(gh(f'repos/{REPO}/releases/{release["id"]}'))
    if not remote['draft'] or remote['target_commitish']!=BASE or {a['name'] for a in remote['assets']}!=set(local):raise ValueError('Remote draft inventory/base mismatch')
    verified=[]
    for asset in remote['assets']:
        row=local[asset['name']]
        with tempfile.NamedTemporaryFile() as tmp:
            subprocess.run(['gh','api',f'repos/{REPO}/releases/assets/{asset["id"]}','-H','Accept: application/octet-stream'],stdout=tmp,stderr=subprocess.PIPE,check=True)
            tmp.flush()
            if Path(tmp.name).stat().st_size!=row['bytes'] or digest(tmp.name)!=row['sha256']:raise ValueError('Remote download/hash mismatch')
            count=1 if row['manifest'] is None else verify(tmp.name,json.loads((out/row['manifest']).read_text()))
        verified.append(dict(asset_id=asset['id'],filename=asset['name'],bytes=row['bytes'],sha256=row['sha256'],remote_download_hash_verified=True,regular_member_hashes_verified=count))
    receipt=dict(draft_id=remote['id'],draft=True,base=BASE,url=remote['html_url'],published=False,
        repository=repository['full_name'],repository_private=repository['private'],authenticated_actor=actor,assets=verified)
    emit(receipt_path,receipt);return receipt
