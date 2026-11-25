"""One unpublished root-relative capsule with frozen input bytes, no publication."""
from pathlib import Path,PurePosixPath
import gzip
import hashlib
import io
import json
import subprocess
import tarfile
import tempfile
import numpy as np
from continuation import HERE,RESULTS,TRACES,load,save

ROOT = HERE.parents[2]
PREFIX = HERE.relative_to(ROOT)
FOCUS = ('README.md','method.md','keep-drop.md','patterns.md','rate_mapping.md',
    'engine.py','patterns.py','run.py','rate_mapping.py','analyze.py','verify.py','package.py',
    'results.json','comparison.csv','trajectories.npz','launches.npz','consumer_proof.json',
    'focusing.png','input_manifest.json','large_geometry.py','continuation.py','exterior.py',
    'continuation_verify.py','continuation_package.py','continuation_results.json',
    'continuation_traces.npz','continuation_states.npz')


def compact(data,traces):
    """All-step scoring stays in candidate rows; retain selected all-step traces."""
    if data.get('archive_trace_policy'):
        return
    for key,case in data['cases'].items():
        for row in case['candidates']:
            if row.get('trace_compacted'):
                continue
            row['peak_integration_step'] = int(round(row['peak_time_fs']/case['dt_fs']))
            name = key+'__'+row['candidate']
            trace = traces[name]
            assert trace[row['peak_index'],0]==row['peak_time_fs']
            assert trace[row['peak_index'],1]==row['peak_relative_kinetic_eV']
            if row.get('selected'):
                row['saved_trace_policy'] = 'all integration steps'
                continue
            stride = max(1,int(round(2/case['dt_fs'])))
            # Keep the peak, extrema and endpoint in addition to the 2-fs grid.
            events = [int(round(e['time_fs']/case['dt_fs'])) for e in
                      row.get('all_step_original_support_transitions',[])]
            chosen = np.unique(np.r_[np.arange(0,len(trace),stride),row['peak_index'],
                np.argmax(trace[:,4]),np.argmin(trace[:,5]),len(trace)-1,events]).astype(int)
            traces[name] = trace[chosen]
            row['peak_index'] = int(np.flatnonzero(chosen==row['peak_index'])[0])
            row['trace_compacted'] = True
            row['saved_trace_policy'] = '2-fs grid plus actual selected-in-candidate peak, radius/bond extrema and endpoint; summary used all steps'
            assert traces[name][row['peak_index'],0]==row['peak_time_fs']
    data['archive_trace_policy'] = 'selected traces every step; raw nonselected candidate summaries every step, compact traces as labelled; exact deterministic source reconstructs full traces'
    save(data,traces)


def package():
    data,traces = load()
    if 'verification' not in data:
        raise RuntimeError('settled changed-consumer verification must precede packaging')
    compact(data,traces)
    inputs = json.loads((HERE/'input_manifest.json').read_text())['inputs']
    payload = {}
    for name in FOCUS:
        payload[str(PREFIX/name)] = (HERE/name).read_bytes()
    for item in inputs:
        content = (ROOT/item['path']).read_bytes()
        assert len(content)==item['bytes'] and hashlib.sha256(content).hexdigest()==item['sha256']
        payload[item['path']] = content
    records = [dict(path=p,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()) for p,b in payload.items()]
    manifest = dict(schema='pdd-ordered-continuation-manifest-v1',publication='unpublished local only',
        base='accepted focusing main 01d56ff2318bd2fbb334a22e55a4f5fe92e6c9f6',
        plan='coherencev1-g11-focusing-continuation-20261125',
        files=records,manifest_self_excluded=True,
        portable_input_layout='Exact prior bytes included at original imported root-relative paths; no live model input was rebound',
        archive='continuation.tar.gz',archive_hash='reported separately to avoid circular archive identity')
    manifest_bytes = (json.dumps(manifest,indent=2)+'\n').encode()
    (HERE/'continuation_manifest.json').write_bytes(manifest_bytes)
    payload[str(PREFIX/'continuation_manifest.json')] = manifest_bytes
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw,mode='w') as tf:
        for name,content in payload.items():
            item = tarfile.TarInfo(name); item.size=len(content); item.mode=0o644; item.mtime=0
            tf.addfile(item,io.BytesIO(content))
    compressed = gzip.compress(raw.getvalue(),compresslevel=9,mtime=0)
    destination = HERE/'continuation.tar.gz'
    destination.write_bytes(compressed)
    with tempfile.TemporaryDirectory(prefix='pdd-continuation-receiver-') as temp:
        root = Path(temp)
        with tarfile.open(fileobj=io.BytesIO(compressed),mode='r:gz') as tf:
            for member in tf.getmembers():
                path = PurePosixPath(member.name)
                if path.is_absolute() or '..' in path.parts or not member.isfile() or member.name not in payload:
                    raise ValueError('unsafe or unlisted capsule member')
                target = root/member.name
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes(tf.extractfile(member).read())
        for name,content in payload.items():
            assert (root/name).read_bytes()==content
        subprocess.run(['python','-c','import continuation, exterior, large_geometry; print("portable archived consumer imports OK")'],
            cwd=root/PREFIX,check=True)
    print(json.dumps(dict(path=str(destination),files=len(payload),bytes=len(compressed),
        sha256=hashlib.sha256(compressed).hexdigest(),verification='safe fresh extraction, every byte equal, cold archived consumer imports'),indent=2))


if __name__=='__main__':
    package()
