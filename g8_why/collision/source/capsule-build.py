"""Seal the approved G8 collision page using the existing ORIGINAL/EXPORTED builder.

Inputs are already-computed artifacts. No Native calculation, physics fitting,
network request, source publication or release publication is performed here.
"""
import argparse
import gzip
import io
import json
from pathlib import Path
import re
import sys
import tarfile
import time


def build(spec, output, archive, receipt):
    started = time.process_time()
    sys.path.insert(0,spec['existing_builder_directory'])
    from build_capsule import Builder, digest, json_bytes, PRIVATE_PATH
    if output.exists() or archive.exists():
        raise ValueError('Do not overwrite a delivered capsule')
    output.mkdir(parents=True)
    builder = Builder(output)
    for item in spec['files']:
        source = Path(item['source']['path'])
        actual = digest(source.read_bytes())
        expected = {key:item['source'][key] for key in ('bytes','sha256')}
        if actual != expected:
            raise ValueError('Original source changed:'+source.name)
        builder.add(source,item['destination'],exact=item.get('exact',False))
    index = {}
    for name,item in builder.files.items():
        index.setdefault(item['ORIGINAL']['sha256'],[]).append(name)
    manifest = dict(schema='g8-collision-capsule-v1',authority_entry='entry_3259',
        publishing=False,source_commit=spec['source_commit'],files=builder.files,
        original_sha256_index=index,closed_omissions=spec['closed_omissions'],
        budget=spec['budget'],reader='python collision/reader.py .',
        namespace='collision.md and collision/**; other mechanism pages are outside this capsule',
        native_invocations=0,new_grid_runs=0,
        provenance='ORIGINAL byte counts/SHA256 bind original files. EXPORTED differences are private path values only; numeric arrays/scalars and complete findings are unchanged. Original private paths cannot be recovered. Manifest is externally bound by archive SHA256, not a signature.',
        excluded='No PDFs, GPWs, raw GPAW user headers, provider/session payloads or launch/watch administrative logs. Full nuclear eigenvector matrix is retained privately; only its finite-local pilot joint measure is needed for this page.')
    manifest_path = output/'collision/manifest.json'
    manifest_path.parent.mkdir(exist_ok=True)
    manifest_path.write_bytes(json_bytes(manifest))
    for path in output.rglob('*'):
        if path.is_symlink():
            raise ValueError('Links are not capsule members')
        if path.is_file() and path.suffix != '.npz':
            text = path.read_text()
            if PRIVATE_PATH.search(text) or re.search(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}',text):
                raise ValueError('Unremoved private path/email:'+path.name)
    with archive.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed:
            with tarfile.open(fileobj=compressed,mode='w') as tar:
                for path in sorted(output.rglob('*')):
                    if path.is_file():
                        data = path.read_bytes()
                        info = tarfile.TarInfo('g8_why/'+path.relative_to(output).as_posix())
                        info.size,info.mode,info.mtime = len(data),0o644,0
                        tar.addfile(info,io.BytesIO(data))
    result = dict(schema='g8-collision-capsule-seal-v1',archive=digest(archive.read_bytes()),
        manifest=digest(manifest_path.read_bytes()),members=len(builder.files)+1,
        original_exported_accounting=True,privacy_scan_passed=True,
        native_invocations=0,new_grid_runs=0,process_cpu_s=time.process_time()-started)
    receipt.write_bytes(json_bytes(result))
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args = parser.parse_args()
    build(json.loads(args.input.read_bytes()),args.output,args.archive,args.receipt)
