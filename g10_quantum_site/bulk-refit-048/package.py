"""Freeze only this refit's named payload and its exact coupled replay inputs."""
from pathlib import Path
import hashlib
import json
import tarfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PREFIX=HERE.relative_to(ROOT)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    owned=[]
    for path in sorted(HERE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts:continue
        if path.name in ('capsule.json','capsule.tar.gz','ARCHIVE_FILES.txt','SHA256SUMS'):continue
        if '_coarse_' in path.name or path.name=='summary_coarse.json':continue
        owned.append(path)
    deps=[ROOT/'g10_quantum_site'/name for name in ('model.py','coupled_momentum.py','output/coupled/trajectory_summary.json','output/coupled/PdD_300K_paths.npz','output/coupled/PdD_77K_paths.npz')]
    deps.append(ROOT/'g3_loading/x1/dense_Pd/k3/final-x1-dense-readout-20261114T0213Z.json')
    files=sorted(owned+deps,key=lambda p:str(p.relative_to(ROOT)))
    entries=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
    (HERE/'ARCHIVE_FILES.txt').write_text('\n'.join(v['path'] for v in entries)+f'\n{PREFIX}/ARCHIVE_FILES.txt\n{PREFIX}/SHA256SUMS\n')
    (HERE/'SHA256SUMS').write_text(''.join(f'{v["sha256"]}  {v["path"]}\n' for v in entries)+f'{sha(HERE/"ARCHIVE_FILES.txt")}  {PREFIX}/ARCHIVE_FILES.txt\n')
    with tarfile.open(HERE/'capsule.tar.gz','w:gz') as tar:
        for path in files+[HERE/'ARCHIVE_FILES.txt',HERE/'SHA256SUMS']:
            info=tar.gettarinfo(str(path),arcname=str(path.relative_to(ROOT)))
            info.uid=info.gid=0;info.uname=info.gname='';info.mtime=0
            with path.open('rb') as f:tar.addfile(info,f)
    capsule=dict(filename='capsule.tar.gz',sha256=sha(HERE/'capsule.tar.gz'),bytes=(HERE/'capsule.tar.gz').stat().st_size,payload_files=len(files)+2,payload_manifest='SHA256SUMS',archive_members='ARCHIVE_FILES.txt',input_archive_sha256='d488fd0e11c5b2c82cd8a438e6174e49a56353b627bd7054641c7c51f133a74f',selected_signed_json_sha256=sha(HERE/'inputs/bulk-wells.json'),governing_plan='coherencev1-g10-pdd-refit048-20261124',template_revision='14c5716f18cd5b1b8ee77624079ed31a0879c64388a625f4db88162dda02faf4',pinned_worktree='7349c26e626d41262f9ad66ab94e9e266373425c',permission='unpublished local capsule only; no push or release')
    (HERE/'capsule.json').write_text(json.dumps(capsule,indent=2)+'\n')
    print(json.dumps(capsule,indent=2))


if __name__=='__main__':main()
