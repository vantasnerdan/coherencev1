"""One local unpublished capsule; old results/paths stay pinned prerequisites."""
import gzip
import hashlib
import json
import tarfile
import well

HERE=well.HERE
ROOT=HERE.parents[1]
META={'capsule.json','capsule.tar.gz','ARCHIVE_FILES.txt','SHA256SUMS','prerequisites.json'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def record(path):
    return dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path))

def main():
    prerequisites=[well.PRIOR/name for name in ('well.py','solve.py','coupled_sensitivity.py',
        'inputs/old-readout.json','inputs/bulk-wells.json','results/summary_main.json','results/summary_check.json')]
    prerequisites += [well.PRIOR/f'results/D_{T}_{v}_{tag}_distributions.npz'
        for T in (300,77) for v in ('old','quadratic','quartic') for tag in ('main','check')]
    prerequisites += [ROOT/p for p in ('g10_quantum_site/model.py','g10_quantum_site/coupled_momentum.py',
        'g10_quantum_site/output/coupled/trajectory_summary.json',
        'g10_quantum_site/output/coupled/PdD_300K_paths.npz','g10_quantum_site/output/coupled/PdD_77K_paths.npz',
        'g10_quantum_site/output/coupled/NiH_300K_paths.npz','g10_quantum_site/output/coupled/NiH_77K_paths.npz',
        'g3_loading/x1/dense_Pd/k3/final-x1-dense-readout-20261114T0213Z.json',
        'g3_loading/nih_x1/dense_Ni/k3/readout-20261114T0917Z.json',
        'g3_loading/nih_x1_quartic/dense_Ni/k3/quartic/readout-100-final-20261114T1001Z.json',
        'g2_refit/computed_inputs.json')]
    prerequisites += [ROOT/p for p in ('g10_quantum_site/exact.py',
        'g10_quantum_site/output/exact/NiH_77K_joint.npz',
        'g10_quantum_site/output/exact/NiH_300K_joint.npz')]
    prerequisite_records=[record(p) for p in sorted(set(prerequisites))]
    script=next(p for p in prerequisite_records if p['path']=='g10_quantum_site/bulk-refit-048/coupled_sensitivity.py')
    script.update(
        bytes=16751,sha256='5a7b84c13be1171adf5603a64012754bafbef1897744c4ae20849a42a79c00ee',
        role='exact pre-portability source consumed for the scientific result; available in prior asset 789149873',
        accepted_replay_equivalent=dict(
            bytes=16709,sha256='10988fafc1061660f234101526288d885d5705af3bcc66caf1f40b9f45ddb476',
            accepted_main='4f1fd691af7d797b4bea7dfedc8103154d0b36cf',
            difference='one reporting-string substitution: (an interpreter path under the author home directory) -> python; callbacks/operator unchanged, as observed by Main',
            identity_basis='SHA/size derived from the consumed source by the exact single substitution reported by Main',
            replay_selection='keep the accepted portable script; do not overwrite it from the prior capsule'))
    owned=[];generated=[]
    for path in sorted(HERE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.name in META:continue
        if path.suffix not in ('.py','.md','.csv','.json','.png','.npz'):
            raise ValueError(f'Unclassified delivery member {path}')
        if path.stat().st_size>1_000_000:
            if path.suffix!='.npz':raise ValueError(f'Large non-array member needs explicit decision: {path}')
            generated.append(record(path));continue
        owned.append(path)
    provenance=dict(source='authenticated Giuliano Fleet entry_3684; asset 789418823',
        source_archive_bytes=94239,source_archive_sha256='8085c5623dfbbf137a097249e6637c42182565dcb3760ee544e36df1f90553f4',
        signed_json_sha256=sha(well.INPUT/'bulk-wells.json'),
        prior_accepted_revision='ed06115',repository_base='7349c26e626d41262f9ad66ab94e9e266373425c',
        external_prerequisites=prerequisite_records,
        prior_archive_overlay=dict(
            asset_id=789149873,archive_filename='g10-pdd-refit-048.tar.gz',
            archive_sha256='52c9b8fc0973d274244949d00ac1eadddac9997b153f8359cfaf023e8e48d9aa',
            asset_url='https://api.github.com/repos/vantasnerdan/coherencev1/releases/assets/789149873',
            accepted_main='4f1fd691af7d797b4bea7dfedc8103154d0b36cf',
            members=[f'g10_quantum_site/bulk-refit-048/results/D_{T}_{v}_main_distributions.npz'
                     for T in (77,300) for v in ('old','quadratic','quartic')],
            scope='extract only these six arrays into a separate replay tree; no regeneration, old/current mutation or copies in this new capsule',
            instructions='replay.md: Prerequisites without duplicating old data'),
        regenerable_generated_outputs=generated,
        regeneration='Exact portable Python commands are in replay.md and nih/README.md.',
        rule='external prerequisite bytes are not duplicated; >1 MB generated NPZ arrays are omitted and regenerable, not required for supplied table/figure consumption')
    (HERE/'prerequisites.json').write_text(json.dumps(provenance,indent=2)+'\n')
    files=sorted(owned+[HERE/'prerequisites.json']+[p for p in well.INPUT.rglob('*') if p.is_file()],key=lambda p:str(p.relative_to(ROOT)))
    prefix=HERE.relative_to(ROOT)
    names=[str(p.relative_to(ROOT)) for p in files]+[str(prefix/'ARCHIVE_FILES.txt'),str(prefix/'SHA256SUMS')]
    (HERE/'ARCHIVE_FILES.txt').write_text('\n'.join(names)+'\n')
    (HERE/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(ROOT)}\n' for p in files+[HERE/'ARCHIVE_FILES.txt']))
    # Freeze metadata and gzip timestamp: same bytes yield the same local archive.
    with (HERE/'capsule.tar.gz').open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed:
            with tarfile.open(fileobj=compressed,mode='w') as tar:
                for path in files+[HERE/'ARCHIVE_FILES.txt',HERE/'SHA256SUMS']:
                    info=tar.gettarinfo(str(path),arcname=str(path.relative_to(ROOT)))
                    info.uid=info.gid=0;info.uname=info.gname='';info.mtime=0
                    with path.open('rb') as f:tar.addfile(info,f)
    capsule=dict(filename='capsule.tar.gz',sha256=sha(HERE/'capsule.tar.gz'),bytes=(HERE/'capsule.tar.gz').stat().st_size,
        payload_files=len(names),archive_members='ARCHIVE_FILES.txt',payload_hashes='SHA256SUMS',
        input_archive_sha256=provenance['source_archive_sha256'],selected_signed_json_sha256=provenance['signed_json_sha256'],
        external_prerequisites='prerequisites.json',governing_plan='coherencev1-g10-bulk080-20261125',
        template_revision='14c5716f18cd5b1b8ee77624079ed31a0879c64388a625f4db88162dda02faf4',
        permission='unpublished local capsule only; no push, publication, release mutation or public delivery')
    (HERE/'capsule.json').write_text(json.dumps(capsule,indent=2)+'\n')
    print(json.dumps(capsule,indent=2))

if __name__=='__main__':
    main()
