"""Package only new delivery and exact necessary external replay bindings."""
from pathlib import Path, PurePosixPath
import gzip
import hashlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREFIX = HERE.relative_to(ROOT)
OWN = ('README.md','method.md','adjoint.py','nonlinear.py','nonlinear_select.py',
       'nonlinear_verify.py','nonlinear_package.py','nonlinear_proof.json',
       'nonlinear_probe.json','nonlinear_states.npz','nonlinear_geometry.json',
       'nonlinear_selection.json','nonlinear_selection.npz')
BINDINGS = ('engine.py','patterns.py','exterior.py','large_geometry.py',
            'continuation_results.json','continuation_states.npz','input_manifest.json')


def package():
    from nonlinear import CASES
    payload, roles = {}, {}
    def add(path, role):
        relative = str(path.relative_to(ROOT))
        payload[relative] = path.read_bytes()
        roles[relative] = role
    for name in OWN:
        add(HERE/name,'owned delivery' if name != 'adjoint.py' else 'integrated peer source')
    for key in CASES:
        for suffix in ('.json','.npz'):
            add(HERE/('nonlinear_'+key+suffix),'new finite-model output')
    for name in BINDINGS:
        add(HERE/name,'immutable exact replay binding')
    inputs = json.loads((HERE/'input_manifest.json').read_text())['inputs']
    for item in inputs:
        if item['path'].startswith('g11_cmns/collective/contact/'):
            continue
        path = ROOT/item['path']
        content = path.read_bytes()
        assert len(content) == item['bytes'] and hashlib.sha256(content).hexdigest() == item['sha256']
        add(path,'immutable imported physical model input')
    records = [dict(path=p,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),role=roles[p])
               for p,b in payload.items()]
    manifest = dict(schema='pdd-nonlinear-delivery-manifest-v1',publication='unpublished local only',
        accepted_baseline='fcf949c21784331e44181ec81c5580cb084db51a',
        plan='coherencev1-g11-nonlinear-focusing-20261125',files=records,
        trace_columns=json.loads((HERE/'continuation_results.json').read_text())['trace_columns'],
        trace_sampling='every integration step; first unsupported sample retained as crossing evidence, never scored',
        state_units=dict(x='equilibrium-relative Cartesian Angstrom',v='Cartesian Angstrom/fs'),
        score_semantics=dict(K='full relative kinetic eV',inward='strictly positive closing radial kinetic eV',
            distance='reference R minus actual target R, Angstrom',loss_K='full relative kinetic eV under one 4100-fs amplitude decay'),
        supported='scored strictly before first all-site radius or original-NN separation exit; not measured angular coverage',
        excluded='all unrelated old results, contact/rate data, prior archives and checkpoints',
        external_replay='exact necessary source, accepted seed bytes and imported inputs only',
        archive='nonlinear.tar.gz',archive_hash='reported separately; no circular identity')
    encoded = (json.dumps(manifest,indent=2)+'\n').encode()
    (HERE/'nonlinear_manifest.json').write_bytes(encoded)
    payload[str(PREFIX/'nonlinear_manifest.json')] = encoded
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw,mode='w') as archive:
        for path,content in payload.items():
            entry = tarfile.TarInfo(path)
            entry.size,entry.mode,entry.mtime = len(content),0o644,0
            archive.addfile(entry,io.BytesIO(content))
    compressed = gzip.compress(raw.getvalue(),compresslevel=9,mtime=0)
    (HERE/'nonlinear.tar.gz').write_bytes(compressed)
    with tempfile.TemporaryDirectory(prefix='pdd-nonlinear-receiver-') as temp:
        receiver = Path(temp)
        with tarfile.open(fileobj=io.BytesIO(compressed),mode='r:gz') as archive:
            for item in archive.getmembers():
                path = PurePosixPath(item.name)
                assert item.isfile() and not path.is_absolute() and '..' not in path.parts
            archive.extractall(receiver,filter='data')
        assert all((receiver/path).read_bytes() == content for path,content in payload.items())
        # Cold imports and saved consumer readout in the actual extracted tree;
        # no duplicate simulation or original rate import.
        command = """import json,numpy as np; import nonlinear as n; from large_geometry import LargeLattice
from nonlinear_select import selected_result
lattices={}
for key,(m,e,dt,source) in n.CASES.items():
 r=json.loads((n.HERE/('nonlinear_'+key+'.json')).read_text())
 with np.load(n.HERE/('nonlinear_'+key+'.npz')) as a:
  if m not in lattices: lattices[m]=LargeLattice(m)
  lattice=lattices[m]
  for name,search in r['searches'].items():
   assert n.constraints(lattice,a[name+'_v'],e)['valid']
   final=search['final']; exit_fs=final['first_original_exit_fs']
   result=(None,None,a[name+'_trace'],-1. if exit_fs is None else exit_fs,{None:0,'site_radius':1,'bond_distance':2}[final['exit_kind']])
   assert n.summarize(result,'K' if name=='loss_K' else name,dt)==final
selected=json.loads((n.HERE/'nonlinear_selection.json').read_text())
with np.load(n.HERE/'nonlinear_selection.npz') as a:
 for key,case in selected.items():
  assert case['loss_comparison']['posthoc_damped_selected_conservative']['non_dynamical_diagnostic']
  with np.load(n.HERE/('nonlinear_'+key+'.npz')) as raw:
   for name,outcome in case['selected'].items():
    launch,result=selected_result(key,name,outcome,raw,a)
    assert n.constraints(lattices[n.CASES[key][0]],launch,case['energy_per_D_eV'])['valid']
    assert n.summarize(result,'K' if name=='loss_K' else name,case['dt_fs'])==outcome['readout']
print('cold extracted nonlinear consumer: exact output/constraints')"""
        subprocess.run([sys.executable,'-c',command],cwd=receiver/PREFIX,check=True)
    identity = dict(path=str(PREFIX/'nonlinear.tar.gz'),bytes=len(compressed),
        sha256=hashlib.sha256(compressed).hexdigest(),members=len(payload),
        verdict='safe fresh extraction, every archived byte equal, cold extracted saved consumer passed',
        publication='none; local unpublished capsule')
    (HERE/'nonlinear_archive_identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    print(json.dumps(identity,indent=2),flush=True)


if __name__ == '__main__':
    package()
