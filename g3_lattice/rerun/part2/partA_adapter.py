"""Prepare collective-contrast fields only; no geometry/DFT/lattice dynamics.
Input strain_law is the producer's ABSOLUTE Biot strain relative to4.084A:
E_abs=(a_geometry/4.084)*(I+E_geometry)-I. Do NOT rescale it again.
DeltaK_lab=R[A(E_abs)-A((a_bulk_loading/4.084-1)I)]R.T.
Transfer to local Kii assumes only the EXTRA elastic delta(sumJ)=0;
composition/vacancy dependence belongs to the separately supplied loading law.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from lattice import CELL
from new_inputs.strain.strain_tensor_law import calibration_mask,stiffness_tensor


def collective_contrast(absolute_biot,rotation,model,bulk_a_A,*,extrapolate=False):
    """Material/lab collective tensors and the local-additive frozen-J contrast."""
    E=np.asarray(absolute_biot);R=np.asarray(rotation)
    reference=model['reference_a_A']
    bulk_strain=(bulk_a_A/reference-1)*np.eye(3)
    absolute=stiffness_tensor(E,model,extrapolate=extrapolate)
    bulk=stiffness_tensor(bulk_strain,model,extrapolate=extrapolate)
    contrast=np.einsum('...ia,...ab,...jb->...ij',R,absolute-bulk,R)
    lab=np.einsum('...ia,...ab,...jb->...ij',R,absolute,R)
    return contrast,lab,bulk


def reference_rows(contract):
    """Only actual delivered composition-reference rows; no guessed interpolation."""
    return [{'x':1.,'a_reference_A':contract['a1_A']}]+[
        {'x':r['x'],'a_reference_A':r['a_reference_A']} for r in contract['rows']]


def prepare_grid(source,output,model,reference,*,extrapolate=False):
    with np.load(source) as data:
        if not np.allclose(data['CELL_dimensionless'],CELL,rtol=0,atol=1e-12):
            raise ValueError('Registered FCC basis does not match the consumer')
        if not np.allclose(data['O_to_cubic_lab'],np.eye(3),rtol=0,atol=1e-12):
            raise ValueError('This delivered-map contract requires cubic Cartesian lab axes')
        # First two axes are Cartesian tensor indices; remaining axes are FCC i,j,k.
        absolute=np.moveaxis(data['strain_law'],(0,1),(-2,-1))
        geometry=np.moveaxis(data['strain_geometry'],(0,1),(-2,-1))
        rotation=np.moveaxis(data['R'],(0,1),(-2,-1))
        geometry_a=float(data['a_Pd_A'])
        expected=geometry_a/model['reference_a_A']*(np.eye(3)+geometry)-np.eye(3)
        if not np.allclose(absolute,expected,rtol=0,atol=1e-12):
            raise ValueError('Absolute Biot reference is not the declared4.084A producer contract')
        delta,lab,bulk=collective_contrast(absolute,rotation,model,reference['a_reference_A'],extrapolate=extrapolate)
        eligible=calibration_mask(absolute,model)
        output.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(output,delta_K_eV_A2=np.moveaxis(delta,(-2,-1),(0,1)),
            PartA_absolute_Gamma_lab_eV_A2=np.moveaxis(lab,(-2,-1),(0,1)),
            absolute_Biot_material=np.moveaxis(absolute,(-2,-1),(0,1)),R=data['R'],
            within_PartA_coordinate_box=eligible,CELL_dimensionless=data['CELL_dimensionless'],
            O_to_cubic_lab=data['O_to_cubic_lab'],a_Pd_A=data['a_Pd_A'],
            origin_D_A=data['origin_D_A'],reference_positions_A=data['reference_positions_A'],
            registered_positions_A=data['registered_positions_A'],
            bulk_loading_fraction=np.array(reference['x']),bulk_loading_a_A=np.array(reference['a_reference_A']),
            law_reference_a_A=np.array(model['reference_a_A']),PartA_bulk_loading_Gamma_material_eV_A2=bulk,
            extra_elastic_delta_sumJ_assumption=np.array('frozen at zero; loading-dependent J belongs to baseline'),
            operator=np.array('Kii_load + delta_K; NOT absolute collective Gamma as onsite'),
            strain_source=np.array('producer strain_law: absolute Biot U-I relative to4.084A'))
    return {'shape_tensor':list(np.moveaxis(delta,(-2,-1),(0,1)).shape),
        'bulk_loading_fraction':reference['x'],'bulk_loading_a_A':reference['a_reference_A'],
        'geometry_a_A':geometry_a,'law_reference_a_A':model['reference_a_A'],
        'fraction_outside_PartA_coordinate_box':float(1-np.mean(eligible)),
        'explicit_extrapolation':extrapolate,
        'tensor_symmetry_defect_eV_A2':float(np.max(abs(delta-np.swapaxes(delta,-1,-2)))),
        'maximum_absolute_contrast_eV_A2':float(np.max(abs(delta))),
        'source_SHA256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'output_SHA256':hashlib.sha256(output.read_bytes()).hexdigest(),
        'limits':'Conditional pure-metal EMT registered-D response; no hydrogen geometry or validated out-of-box hydride prediction; no clipping.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True,type=Path)
    ap.add_argument('--out',type=Path,default=Path('partA_prepared_maps'))
    ap.add_argument('--model',type=Path,default=Path('new_inputs/strain/strain-law.json'))
    ap.add_argument('--references',type=Path,default=Path('new_inputs/strain/loading-strain-reference-contract.json'))
    ap.add_argument('--extrapolate',action='store_true')
    args=ap.parse_args()
    model=json.loads(args.model.read_text());contract=json.loads(args.references.read_text())
    if model['reference_a_A']!=contract['a1_A']:
        raise ValueError('PartA and composition reference differ')
    rows=reference_rows(contract)
    # One annotated693 known case through THIS adapter: loading bulk must give zero.
    reference=next(r for r in rows if r['x']==.875)
    E=(reference['a_reference_A']/model['reference_a_A']-1)*np.eye(3)
    delta,absolute,bulk=collective_contrast(E,np.eye(3),model,reference['a_reference_A'],extrapolate=args.extrapolate)
    error=float(np.max(abs(delta)))
    if error!=0:raise ArithmeticError('Loading bulk acquired a false elastic contrast')
    wrong=float(np.max(abs(absolute-model['K0_eV_A2']*np.eye(3))))
    results=[]
    for reference in rows:
        label='x'+str(reference['x']).replace('.','p')
        for case in ['pd389','pd4084']:
            source=args.source/case/'fcc_octahedral_grid.npz'
            output=args.out/label/case/'fcc_octahedral_grid.npz'
            row=prepare_grid(source,output,model,reference,extrapolate=args.extrapolate)
            row.update(case=case,path=str(output.relative_to(args.out)))
            results.append(row)
    manifest={'source':'Giuliano547 registered FCC grid and640 PartA;693 composition reference',
        'operator':'delta_K=R[A(E_absolute)-A(E_bulk(a(x)))]R.T, added to loading local Kii; extra elastic deltaSumJ=0',
        'bulk_loading_reference_domain':contract['loading_reference_domain'],
        'loading_transfer_assumption':contract['local_transfer_assumption'],
        'registration':'Delivered FCC D-site grid; periodic i,j; open four-layer k, k=3 interface; actual reference/registered positions and polar frames retained',
        'known_case693':{'x':reference_rows(contract)[1]['x'],'correct_loading_bulk_max_abs_delta':error,
            'wrong_subtract_Azero_max_abs_eV_A2':wrong,'purpose':'Reject double-counting the composition-induced bulk curvature'},
        'law_SHA256':hashlib.sha256(args.model.read_bytes()).hexdigest(),
        'reference_contract_SHA256':hashlib.sha256(args.references.read_bytes()).hexdigest(),
        'results':results,'execution':'field algebra only; no dynamics, eigenmode search, geometry or DFT replay'}
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps(manifest,indent=2,allow_nan=False))


if __name__=='__main__':main()
