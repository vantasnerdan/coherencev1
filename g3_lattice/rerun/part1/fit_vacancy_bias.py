"""Fit the authorized central bias law ONLY after all three zero-force cells exist.
No DFT, displacement dynamics, per-site vector interpolation, or b-only fit.
C2=0 is director819's assumption, never a measured/pseudoinverse coefficient.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from lattice import CELL
from measured_environment import build_template
from vacancy_bias import central_vacancy_features


def cell_rows(directory):
    zero_path=directory/'zero.json';config_path=directory/'config.json'
    zero=json.loads(zero_path.read_text());config=json.loads(config_path.read_text())
    if zero['profile']['config_id']!=config['config_id']:
        raise ValueError('Zero-force/configuration identity mismatch')
    native=build_template(zero,config,CELL)
    occupied=native['occupied'];original=native['original_atom_index']
    features,contract=central_vacancy_features(occupied,CELL,(1,2,3))
    geometry=zero['profile']['geometry']
    measured_rows=[i for i,z in enumerate(geometry['atomic_numbers']) if z==1]
    picks=[tuple(np.argwhere(original==geometry['original_atom_indices'][i])[0]) for i in measured_rows]
    design_all=np.concatenate([features[(slice(None),slice(None))+index].T for index in picks],axis=0)
    target=np.asarray(zero['forces_eV_A'])[measured_rows]
    return {'id':config['config_id'],'loading_x':zero['profile']['loading_x'],'a_reference_A':zero['profile']['a_A'],
        'original_H_indices':[geometry['original_atom_indices'][i] for i in measured_rows],
        'X':design_all[:,[0,2]],'Y':target.reshape(-1),'Y_vectors':target,
        'second_shell_feature_max_abs':float(np.max(abs(design_all[:,1]))),
        'zero_SHA256':hashlib.sha256(zero_path.read_bytes()).hexdigest(),
        'config_SHA256':hashlib.sha256(config_path.read_bytes()).hexdigest(),'contract':contract}


def solve(cells):
    X=np.concatenate([c['X'] for c in cells]);Y=np.concatenate([c['Y'] for c in cells])
    beta,_,rank,singular=np.linalg.lstsq(X,Y,rcond=None)
    return beta,int(rank),singular


def errors(residual):
    residual=np.asarray(residual).reshape(-1,3)
    return {'D_vectors':len(residual),'vector_RMS_eV_A':float(np.sqrt(np.mean(np.sum(residual**2,axis=1)))),
        'component_RMS_eV_A':float(np.sqrt(np.mean(residual**2))),
        'maximum_vector_residual_eV_A':float(np.max(np.linalg.norm(residual,axis=1)))}


def fit(cases):
    if len(cases)!=3:raise ValueError('All THREE b/c/d cells are required; no provisional one-cell fit')
    cells=[cell_rows(p) for p in cases]
    if len({c['id'] for c in cells})!=3:raise ValueError('Three distinct measured cells required')
    beta,rank,singular=solve(cells)
    if rank!=2:raise ValueError('Shell1/3 coefficients are not jointly identifiable from supplied cells')
    per_cell=[]
    for cell in cells:
        per_cell.append(dict(case_id=cell['id'],loading_x=cell['loading_x'],a_reference_A=cell['a_reference_A'],
            **errors(cell['X']@beta-cell['Y'])))
    total=errors(np.concatenate([c['X']@beta-c['Y'] for c in cells]))
    heldout=[]
    for index,cell in enumerate(cells):
        training=[other for j,other in enumerate(cells) if j!=index]
        b,r,s=solve(training)
        row={'held_out_cell':cell['id'],'training_cells':[c['id'] for c in training],
            'training_rank':r,'training_singular_values':s.tolist()}
        if r==2:
            row.update(coefficients13_eV_A=b.tolist(),**errors(cell['X']@b-cell['Y']))
        else:
            row.update(coefficients13_eV_A=None,vector_RMS_eV_A=None,component_RMS_eV_A=None,
                reason='Training cells do not identify the two coefficients; no pseudoinverse-zero accuracy is claimed')
        heldout.append(row)
    return {'schema':'central-vacancy-bias-v1','authority':'Axiom812/819',
        'force_law':'F0_i=sum_s beta_s sum_{vacant physical reference images j in shell s} unit(R_ij)',
        'shell_coefficients_eV_A':{'1':float(beta[0]),'2':0.,'3':float(beta[1])},
        'second_shell_policy':'FIXED ZERO ASSUMPTION by819;2a-periodic source cells cancel the odd-force feature and cannot measure beta2',
        'fitted_shells':[1,3],'fit_rank':rank,'fit_singular_values':singular.tolist(),
        'weighting':'Equal weight per Cartesian force component, each measured D ONCE; no duplicate primitive-cover rows',
        'RMS_definition':'vector_RMS=sqrt(sum_i ||Fpred-Fmeas||²/N_D); component_RMS=vector_RMS/sqrt(3)',
        'fit_error':total,'per_cell_error':per_cell,'leave_one_cell_out':heldout,
        'accuracy_scope':'Calibration/held-out errors on these cells, NOT guaranteed arbitrary-arrangement accuracy or confidence bounds',
        'source_state':'Registered unrelaxed sites; law is the constant F0 term, not a relaxed minimum or a new displacement-dependent spring',
        'strain_transfer':'Directions are reference lattice lines; extra strain dependence of F0 is unmeasured and frozen',
        'sources':[{'case_id':c['id'],'loading_x':c['loading_x'],'a_reference_A':c['a_reference_A'],
            'H_vectors':len(c['Y_vectors']),'original_H_indices':c['original_H_indices'],
            'zero_SHA256':c['zero_SHA256'],'config_SHA256':c['config_SHA256'],
            'second_shell_feature_max_abs':c['second_shell_feature_max_abs']} for c in cells],
        'feature_contract':cells[0]['contract'],
        'single_loading_sensitivity':'One non2a-periodic large occupancy; SAME seed/other settings; beta2=0,+beta1/2,-beta1/2. Scenarios, not calibrated bounds; no additional cases on this point.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',type=Path,action='append',required=True)
    ap.add_argument('--out',type=Path,default=Path('vacancy_bias_law.json'))
    args=ap.parse_args();result=fit(args.case)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
