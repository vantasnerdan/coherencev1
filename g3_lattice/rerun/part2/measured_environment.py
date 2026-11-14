"""Exact source-periodic occupancy/F0 input, NOT a vector bias law.
The supplied zero-SCF geometry is the displacement origin, not a relaxed minimum.
Only the actual measured configuration (and integer periodic copies) is supported.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


FORCE_FRAME='Cartesian cubic [100],[010],[001]'


def build_template(zero,config,CELL):
    profile=zero['profile'];geometry=profile['geometry'];a=profile['a_A']
    positions=np.asarray(geometry['reference_positions_A'],float)
    if not np.allclose(positions,np.asarray(zero['positions_A']),rtol=0,atol=1e-12):
        raise ValueError('Zero-SCF and registered reference positions differ')
    if not all(geometry['pbc']) or not np.allclose(geometry['cell_A'],2*a*np.eye(3),rtol=0,atol=1e-12):
        raise ValueError('This actual-data adapter requires the delivered periodic2x2x2 conventional cell')
    if geometry['original_atom_indices']!=config['geometry']['original_atom_indices']:
        raise ValueError('Original-to-retained atom identities differ')
    atomic=np.asarray(geometry['atomic_numbers']);H=np.flatnonzero(atomic==1)
    original=np.asarray(geometry['original_atom_indices'],int)
    forces=np.asarray(zero['forces_eV_A'],float)
    if forces.shape!=positions.shape or not np.isfinite(forces).all():
        raise ValueError('Full zero-SCF Cartesian force array is required')
    coordinates=2*positions[H]/a
    rounded=np.rint(coordinates).astype(int)
    if np.max(abs(coordinates-rounded))>1e-10:
        raise ValueError('Measured sites do not lie on the registered half-grid')
    lookup={tuple(key%4):int(row) for key,row in zip(rounded,H)}
    if len(lookup)!=len(H):raise ValueError('Duplicate retained H identity in source cell')
    origin=positions[H[0]].copy();origin2=np.rint(2*origin/a).astype(int)
    basis2=2*np.asarray(CELL)
    if not np.array_equal(basis2,np.rint(basis2)):
        raise ValueError('Unexpected FCC basis')
    index=np.indices((4,4,4))
    keys=origin2[:,None,None,None]+np.einsum('ab,bijk->aijk',basis2.astype(int),index)
    occupied=np.zeros((4,4,4),bool);bias=np.zeros((3,4,4,4))
    original_grid=np.full((4,4,4),-1,int);retained_grid=np.full((4,4,4),-1,int)
    for ijk in np.ndindex((4,4,4)):
        row=lookup.get(tuple(keys[(slice(None),)+ijk]%4))
        if row is not None:
            occupied[ijk]=True;bias[(slice(None),)+ijk]=forces[row]
            original_grid[ijk]=original[row];retained_grid[ijk]=row
    counts=[int(np.sum(retained_grid==row)) for row in H]
    if counts!=[2]*len(H):raise ValueError('Primitive cover did not preserve each source H exactly twice')
    return dict(occupied=occupied,force_bias_eV_A=bias,original_atom_index=original_grid,
        retained_atom_index=retained_grid,CELL_dimensionless=np.asarray(CELL),origin_D_A=origin,
        registered_a_A=np.array(a),loading_fraction=np.array(profile['loading_x']),
        force_frame=np.array(FORCE_FRAME),config_id=np.array(profile['config_id']),
        profile_id=np.array(profile['profile_id']),vacancy_original_indices=np.array(config['vacancy_original_indices']),
        source_H_atoms=np.array(len(H)),source_cell_A=np.asarray(geometry['cell_A']),
        template_scope=np.array('Two exact conventional-cell copies on4x4x4 FCC primitive grid; NOT arbitrary-arrangement/strain bias prediction'))


def load_measured_environment(path,shape,expected_loading,CELL,shift=(0,0,0)):
    with np.load(path) as source:
        occupied=source['occupied'];force=source['force_bias_eV_A']
        original=source['original_atom_index']
        if force.shape!=(3,)+occupied.shape or occupied.dtype!=np.bool_ or not np.isfinite(force).all():
            raise ValueError('Measured occupancy and complete Cartesian bias must share one native grid')
        if source['force_frame'].item()!=FORCE_FRAME or not np.array_equal(source['CELL_dimensionless'],CELL):
            raise ValueError('Bias/frame/FCC convention mismatch')
        if float(source['loading_fraction'])!=expected_loading:
            raise ValueError('Measured template and requested loading differ')
        if np.any(force[:,~occupied]!=0):
            raise ValueError('Vacant sites are absent DOFs, not measured force-bearing sites')
        if len(shape)!=3 or any(n%m for n,m in zip(shape,occupied.shape)):
            raise ValueError('Target must contain integer periods of the measured occupancy/bias template')
        if len(shift)!=3 or any(int(x)!=x for x in shift):
            raise ValueError('Measured template shift is an integer primitive translation')
        metadata=dict(config_id=source['config_id'].item(),profile_id=source['profile_id'].item(),
            registered_a_A=float(source['registered_a_A']),loading_fraction=float(source['loading_fraction']),
            origin_D_A=source['origin_D_A'].tolist(),primitive_shift=list(shift),
            source_H_atoms=int(source['source_H_atoms']),native_occupied=int(occupied.sum()),
            scope=source['template_scope'].item(),force_frame=FORCE_FRAME,
            origin='Actual zero-SCF registered sites; NOT relaxed/recentered',
            source_SHA256=hashlib.sha256(Path(path).read_bytes()).hexdigest())
        repetitions=tuple(n//m for n,m in zip(shape,occupied.shape))
        occupied=np.tile(np.roll(occupied,shift,axis=(0,1,2)),repetitions)
        force=np.tile(np.roll(force,shift,axis=(1,2,3)),(1,)+repetitions)
        original=np.tile(np.roll(original,shift,axis=(0,1,2)),repetitions)
    return occupied,force,original,metadata


def normalized_bias(force_eV_A,length_A,Eref_eV):
    """Force conversion is length/Eref, NOT the curvature factor length²/Eref."""
    if length_A<=0 or Eref_eV<=0:raise ValueError('Positive physical reference units required')
    return np.asarray(force_eV_A)*(length_A/Eref_eV)


def main():
    from lattice import CELL
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('zero',type=Path);ap.add_argument('config',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args()
    zero=json.loads(args.zero.read_text());config=json.loads(args.config.read_text())
    fields=build_template(zero,config,CELL)
    fields['zero_source_SHA256']=np.array(hashlib.sha256(args.zero.read_bytes()).hexdigest())
    fields['config_source_SHA256']=np.array(hashlib.sha256(args.config.read_bytes()).hexdigest())
    np.savez_compressed(args.output,**fields)
    print(json.dumps({'config_id':zero['profile']['config_id'],'source_H_atoms':int(fields['source_H_atoms']),
        'native_occupied':int(fields['occupied'].sum()),'native_sites':int(fields['occupied'].size),
        'source_origin':'registered unrelaxed zero-SCF','output_SHA256':hashlib.sha256(args.output.read_bytes()).hexdigest()},indent=2))


if __name__=='__main__':main()
