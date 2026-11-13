#!/usr/bin/env python3
"""B/C only: EMT metal geometry and explicitly provisional collective-Gamma contrast.
Run: .venv/bin/python mismatch.py relax --case pd389 --steps 1200
Then: .venv/bin/python mismatch.py maps
No hydrogen energies are evaluated by EMT. See README.md for conventions.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
import resource
import time
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.neighborlist import neighbor_list
from ase.optimize import FIRE

ROOT = Path(__file__).resolve().parent
CASES = {'pd389': (3.89, 10, 11), 'pd4084': (4.084, 6, 7)}
A_NI = 3.52
A_LAW = 4.084
FMAX = 0.02

def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def build(case):
    a, n, m = CASES[case]
    L = n * a
    positions, symbols, layers = [], [], []
    for species, count, scale, zbase, start in (
            ('Pd', n, a, 10.0, 0),
            ('Ni', m, L / m, 10 + 5*a/2 + (a+A_NI)/4, 6)):
        for layer in range(6 if species == 'Pd' else 4):
            basis = ((0, 0), (.5, .5)) if layer % 2 == 0 else ((.5, 0), (0, .5))
            for ix in range(count):
                for iy in range(count):
                    for bx, by in basis:
                        positions.append(((ix+bx)*scale, (iy+by)*scale,
                                          zbase + layer*(a if species == 'Pd' else A_NI)/2))
                        symbols.append(species)
                        layers.append(start + layer)
    height = max(p[2] for p in positions) + 10
    atoms = Atoms(symbols, positions=positions, cell=(L, L, height), pbc=(True, True, False))
    atoms.set_array('layer', np.array(layers))
    atoms.set_constraint(FixAtoms(mask=np.array(layers) < 2))
    return atoms

def relax(case, steps):
    out = ROOT / case
    out.mkdir(exist_ok=True)
    reference_path = out / 'reference.extxyz'
    if reference_path.exists():
        reference = read(reference_path)
    else:
        reference = build(case)
        write(reference_path, reference)
    final_path = out / 'relaxed.extxyz'
    atoms = read(final_path) if final_path.exists() else reference.copy()
    atoms.set_constraint(FixAtoms(mask=atoms.arrays['layer'] < 2))
    atoms.calc = EMT()
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    initial_energy = atoms.get_potential_energy()
    optimizer = FIRE(atoms, logfile=str(out / 'relax.log'), maxstep=0.10)
    converged = optimizer.run(fmax=FMAX, steps=steps)
    forces = atoms.get_forces()
    write(final_path, atoms)
    record = {
        'case': case, 'atoms': len(atoms), 'Pd_atoms': int(sum(atoms.numbers == 46)),
        'Ni_atoms': int(sum(atoms.numbers == 28)), 'Pd_a_A': CASES[case][0],
        'Ni_natural_a_A': A_NI, 'Ni_imposed_a_A': atoms.cell[0, 0]/CASES[case][2],
        'Ni_residual_inplane_strain': atoms.cell[0, 0]/CASES[case][2]/A_NI - 1,
        'cell_A': atoms.cell.array.tolist(), 'pbc': atoms.pbc.tolist(),
        'frame': ['[100]', '[010]', '[001]'], 'layers_Pd_Ni': [6, 4],
        'fixed_Pd_layers_from_bottom': [0, 1], 'registry': 'zero lateral offset; alternating conventional fcc(001) bases; no twist',
        'optimizer': 'ASE FIRE, maxstep=0.10 A; fixed cell', 'fmax_target_eV_A': FMAX,
        'converged': bool(converged), 'steps_this_segment': optimizer.nsteps,
        'max_mobile_force_eV_A': float(np.linalg.norm(forces, axis=1).max()),
        'initial_energy_this_segment_eV': float(initial_energy),
        'final_energy_eV': float(atoms.get_potential_energy()),
        'wall_seconds': time.perf_counter()-start_wall,
        'cpu_seconds': time.process_time()-start_cpu,
        'peak_process_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        'threads': 1, 'status': 'local minimum candidate only; no global-minimum or deposition-kinetics claim'}
    history_path = out / 'relaxation.json'
    history = json.loads(history_path.read_text()) if history_path.exists() else []
    history.append(record)
    dump(history_path, history)
    print(json.dumps(record), flush=True)

def local_F(reference, atoms):
    """Reference-neighbour affine least squares; preserve each reference image branch."""
    pd = reference[reference.numbers == 46]
    current = atoms.positions[atoms.numbers == 46]
    a = float(reference.cell[0, 0])/round(np.sqrt(len(pd)/12))
    ii, jj, ss = neighbor_list('ijS', pd, .80*a)
    X = pd.positions[jj] + ss @ pd.cell.array - pd.positions[ii]
    Y = current[jj] + ss @ atoms.cell.array - current[ii]
    F, residual, count = [], [], []
    for i in range(len(pd)):
        select = ii == i
        x, y = X[select], Y[select]
        f = np.linalg.lstsq(x, y, rcond=None)[0].T
        if np.linalg.matrix_rank(x) != 3 or np.linalg.det(f) <= 0:
            raise ValueError('Invalid local reference geometry or inverted fitted deformation')
        F.append(f)
        residual.append(np.sqrt(np.mean(np.sum((y-x@f.T)**2, axis=1))))
        count.append(len(x))
    return np.array(F), np.array(residual), np.array(count), (ii, jj, ss)

def polar(F):
    u, s, vt = np.linalg.svd(F)
    R = u @ vt
    stretch = np.einsum('...ji,...j,...jk->...ik', vt, s, vt)
    return R, stretch

def law(F, a_pd, source):
    """Normal-only cubic closure in the local reference crystal; shear coefficient UNKNOWN."""
    row = next(v for v in source['model_inputs'] if v['material'] == 'PdD')
    k0 = row['K_zero_eV_A2']
    A = row['strain']['parallel']['dK_depsilon_eV_A2']
    B = row['strain']['perpendicular']['dK_depsilon_eV_A2']
    R, U = polar(F)
    strain_geometry = U - np.eye(3)
    strain_law = U * (a_pd/A_LAW) - np.eye(3)
    diagonal = np.diagonal(strain_law, axis1=-2, axis2=-1)
    kval = k0 + B*diagonal.sum(axis=-1)[..., None] + (A-B)*diagonal
    K = np.einsum('...ik,...k,...jk->...ij', R, kval, R)
    offset = a_pd/A_LAW-1
    kbulk_matched = k0 + (A+2*B)*offset
    return {'F': F, 'R': R, 'strain_geometry': strain_geometry, 'strain_law': strain_law,
            'K_Gamma_eV_A2': K, 'DeltaK_law_bulk_eV_A2': K-k0*np.eye(3),
            'DeltaK_matched_bulk_eV_A2': K-kbulk_matched*np.eye(3),
            'K_matched_bulk_eV_A2': np.array(kbulk_matched)}

def spectral_fields(K):
    values, vectors = np.linalg.eigh(K)
    gap = values[..., 1]-values[..., 0]
    return {'eigenvalues_eV_A2': values, 'lowest_eigenvalue_eV_A2': values[..., 0],
            'lowest_direction_lab': vectors[..., :, 0],
            'lowest_direction_degenerate': gap < 1e-6}

def affine_check():
    angle = .37
    R = np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    U = np.array([[1.02,.012,0],[.012,.99,0],[0,0,1.005]])
    target = R@U
    reference = build('pd4084')
    deformed = reference.copy()
    deformed.positions = reference.positions@target.T + np.array([3.,-2.,7.])
    deformed.set_cell(reference.cell.array@target.T, scale_atoms=False)
    got, residual, _, _ = local_F(reference, deformed)
    rg, ug = polar(got)
    rigid = reference.copy()
    rigid.positions = reference.positions@R.T + np.array([3.,-2.,7.])
    rigid.set_cell(reference.cell.array@R.T, scale_atoms=False)
    rigid_F, _, _, _ = local_F(reference, rigid)
    rr, ur = polar(rigid_F)
    errors = {'affine_F_max_error': float(abs(got-target).max()),
              'stretch_max_error': float(abs(ug-U).max()),
              'rotation_max_error': float(abs(rg-R).max()),
              'rigid_strain_max_error': float(abs(ur-np.eye(3)).max())}
    errors['affine_neighbour_residual_A'] = float(residual.max())
    if max(errors.values()) > 1e-12:
        raise ValueError(errors)
    return errors

def octahedral_contract(reference, atoms, metal_ids, Fs):
    """Complementary fcc octahedral sublattice, four-coplanar-metal interpolation."""
    a, L = CASES[reference.info['case']][0], reference.cell[0,0]
    pd_ref = reference.positions[metal_ids]
    layer = reference.arrays['layer'][metal_ids]
    dref = pd_ref.copy()
    dref[:,0] = (dref[:,0]+a/2) % L
    displacement = atoms.positions[metal_ids]-pd_ref
    supports, dpos, dF = [], [], []
    for i, point in enumerate(dref):
        delta = pd_ref-point
        delta[:,:2] -= np.rint(delta[:,:2]/L)*L
        dist = np.linalg.norm(delta, axis=1)
        dist[layer != layer[i]] = np.inf
        neighbours = np.argsort(dist)[:4]
        if not np.allclose(dist[neighbours], a/2, atol=1e-8):
            raise ValueError('Octahedral registration does not have four nearest coplanar metals')
        supports.append(neighbours)
        dpos.append(point+displacement[neighbours].mean(axis=0))
        dF.append(Fs[neighbours].mean(axis=0))
    return dref, np.array(dpos), np.array(dF), np.array(supports)

def pattern_metrics(xy, fields, L):
    values = fields-fields.mean(axis=0)
    modes = []
    for mx in range(-4,5):
        for my in range(-4,5):
            if mx == 0 and my == 0:
                continue
            amp = np.mean(values*np.exp(-2j*np.pi*(mx*xy[:,0]+my*xy[:,1])/L)[:,None,None], axis=0)
            modes.append((float(np.linalg.norm(amp)), mx, my))
    amp, mx, my = max(modes)
    return {'rms_modulation_Frobenius': float(np.sqrt(np.mean(np.sum(values**2,axis=(1,2))))),
            'dominant_mode_mx_my': [mx,my], 'dominant_mode_amplitude': amp,
            'dominant_mode_wavelength_A': float(L/np.hypot(mx,my))}

def morphology(xy, selected, a, n):
    """Threshold-set nearest-metal graph on a periodic plane; winding rank is explicit."""
    q = np.rint(2*xy/a).astype(int) % (2*n)
    lookup = {tuple(v): i for i,v in enumerate(q)}
    remaining = set(np.flatnonzero(selected).tolist())
    components = []
    for root in list(remaining):
        if root not in remaining:
            continue
        lifts = {root: q[root].copy()}
        remaining.remove(root)
        queue, windings = [root], []
        for i in queue:
            for step in (np.array([1,1]),np.array([1,-1]),np.array([-1,1]),np.array([-1,-1])):
                j = lookup[tuple((q[i]+step) % (2*n))]
                if not selected[j]:
                    continue
                candidate = lifts[i]+step
                if j not in lifts:
                    lifts[j] = candidate
                    remaining.remove(j)
                    queue.append(j)
                elif np.any(candidate != lifts[j]):
                    windings.append((candidate-lifts[j])//(2*n))
        rank = int(np.linalg.matrix_rank(np.array(windings))) if windings else 0
        components.append({'sites':len(queue),'winding_rank':rank})
    rank = max((v['winding_rank'] for v in components),default=-1)
    kind = {-1:'none',0:'isolated spots or finite patches',1:'periodic lines',2:'percolating two-dimensional grid/region'}[rank]
    if np.all(selected):
        kind = 'entire plane, not isolated soft lines'
    return {'fraction':float(np.mean(selected)), 'classification':kind,
            'components': sorted(components,key=lambda v:-v['sites'])}

def export_fcc(out, d, a, n):
    """Register octahedral sites to the integration owner's FCC primitive layer grid."""
    CELL = np.array([[1,-1,1],[1,1,0],[0,0,1]],dtype=float)/2
    k = d['Pd_layer_from_bottom']-2
    x, y = d['reference_positions_A'][:,:2].T/a
    raw_ij = np.column_stack((x-.5+y-k/2, y-x+.5+k/2))
    if np.max(abs(raw_ij-np.rint(raw_ij))) > 1e-8:
        raise ValueError('Octahedral/FCC registration is not integral')
    ij = np.rint(raw_ij).astype(int)%n
    keys = ['K_Gamma_eV_A2','DeltaK_law_bulk_eV_A2','DeltaK_matched_bulk_eV_A2',
            'F','R','strain_geometry','strain_law','shear_magnitude_geometry']
    grid = {}
    disagreement = {}
    for key in keys:
        values = d[key]
        tail = values.shape[1:]
        field = np.empty(tail+(n,n,4))
        maximum = 0.
        for i in range(n):
            for j in range(n):
                for z in range(4):
                    rows = (ij[:,0]==i)&(ij[:,1]==j)&(k==z)
                    if np.sum(rows) != 2:
                        raise ValueError('Expected exactly two conventional-box copies')
                    pair = values[rows]
                    maximum = max(maximum,float(np.max(abs(pair[0]-pair[1]))))
                    field[(...,i,j,z)] = pair.mean(axis=0)
        if maximum > 1e-7:
            raise ValueError(f'Cannot fold inequivalent primitive copies: {key} {maximum}')
        grid[key] = field
        disagreement[key] = maximum
    indices = np.indices((n,n,4))
    origin = np.array([a/2,0,10+a])
    grid['reference_positions_A'] = np.einsum('ab,bijk->aijk',a*CELL,indices)+origin[:,None,None,None]
    displacement = d['registered_positions_A']-d['reference_positions_A']
    disp_grid = np.empty((3,n,n,4))
    for i in range(n):
        for j in range(n):
            for z in range(4):
                rows = (ij[:,0]==i)&(ij[:,1]==j)&(k==z)
                disp_grid[:,i,j,z] = displacement[rows].mean(axis=0)
    grid['registered_positions_A'] = grid['reference_positions_A']+disp_grid
    grid['CELL_dimensionless'] = CELL
    grid['origin_D_A'] = origin
    grid['a_Pd_A'] = np.array(a)
    grid['O_to_cubic_lab'] = np.eye(3)
    grid['K_matched_bulk_eV_A2'] = d['K_matched_bulk_eV_A2']
    np.savez_compressed(out/'fcc_octahedral_grid.npz',**grid)
    return {'shape_tensor': [3,3,n,n,4], 'origin_D_A':origin.tolist(),
            'CELL_dimensionless':CELL.tolist(),'O_to_cubic_lab':np.eye(3).tolist(),
            'k_to_original_Pd_layer':'k+2; k=3 is interface; k=0 deepest mobile plane',
            'duplicate_copy_max_abs_disagreement':disagreement,
            'registration':'r_D=origin_D+aPd*CELL@(i,j,k); periodic i,j of length n; open k',
            'source':'octahedral_maps.npz, not metal_maps.npz'}

def maps():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    source_path = ROOT/'adopted-strains.json'
    source = json.loads(source_path.read_text())
    row = next(v for v in source['model_inputs'] if v['material'] == 'PdD')
    k0 = row['K_zero_eV_A2']
    A = row['strain']['parallel']['dK_depsilon_eV_A2']
    B = row['strain']['perpendicular']['dK_depsilon_eV_A2']
    summary = {'law_source': 'Fleet entry_464', 'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
               'law_reference_a_A': A_LAW, 'K0_Gamma_eV_A2': k0,
               'normal_closure_A_B_eV_A2': [A,B],
               'predicted_hydro_derivative_eV_A2': A+2*B,
               'source_hydro_derivative_eV_A2': row['strain']['hydro']['dK_depsilon_eV_A2'],
               'units': {'coordinates':'angstrom','strain':'dimensionless Biot stretch minus identity','stiffness':'eV/angstrom^2'},
               'crystal_and_lab_frame':['[100]','[010]','[001]'],
               'site_kind':{'metal_maps.npz':'Pd metal nuclei','octahedral_maps.npz':'registered octahedral interstitial centers, NOT metal nuclei'},
               'frozen_J_contrast_contract':'For source bulk a0=4.084 A: Kii_effective=Kii_source+DeltaK_law_bulk, J unchanged. DeltaK_matched_bulk isolates contrast about same-aPd cubic baseline; background compression is separately explicit.',
               'missing_physics':['independent cubic shear-response coefficient','strained onsite/inter-site split','hydride metal geometry response','part A law'],
               'shear_response_coefficient': None, 'affine_rigid_motion_check': affine_check(), 'cases': {}}
    fig, axes = plt.subplots(2,4,figsize=(11,5.4),constrained_layout=True)
    softfig, softaxes = plt.subplots(2,4,figsize=(11,5.4),constrained_layout=True)
    depthfig, daxes = plt.subplots(1,2,figsize=(8,3.4),constrained_layout=True)
    for ci, (case,(a,n,m)) in enumerate(CASES.items()):
        out = ROOT/case
        reference, atoms = read(out/'reference.extxyz'), read(out/'relaxed.extxyz')
        reference.info['case'] = case
        F_all, residual, neighbours, graph = local_F(reference,atoms)
        all_pd = np.flatnonzero(reference.numbers == 46)
        layer_all = reference.arrays['layer'][all_pd]
        chosen = layer_all >= 2
        ids, Fs = all_pd[chosen], F_all[chosen]
        layer = layer_all[chosen]
        result = law(Fs,a,source)
        result.update(spectral_fields(result['K_Gamma_eV_A2']))
        result.update({'metal_atom_index': ids, 'Pd_layer_from_bottom': layer,
                       'depth_layer_from_interface': 5-layer,
                       'reference_positions_A': reference.positions[ids],
                       'relaxed_positions_A': atoms.positions[ids],
                       'nonaffinity_RMS_A': residual[chosen], 'reference_neighbor_count': neighbours[chosen]})
        ii,jj,ss = graph
        np.savez_compressed(out/'reference_bonds.npz', Pd_atom_index=all_pd,
                            center_Pd_index=ii, neighbor_Pd_index=jj, periodic_image_shift=ss)
        dref,dpos,dF,support = octahedral_contract(reference,atoms,ids,Fs)
        d = law(dF,a,source)
        d.update(spectral_fields(d['K_Gamma_eV_A2']))
        d.update({'reference_positions_A': dref, 'registered_positions_A': dpos,
                  'metal_support_row_indices': support, 'metal_support_atom_indices': ids[support],
                  'metal_support_weights': np.full(support.shape,.25),
                  'Pd_layer_from_bottom': layer, 'depth_layer_from_interface': 5-layer})
        d['shear_magnitude_geometry'] = np.sqrt(2*np.sum(d['strain_geometry'][:,np.triu_indices(3,1)[0],np.triu_indices(3,1)[1]]**2,axis=1))
        np.savez_compressed(out/'octahedral_maps.npz', **d)
        fcc_contract = export_fcc(out,d,a,n)
        Kmin = result['lowest_eigenvalue_eV_A2']
        kmatched = float(result['K_matched_bulk_eV_A2'])
        layers = []
        for li in range(6):
            select = layer_all == li
            _,stretch = polar(F_all[select])
            metric = pattern_metrics(reference.positions[all_pd[select],:2],stretch-np.eye(3),a*n)
            metric.update({'layer_from_interface': 5-li,'nominal_depth_A': (5-li)*a/2,
                           'fixed': li < 2, 'mean_strain_geometry': (stretch-np.eye(3)).mean(axis=0).tolist(),
                           'nonaffinity_RMS_A': float(np.sqrt(np.mean(residual[select]**2)))})
            layers.append(metric)
        for li in range(4):
            mask = (5-layer) == li
            xy = reference.positions[ids[mask],:2]
            shown_depth = 0 if li < 2 else 3
            shown = (5-layer) == shown_depth
            shown_xy = reference.positions[ids[shown],:2]
            strain = result['strain_geometry'][shown]
            shear = np.sqrt(2*np.sum(strain[:,np.triu_indices(3,1)[0],np.triu_indices(3,1)[1]]**2,axis=1))
            strain_scalar = shear if li % 2 else np.trace(strain,axis1=1,axis2=2)
            label = 'shear Frobenius' if li % 2 else 'tr(U-I)'
            artist=axes[ci,li].scatter(shown_xy[:,0],shown_xy[:,1],c=strain_scalar,s=9,cmap='coolwarm')
            fig.colorbar(artist,ax=axes[ci,li],label=label)
            axes[ci,li].set(title=f'{case}, depth {shown_depth}',xlabel='[100] Å',ylabel='[010] Å',aspect='equal')
            artist=softaxes[ci,li].scatter(xy[:,0],xy[:,1],c=Kmin[mask],s=9,cmap='viridis',vmin=Kmin.min(),vmax=Kmin.max())
            softfig.colorbar(artist,ax=softaxes[ci,li],label='Kmin eV/A²')
            softaxes[ci,li].set(title=f'{case}, Pd depth {li}',xlabel='[100] Å',ylabel='[010] Å',aspect='equal')
        layers.sort(key=lambda x:x['layer_from_interface'])
        depth = np.array([v['nominal_depth_A'] for v in layers])
        rms = np.array([v['rms_modulation_Frobenius'] for v in layers])
        daxes[0].plot(depth,rms,'o-',label=case)
        daxes[1].plot(depth,[v['nonaffinity_RMS_A'] for v in layers],'o-',label=case)
        low = int(Kmin.argmin())
        normal_law = np.diagonal(result['strain_law'],axis1=1,axis2=2)
        result['shear_magnitude_geometry'] = np.sqrt(2*np.sum(result['strain_geometry'][:,np.triu_indices(3,1)[0],np.triu_indices(3,1)[1]]**2,axis=1))
        np.savez_compressed(out/'metal_maps.npz', **result)
        threshold_morphology = []
        for depth_index in range(4):
            select = (5-layer)==depth_index
            xy = reference.positions[ids[select],:2]
            threshold_morphology.append({
                'depth_layer_from_interface':depth_index,
                'half_law_bulk':morphology(xy,Kmin[select]<k0/2,a,n),
                'half_matched_bulk':morphology(xy,Kmin[select]<kmatched/2,a,n),
                'octahedral_half_law_bulk':morphology(dref[select,:2],d['lowest_eigenvalue_eV_A2'][select]<k0/2,a,n)})
        details = {'metal_sites': len(ids),'octahedral_sites': len(dref),
                   'K_matched_cubic_bulk_eV_A2': kmatched,
                   'uniform_law_reference_linear_offset': a/A_LAW-1,
                   'Kmin_range_eV_A2': [float(Kmin.min()),float(Kmin.max())],
                   'fraction_below_half_law_bulk': float(np.mean(Kmin<k0/2)),
                   'fraction_below_half_matched_bulk': float(np.mean(Kmin<kmatched/2)),
                   'fraction_negative': float(np.mean(Kmin<0)),
                   'octahedral_fraction_below_half_law_bulk': float(np.mean(d['lowest_eigenvalue_eV_A2']<k0/2)),
                   'octahedral_Kmin_range_eV_A2': [float(d['lowest_eigenvalue_eV_A2'].min()),float(d['lowest_eigenvalue_eV_A2'].max())],
                   'normal_law_strain_component_range': [float(normal_law.min()),float(normal_law.max())],
                   'fraction_outside_1percent_normal_component_box': float(np.mean(np.max(abs(normal_law),axis=1)>.01)),
                   'max_abs_geometry_shear_Biot': float(np.max(abs(result['strain_geometry'][:,np.triu_indices(3,1)[0],np.triu_indices(3,1)[1]]))),
                   'degenerate_lowest_fraction': float(np.mean(result['lowest_direction_degenerate'])),
                   'softest_metal_atom_index': int(ids[low]),'softest_reference_position_A': reference.positions[ids[low]].tolist(),
                   'softest_direction_lab': result['lowest_direction_lab'][low].tolist(),
                   'imposed_square_period_A': a*n,'coincidence_primitive_diagonal_translation_A': a*n/np.sqrt(2),
                   'layers': layers,'mobile_layer3_to_interface_modulation_ratio': float(rms[3]/rms[0]),
                   'threshold_morphology_by_depth': threshold_morphology,
                   'fcc_consumer_contract': fcc_contract,
                   'penetration_interpretation': 'finite slab, bottom two layers fixed; layer-resolved RMS is evidence, not a bulk penetration constant',
                   'relaxation': json.loads((out/'relaxation.json').read_text())}
        summary['cases'][case]=details
        print(case,json.dumps({k:v for k,v in details.items() if k not in ('layers','relaxation')}),flush=True)
    for ax in daxes:
        ax.set(xlabel='Nominal depth below top Pd layer (Å)')
        ax.legend()
    daxes[0].set(ylabel='Strain modulation RMS (dimensionless)')
    daxes[1].set(ylabel='Neighbour affine residual RMS (Å)')
    fig.savefig(ROOT/'strain_maps.png',dpi=130)
    softfig.savefig(ROOT/'stiffness_maps.png',dpi=130)
    depthfig.savefig(ROOT/'penetration.png',dpi=130)
    plt.close('all')
    dump(ROOT/'summary.json',summary)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['relax','maps'])
    parser.add_argument('--case',choices=CASES,default='pd389')
    parser.add_argument('--steps',type=int,default=1200)
    args=parser.parse_args()
    if args.action=='relax':
        relax(args.case,args.steps)
    else:
        maps()
