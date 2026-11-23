"""Saved-force bulk models. Units: eV, angstrom, atomic mass units.

K is local Kii, never Gamma. J uses means per occupied folded bond once.
Partial cells are registered, unrelaxed closures, not invented DFT structures.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import physical_constants, hbar, electron_volt, atomic_mass

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
KB = physical_constants['Boltzmann constant in eV/K'][0]
MASS = {'H': physical_constants['proton mass'][0]/atomic_mass,
        'D': physical_constants['deuteron mass'][0]/atomic_mass,
        'T': physical_constants['triton mass'][0]/atomic_mass}
PD = 'g3_loading/x1/dense_Pd/k3/final-x1-dense-readout-20261114T0213Z.json'
NI = 'g3_loading/nih_x1/dense_Ni/k3/readout-20261114T0917Z.json'
NIQ = 'g3_loading/nih_x1_quartic/dense_Ni/k3/quartic/readout-100-final-20261114T1001Z.json'
PARTIAL = [('x0875_spread','x0875-spread4','clean','0421'),
           ('x0875_spread','x0875-spread4','near','0500'),
           ('x0750_bunched','x0750-bunched8','clean','0608'),
           ('x0750_bunched','x0750-bunched8','near','0647'),
           ('x0875_bunched','x0875-bunched4','near','0746')]


def load(path):
    return json.loads((ROOT/path).read_text())


def kinetic_coefficient(mass_amu):
    """hbar^2/(2m), eV A^2; not the old hardcoded D coefficient."""
    return hbar*hbar/(2*mass_amu*atomic_mass*electron_volt*1e-20)


def folded_graph(positions, cell, shells):
    delta = positions[:,None]-positions[None,:]
    delta -= np.rint(delta/np.diag(cell))*np.diag(cell)
    distance = np.linalg.norm(delta,axis=-1)
    J = np.zeros(distance.shape)
    for radius, mean in shells:
        J[np.isclose(distance,radius,atol=2e-5,rtol=0)] = mean
    np.fill_diagonal(J,0.)
    return J


def full(material='Pd', variant='joint'):
    d = load(PD if material=='Pd' else NI)
    c, g = d['coupling'], d['profile']['geometry']
    ids = np.array(g['H_indices'])
    positions = np.array(g['reference_positions_A'])[ids]
    shells = [(s['distance_A'],s['sum_pulling_J_xx_eV_A2']/s['distinct_supercell_H_sites']) for s in c['shells']]
    k = c['on_site_k2_eV_A2']
    if material=='Pd':
        key = {'joint':'k4_eV_A4','energy':'K4_energy_only_eV_A4','force':'K4_force_only_eV_A4'}[variant]
        if variant=='joint':
            a,q = [d['directions'][axis][key] for axis in ('100','111')]
        else:
            summary=load('g3_loading/x1/dense-x1-summary.json')
            a,q=[summary['quartics'][axis][key] for axis in ('100','111')]
    else:
        qd=load(NIQ)['directions']['100']
        a=qd[{'joint':'K4_joint_eV_A4','energy':'K4_energy_eV_A4','force':'K4_force_eV_A4'}[variant]]
        old=load('g2_refit/computed_inputs.json')
        q=next(row['K4_111_eV_A4'] for row in old['materials'] if row['name']=='NiH')
    b=(3*q-a)/2
    return dict(name=material, K=np.tile(np.eye(3)*k,(len(ids),1,1)),
                f=np.zeros((len(ids),3)),a=a,b=b,J=folded_graph(positions,np.array(g['cell_A']),shells),
                positions=positions,original_indices=ids,site=int(np.where(ids==g['moving_index'])[0][0]),
                source=PD if material=='Pd' else NI,support_A=.24,variant=variant,
                axis=np.array([1.,0,0]),closure='scalar shell mean; cubic extension of source-x column',
                sumJ_source=c['sum_neighbour_pulling_J_xx_eV_A2'])


def partial_cases():
    base=full()
    cases=[]
    for directory,config,probe,stamp in PARTIAL:
        path=f'g3_loading/{directory}/loading_Pd/{config}/k3/readout-{probe}-20261114T{stamp}Z.json'
        d=load(path); c=d['coupling']; geom=d['profile']['geometry']; pr=d['probe']
        ids=np.flatnonzero(np.array(geom['atomic_numbers'])==1)
        positions=np.array(geom['reference_positions_A'])[ids]
        original=np.array(geom['original_atom_indices'])[ids]
        axis=np.array(pr['direction_vector']); site=int(np.where(ids==pr['moving_index'])[0][0])
        k=d['response']['K2_eV_A2']
        # Preserve the measured parallel curvature and signed F0 exactly. The
        # transverse curvature and all unprobed sites explicitly use x1 Kii.
        K=base['K'][:len(ids)].copy()
        K[site] += (k-base['K'][0,0,0])*np.outer(axis,axis)
        f=np.array(d['zero_forces_eV_A'])[ids]
        shells=[];unobserved=[]
        reference=load(PD)['coupling']['shells']
        for shell,ref in zip(c['shells'],reference):
            mean=shell['mean_pulling_Jdd_eV_A2']
            if mean is None:
                # Zero occupied neighbours supplies no per-bond measurement.
                # Other graph sites may have this shell: explicit x1 closure,
                # not zero and not an inferred partial-loading force constant.
                mean=ref['sum_pulling_J_xx_eV_A2']/ref['distinct_supercell_H_sites']
                unobserved.append(shell['radius_over_a'])
            shells.append((shell['distance_A'],mean))
        J=folded_graph(positions,np.array(geom['cell_A']),shells)
        cases.append(dict(name=config+'_'+probe,K=K,f=f,a=base['a'],b=base['b'],J=J,
                          positions=positions,original_indices=original,site=site,axis=axis,
                          source=path,support_A=.03,variant='transferred-x1',
                          closure='measured source axis; x1 transverse/unprobed K; measured all-site F0; source-shell J transferred symmetrically',
                          sumJ_source=c['sum_neighbour_Jdd_eV_A2'],vacancies=d['profile']['vacancy_original_indices'],
                          unobserved_source_shells_using_x1=unobserved))
    return cases


def strain_cases():
    base=full(); cases=[]
    # Preserve the force-column local baseline; transfer calibrated Gamma
    # CHANGES under frozen J. The k12 Gamma zero is not the k3 local Kii.
    reference=load('g3_skew/partA/strain-law.json')['K0_eV_A2']
    for name in ('g3_shear100_g004','g3_shear100_g008','g3_biax_m0p0200','g3_biax_0p0200'):
        path=f'g3_skew/partA/tensor-results/{name}.json';d=load(path)
        delta=np.array(d['K_lab_eV_A2'])-reference*np.eye(3)
        K=base['K'][0]+delta
        m=dict(base);m.update(name=name,K=np.tile(K,(32,1,1)),source=path,support_A=0.,
                            Gamma_reference_eV_A2=reference,Gamma_contrast_eV_A2=delta,
                            closure='Kii_bulk + calibrated [Gamma(strain)-Gamma(reference)] with frozen J; transferred x1 quartic; no strained local large-displacement data')
        cases.append(m)
    return cases


def potential(x,model):
    """x (...,sites,3); energy (...) and force with the same shape."""
    x=np.asarray(x); xx=x*x; r2=xx.sum(axis=-1,keepdims=True)
    quadratic=np.einsum('...ia,iab,...ib->...',x,model['K'],x)/2
    field=np.einsum('ij,...ja->...ia',model['J'],x)
    energy=quadratic-np.einsum('...ia,ia->...',x,model['f'])
    energy+=((model['a']-model['b'])/4*(xx*xx).sum(axis=-1)+model['b']/4*r2[...,0]**2).sum(axis=-1)
    energy-=np.sum(x*field,axis=(-2,-1))/2
    force=-np.einsum('iab,...ib->...ia',model['K'],x)+model['f']+field
    force-=((model['a']-model['b'])*xx+model['b']*r2)*x
    return energy,force


def serializable(model):
    return {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in model.items()}


def inventory():
    models=[full('Pd'),full('Ni')]+partial_cases()+strain_cases()
    rows=[]
    for m in models:
        H=np.block([[m['K'][i] if i==j else -m['J'][i,j]*np.eye(3) for j in range(len(m['K']))] for i in range(len(m['K']))])
        row=serializable(m)
        row.update(source_sha256=hashlib.sha256((ROOT/m['source']).read_bytes()).hexdigest(),
                   harmonic_min_eigenvalue_eV_A2=float(np.linalg.eigvalsh(H)[0]),
                   quartic_110_measured=False,quartic_min_direction_eV_A4=min(m['a'],(m['a']+2*m['b'])/3),
                   predicted_K4_110_eV_A4=(m['a']+3*(m['a']+2*m['b'])/3)/4)
        dependencies=[m['source'],'g3_loading/x1/dense-x1-summary.json']
        if m['name']=='Ni':dependencies += [NIQ,'g2_refit/computed_inputs.json']
        else:dependencies += [PD]
        if m['name'].startswith('g3_'):dependencies += ['g3_skew/partA/strain-law.json']
        row['coefficient_sources']=[dict(path=p,sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for p in dict.fromkeys(dependencies)]
        row['sampled_cuts']={'100':[-.03,.03,.12,.24],'111':[.12,.24]} if m['name']=='Pd' else (
            {'100':[-.03,.03,.24],'111':'older-settings g2_input quartic, not recomputed'} if m['name']=='Ni' else
            {'collective_Gamma_only':True,'local_displacement_points':[],'quartic':'transferred x1'} if m['name'].startswith('g3_') else
            {'measured_axis':[-.03,.03],'quartic':'transferred x1, not measured here'})
        rows.append(row)
    return rows
