"""Actual864 baseline response; composition spacing is not extra elastic strain."""
import json
from hashlib import sha256
from math import exp
from pathlib import Path
import numpy as np
from new_inputs.loading_reference import loading_lattice as reference_spacing
from new_inputs.loading_law_864.loading_law import evaluate


def bulk_reference(x,*,extrapolate=False):
    x=float(x)
    if not .5<=x<=1:raise ValueError('Lattice slider is0.5<=x<=1; fitted-source range starts.75')
    if x>=.75:return reference_spacing.a(x)
    if not extrapolate:raise ValueError('Below.75 requires explicit868 continuation of the same reference construction')
    boundary=reference_spacing.a(.75)  # Also authenticates the pinned calibration facts.
    facts=json.loads(reference_spacing.SOURCE.read_text())
    return dict(boundary,x=x,a_angstrom=facts['author_extrapolated_a_x1_A']*exp(facts['dln_a_dx']*(x-1)),
        status='model_extrapolated_below_campaign_data',measured_point=False,
        construction='868-authorized continuation of a=4.084*exp[.044*(x-1)]; uncertainty unknown')


def response_fields(p,occupied,vacant_nn,law_path,a_reference_A,*,extrapolate=False):
    raw=Path(law_path).read_bytes();law=json.loads(raw)
    if law['model_id']!='baseline_two_descriptor':raise ValueError('868 selects the exact two-term baseline, not motif')
    A=p['K4_100_eV_A4'];B=(3*p['K4_111_eV_A4']-A)/2
    zp=(3*A+2*B)*p['zero_point_axis_variance_A2']
    K=np.zeros(occupied.shape)
    J=np.zeros((len(p['shells']),)+occupied.shape)
    outside={name:0 for name in law['models']};table=[]
    for count in np.unique(vacant_nn[occupied]):
        mask=occupied&(vacant_nn==count)
        answer=evaluate(int(count),a_reference_A,law,extrapolate=extrapolate)
        values=answer['values_eV_A2'];number=int(mask.sum())
        K[mask]=values['Kii_dd']+zp
        for index,shell in enumerate(p['shells']):
            # The unit kernel sums to ONE. Use reference folded multiplicity,
            # NOT the measured occupied sum or a second occupancy factor.
            J[index,mask]=shell['reference_multiplicity']*values[shell['loading_observable']]
        for name in answer['extrapolated_observables']:outside[name]+=number
        table.append(dict(n_vac_NN=int(count),occupied_sites=number,**answer))
    total=int(occupied.sum())
    metadata=dict(source_entry='864',model_id=law['model_id'],law_SHA256=sha256(raw).hexdigest(),
        a_reference_A=a_reference_A,descriptor_counts=table,
        extrapolated_fraction_by_observable={name:number/total for name,number in outside.items()},
        fraction_Nvac_above4=float(np.mean(vacant_nn[occupied]>4)),
        direction_closure=law['direction_closure'],uncertainty=law['uncertainty'],
        bond_policy='Per-existing folded-bond mean times reference multiplicity; physical-image sharing by shell_kernel; occupancy applied ONCE; reciprocal endpoint average',
        zero_point_policy='Frozen x1 per-axis T0 Hartree only; thermal motion sampled, not added twice',
        absent_sites='No degree of freedom; zero field storage here is not a measured zero bond/well')
    # At the full x1 anchor the existing constant kernel is exactly identical.
    if occupied.all() and np.all(vacant_nn==0) and a_reference_A==law['reference_a_A']:
        J=None
    return K,A,B,J,metadata


def elastic_tensor(strain,bulk_a_A,law_path,*,map_path=None,extrapolate=False):
    """Compose extra hydro strain with actual Biot U; do not add nonlinear contrasts twice."""
    from partA_adapter import collective_contrast
    from new_inputs.strain.strain_tensor_law import calibration_mask
    model=json.loads(Path(law_path).read_text());shape=strain.shape;identity=np.eye(3)
    if map_path is not None:
        with np.load(map_path) as source:
            if abs(float(source['bulk_loading_a_A'])-bulk_a_A)>1e-10:
                raise ValueError('Prepared map does not subtract this composition reference')
            E=np.moveaxis(source['absolute_Biot_material'],(0,1),(-2,-1))
            R=np.moveaxis(source['R'],(0,1),(-2,-1))
            direct=source['delta_K_eV_A2'].copy()
        if E.shape!=shape+(3,3) or R.shape!=E.shape:
            raise ValueError('Registered map frame/array shape differs from the dynamical sites')
        if not np.any(strain):
            delta=direct
        else:
            E=(1+strain)[...,None,None]*(identity+E)-identity
            contrast,_,_=collective_contrast(E,R,model,bulk_a_A,extrapolate=extrapolate)
            delta=np.moveaxis(contrast,(-2,-1),(0,1))
    elif np.any(strain):
        E=(bulk_a_A/model['reference_a_A']*(1+strain)-1)[...,None,None]*identity
        contrast,_,_=collective_contrast(E,identity,model,bulk_a_A,extrapolate=extrapolate)
        delta=np.moveaxis(contrast,(-2,-1),(0,1))
    else:
        return None,dict(kind='no additional elastic contrast')
    inside=calibration_mask(E,model)
    return delta,dict(kind='PartA absolute-minus-composition contrast',reference_a_A=model['reference_a_A'],
        bulk_loading_a_A=bulk_a_A,fraction_outside_PartA_coordinate_box=float(np.mean(~inside)),
        extra_elastic_delta_sumJ='frozen zero',composition='(1+extra_hydro)*(I+absolute_Biot)-I; R unchanged',
        geometry_scope='reference-coordinate dynamics; map geometry is conditional pure-metal EMT, not a physical hydride equilibrium')
