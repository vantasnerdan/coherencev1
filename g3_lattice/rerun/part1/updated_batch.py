"""One changed-input PdD batch; Ni may be appended only with its adopted input.
All displacements are relative to declared reference sites, not DFT-relaxed minima.
"""
import argparse
import copy
import json
import os
import time
from hashlib import sha256
from pathlib import Path
import numpy as np
from lattice import OBSERVATION_DEFAULTS,run_case,regions
from coherence import characterize
from map_consumer import load_native_map,correlations

ROOT=Path(__file__).resolve().parent
LAW='new_inputs/loading_law_864/loading-law-base.json'
BIAS='new_inputs/vacancy_bias_law.json'


def prepare_parameters(destination):
    parameters=json.loads((ROOT/'lattice_inputs.json').read_text())
    quantum=json.loads((ROOT/'pd_x1_quantum'/'result.json').read_text())
    law=json.loads((ROOT/LAW).read_text())
    for name in ('PdH','PdD'):
        row=next(item for item in quantum['rows'] if item['material']==name)
        values=row['material_inputs'];p=parameters['materials'][name]
        p.update(a_A=values['a_A'],mass_eVc2=values['mass_eVc2'],Kii_eV_A2=values['K2_eV_A2'],
            K4_100_eV_A4=values['K4_100_eV_A4'],K4_111_eV_A4=values['K4_111_eV_A4'],
            zero_point_axis_variance_A2=row['spectrum']['local_mean_square_axis_A2_T0'],
            reference_amplitude_decay_ps=values['reference_amplitude_decay_ps'])
        p.pop('collective_strain_slope_eV_A2')
        for shell in p['shells']:
            radius=float(np.linalg.norm(shell['periodic_alias_offsets_A'][0])/p['a_A'])
            matches=[(key,item) for key,item in law['anchor']['shell_inventory'].items()
                if abs(item['radius_over_a']-radius)<1e-8]
            if len(matches)!=1:raise ValueError('Actual spatial shell cannot be matched to864 inventory')
            key,item=matches[0]
            if len(shell['periodic_alias_offsets_A'])!=item['reference_multiplicity']:
                raise ValueError('Folded source-bond multiplicity changed')
            shell.update(J_sum_eV_A2=item['sum_Jdd_eV_A2'],loading_observable=key,
                reference_multiplicity=item['reference_multiplicity'])
        if abs(sum(s['J_sum_eV_A2'] for s in p['shells'])-values['J_sum_eV_A2'])>1e-12:
            raise ValueError('Same-column J sum was not preserved')
        p['status']='Adopted762 fixed-k3 local well and same-column J;864 directional scalar transfer; x1 quartic/T0 variance frozen, NOT mesh-converged; physical decay is a model input'
    parameters['loading_law']=dict(kind='anchored_directional_affine_v2',hosts=['PdD','PdH'],input=LAW,
        strain_law_input='new_inputs/strain/strain-law.json',source_entry='864',selected='two-term baseline by868',
        no_extra_screening=True,weak_shell_scope='All identified full-fit coefficients retained; poor sparse configuration holdouts are not predictive bounds')
    parameters['model']='Actual762 local cubic well +864 folded-shell response +819 central reference F0 +640 PartA contrast; parametric Langevin reference-coordinate dynamics'
    destination.write_text(json.dumps(parameters,indent=2,allow_nan=False)+'\n')
    return parameters


def prepare_map(name,destination):
    source_path=ROOT/'partA_prepared_maps'/'x0p875'/name/'fcc_octahedral_grid.npz'
    _,metadata=load_native_map(source_path)
    with np.load(source_path) as source:
        native=source['delta_K_eV_A2'].shape[2:]
        repeat=(60//native[0],60//native[1],1)
        if native[0]*repeat[0]!=60 or native[1]*repeat[1]!=60:
            raise ValueError('Actual native map cannot tile this declared square index cover')
        data={key:np.tile(source[key],(1,1)+repeat) for key in ('delta_K_eV_A2','absolute_Biot_material','R')}
        data['within_PartA_coordinate_box']=np.tile(source['within_PartA_coordinate_box'],repeat)
        for key in ('CELL_dimensionless','O_to_cubic_lab','a_Pd_A','bulk_loading_fraction','bulk_loading_a_A','law_reference_a_A'):
            data[key]=source[key].copy()
        data['source_native_shape']=np.array(native);data['repeat']=np.array(repeat)
    np.savez_compressed(destination,**data)
    return dict(metadata,source_native_shape=list(native),repeat=list(repeat),shape=list(data['delta_K_eV_A2'].shape[2:]),
        source_SHA256=sha256(source_path.read_bytes()).hexdigest())


def pd_cases(folder):
    cases=[]
    common=dict(material='PdD',L=32,seed=13,temperature_K=0.,strain=0.,pump_ratio=1.,
        pump_polarization=[1.,0.,0.],loading_extrapolate=True,elastic_extrapolate=True,
        bias_law_input=BIAS,bias_geometry_policy='frozen_reference_under_strain',**OBSERVATION_DEFAULTS)
    def add(label,**updates):
        case=dict(common,label=label);case.update(updates);cases.append(case)
    for x in (1.,.875,.75,.6):
        for h in (0.,.08):add(f'random_x{x:g}_h{h:g}',x=x,arrangement='random',pump_h=h)
    for arrangement in ('ordered','clustered'):
        for h in (0.,.08):add(f'{arrangement}_x0.875_h{h:g}',x=.875,arrangement=arrangement,pump_h=h)
    for h in (0.,.08):add(f'stretch_x0.875_h{h:g}',x=.875,arrangement='random',pump_h=h,strain=.05)
    add('gradient_x0.875_h0.08',x=.875,arrangement='random',pump_h=.08,strain_gradient_per_nm=.02)
    for h in (0.,.08):add(f'thermal_x0.875_h{h:g}',x=.875,arrangement='random',pump_h=h,temperature_K=300.)
    add('detuned_x0.875_h0.08',x=.875,arrangement='random',pump_h=.08,pump_ratio=1.1)
    for source,x in (('loading_b_793',.875),('loading_c_827',.75),('loading_d_850',.875)):
        for kind in ('measured','fitted'):
            add(f'{source}_{kind}',x=x,L=24,arrangement='measured_periodic',pump_h=.08,
                measured_environment_input=f'new_inputs/{source}/measured_environment.npz')
            if kind=='measured':cases[-1].pop('bias_law_input')
    for name in ('pd389','pd4084'):
        target=folder/(name+'_tensor.npz');metadata=prepare_map(name,target)
        for h in (0.,.08):
            add(f'map_{name}_h{h:g}',x=.875,L=metadata['shape'],arrangement='random',pump_h=h,
                geometry_a_A=metadata['geometry_a_A'],boundary_z='open',tensor_input=str(target.relative_to(ROOT)),
                tensor_kind='delta_collective_curvature',map_metadata=metadata)
    for multiplier in (.5,-.5):
        add(f'C2_{multiplier:+g}',x=.875,arrangement='random',pump_h=.08,bias_shell2_multiplier=multiplier)
    if len(cases)!=30:raise AssertionError('Declared representative case set changed')
    return cases


def field_comparison(left,right,key):
    with np.load(left) as first,np.load(right) as second:
        mask=first['occupied']
        if not np.array_equal(mask,second['occupied']):raise ValueError('Comparison changed the occupancy')
        a=first[key][:,mask];b=second[key][:,mask]
    denominator=float(np.linalg.norm(a));other=float(np.linalg.norm(b));inner=np.vdot(a,b)
    aligned=b*np.exp(-1j*np.angle(inner)) if np.iscomplexobj(a) and inner!=0 else b
    return dict(RMS_difference_A=float(np.sqrt(np.mean(np.sum(abs(a-b)**2,axis=0)))),
        relative_norm_difference=float(np.linalg.norm(a-b)/denominator) if denominator else None,
        global_phase_invariant_overlap=float(abs(inner)/(denominator*other)) if denominator and other else None,
        phase_aligned_relative_difference=float(np.linalg.norm(a-aligned)/denominator) if denominator else None)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--group',choices=('pd','ni'),default='pd');ap.add_argument('--parameters',type=Path)
    ap.add_argument('--out',type=Path,default=Path('updated'))
    args=ap.parse_args();os.chdir(ROOT);folder=args.out.resolve();folder.mkdir(parents=True,exist_ok=True)
    if args.parameters is not None:parameters=json.loads(args.parameters.read_text())
    elif args.group=='pd':parameters=prepare_parameters(folder/'parameters.json')
    else:raise ValueError('Ni requires its actual adopted input; no provisional substitution')
    if args.group=='pd':
        cases=pd_cases(folder)
        control=dict(cases[0],label='actual_x1_harmonic_consumer',L=24,steps=64,harmonic_control=True,
            scalability_only=True,pump_h=0.)
        witness,_=run_case(parameters,control)
        (folder/'harmonic_consumer.json').write_text(json.dumps(witness,indent=2,allow_nan=False)+'\n')
        print(json.dumps({'consumer_error':witness['harmonic_control_max_coordinate_error']}),flush=True)
    else:
        cases=[dict(label=f'NiH_updated_h{h:g}',material='NiH',L=24,x=1.,arrangement='random',seed=13,
            pump_h=h,pump_polarization=[1.,0.,0.],temperature_K=0.,**OBSERVATION_DEFAULTS) for h in (0.,.08)]
    (folder/(args.group+'_cases.json')).write_text(json.dumps({'cases':cases},indent=2)+'\n')
    results=[];started=time.perf_counter()
    for case in cases:
        result,fields=run_case(parameters,case)
        a=case.get('geometry_a_A',result['loading_response']['composition_reference']['a_angstrom']
            if result['loading_response']['composition_reference'] is not None else parameters['materials'][case['material']]['a_A'])
        result['regions']=regions(fields,a)
        if case.get('boundary_z','periodic')=='periodic':
            result['half_pump_correlations'],_=characterize(fields['driven_amplitude_A'],fields['occupied'],a)
            result['correlation_filter']='Whole late-record joint DC+Omega+2Omega amplitude, NOT the earlier16-cycle equal-time envelope'
        else:
            period=case['map_metadata']['source_native_shape'][:2]
            with np.load(case['tensor_input']) as tensor_source:
                result['map_response']=correlations(fields,tensor_source['delta_K_eV_A2'],period)
        if case['label'].startswith('C2_'):
            base=next(row for row in results if row['label']=='random_x0.875_h0.08')
            if result['occupation_SHA256']!=base['occupation_SHA256']:
                raise AssertionError('819 sensitivity changed its fixed occupancy')
        path=folder/(case['label']+'.npz');np.savez_compressed(path,**fields)
        results.append(result)
        (folder/(args.group+'_results.json')).write_text(json.dumps(dict(results=results,
            elapsed_s=time.perf_counter()-started,scope='Finite-time reference-model response; no physical material-phase certificate'),indent=2,allow_nan=False)+'\n')
        print(json.dumps({key:result[key] for key in ('label','elapsed_s','state','static_rms_A','driven_demodulated_rms_A','pump_frequency_rms_A','peak_fraction')}),flush=True)
    comparisons={}
    if args.group=='pd':
        for source in ('loading_b_793','loading_c_827','loading_d_850'):
            comparisons[source]={key:field_comparison(folder/(source+'_measured.npz'),folder/(source+'_fitted.npz'),key)
                for key in ('static_mean_A','driven_amplitude_A')}
        for multiplier in (.5,-.5):
            comparisons[f'C2_{multiplier:+g}']={key:field_comparison(folder/'random_x0.875_h0.08.npz',folder/(f'C2_{multiplier:+g}.npz'),key)
                for key in ('static_mean_A','driven_amplitude_A')}
        for label in ('random_x1','random_x0.875','random_x0.75','random_x0.6','ordered_x0.875','clustered_x0.875',
                'stretch_x0.875','thermal_x0.875','map_pd389','map_pd4084'):
            comparisons[label+'_pump_addition']={key:field_comparison(folder/(label+'_h0.npz'),folder/(label+'_h0.08.npz'),key)
                for key in ('static_mean_A','driven_amplitude_A')}
    (folder/(args.group+'_comparisons.json')).write_text(json.dumps(comparisons,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'completed_group':args.group,'cases':len(results),'elapsed_s':time.perf_counter()-started}),flush=True)


if __name__=='__main__':main()
