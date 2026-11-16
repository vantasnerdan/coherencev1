"""Complete signed host-minus-neutral-vacuum cycle; explicit numerical eligibility."""
import argparse
import json
from pathlib import Path
import numpy as np
from read_fingerprint import bound_file
from screen_io import manifest,require,write_json


def read_unit(directory, expected_host):
    directory=Path(directory)
    complete=json.loads((directory/'completed.json').read_bytes())
    require(complete['schema']=='coherence-neutral-cycle-completed-v1','Wrong completion role')
    producer=json.loads(bound_file(directory,complete['producer']).read_bytes())
    plan=json.loads(bound_file(directory,producer['input']).read_bytes())
    for source in producer['producer_source']:bound_file(directory,source)
    expected={c['case']:c for c in plan['cases']}
    data={}
    for reference in complete['records']:
        energy=json.loads(bound_file(directory,reference['energy']).read_bytes())
        config=energy['config'];case=energy['case']
        require(case==reference['case'] and config==expected[case], 'Energy/config source mismatch')
        require(config['host'] in [expected_host,'vacuum'],'Wrong requested host')
        require(energy['actual_neutral_valence_electrons']==config['valence_electrons'],
                'Endpoint has wrong neutral electron count')
        if config['host']=='vacuum':
            symbol=config['endpoint'];g=config['geometry']
            require(symbol in ['H','He'] and g['symbols']==[symbol],
                    'Reference must be isolated neutral atoms, never H2')
            require(g['initial_moments_muB']==[1. if symbol=='H' else 0.],
                    'Wrong neutral vacuum spin reference')
            require(config['valence_electrons']==(1 if symbol=='H' else 2),
                    'Wrong neutral vacuum charge reference')
        require(case not in data,'Duplicate endpoint')
        data[case]=energy
    require(set(data)==set(expected),'Incomplete endpoint family')
    return data,dict(completed=manifest(directory/'completed.json'),producer=complete['producer'])


def host_difference(data,host,size,key):
    d=data[f'{host}-{size}-DDfar'];h=data[f'{host}-{size}-He']
    dc,hc=d['config'],h['config'];dg,hg=dc['geometry'],hc['geometry']
    require(dc['valence_electrons']==hc['valence_electrons'],'Host endpoints have unequal charge totals')
    require(dg['symbols'][:size]==hg['symbols'][:size] and dg['symbols'][size:]==['H','H']
            and hg['symbols'][size:]==['He'],'Wrong united-atom species')
    require(dg['cell_A']==hg['cell_A'] and dg['positions_A'][:size]==hg['positions_A'][:size]
            and dg['positions_A'][size]==hg['positions_A'][size],
            'Host or united-atom site changed between electronic endpoints')
    return d[key]-h[key]


def atomic_reference_offset(calibration,pd):
    require(calibration['schema']=='coherence-neutral-PAW-reference-calibration-v1',
            'Wrong PAW energy calibration role')
    references=calibration['atomic_reference_energy_eV']
    offset=2*references['H']-references['He']
    require(np.isclose(offset,calibration['reference_offset_2H_minus_He_eV'],atol=1e-10,rtol=0),
            'Inconsistent PAW species reference offset')
    for source in calibration['sources']:
        symbol=source['symbol'];energy=pd[f'vacuum-{symbol}-10A']
        require(source['gpw']==energy['gpw'],'Calibration belongs to another atomic wavefunction')
        require(np.isclose(energy['energy_extrapolated_zero_temperature_eV']+references[symbol],
                           calibration['neutral_atom_absolute_energy_eV'][symbol],atol=1e-10,rtol=0),
                'Atomic absolute/relative energy mismatch')
    return offset


def analyse(pd,rh,pd_source,rh_source,dos_rows,calibration):
    key='energy_extrapolated_zero_temperature_eV'
    vacuum={L:2*pd[f'vacuum-H-{L}A'][key]-pd[f'vacuum-He-{L}A'][key] for L in [8,10]}
    box_error=abs(vacuum[10]-vacuum[8]);rows=[]
    offset=atomic_reference_offset(calibration,pd)
    require(np.isclose(vacuum[10]+offset,calibration['neutral_vacuum_adiabatic_difference_eV'],
                       atol=1e-10,rtol=0),'Vacuum calibration differs from the consumed endpoints')
    for host,data in [('Pd_hydro_+0pct',pd),('Pd3Rh_L12',rh)]:
        delta={size:host_difference(data,host,size,key) for size in [4,8]}
        primary=delta[8]-vacuum[10]
        error=abs(delta[8]-delta[4])+box_error
        free=host_difference(data,host,8,'free_energy_eV')
        rows.append(dict(material_id=host,
            Uad_signed_neutral_metal_excess_eV=primary,
            Ue_adiabatic_neutral_reference_calibrated_eV=delta[8]+offset,
            host_DD_minus_He_extrapolated_eV={str(k):v for k,v in delta.items()},
            vacuum_2H_minus_He_extrapolated_eV={str(k):v for k,v in vacuum.items()},
            periodic_cell_and_sampling_change_eV=abs(delta[8]-delta[4]),
            vacuum_box_change_eV=box_error,combined_observed_control_change_eV=error,
            relative_control_change_of_metal_excess=error/max(abs(primary),1e-30),
            relative_control_change_of_Ue_ad=error/max(abs(delta[8]+offset),1e-30),
            eligible_15percent=bool(error<=.15*abs(delta[8]+offset) and abs(delta[8]+offset)>1e-10),
            host_entropy_extrapolation_difference_eV=delta[8]-free,
            control_scope='Observed4/8-metal periodic-cell plus k-grid-shape change and8/10A atom-box change; not a rigorous asymptotic image bound or all-electron physical error. Fixed ions, neutral Born-Oppenheimer D=H proxy.'))
    byid={r['material_id']:r for r in dos_rows}
    x=np.array([byid[r['material_id']]['dos']['by_metal_element']['Pd']['occupied_d_centre_relative_EF_eV'] for r in rows])
    y=np.array([r['Uad_signed_neutral_metal_excess_eV'] for r in rows])
    require(abs(x[1]-x[0])>1e-8,'Degenerate Pd-projected d-centres')
    central=float((y[1]-y[0])/(x[1]-x[0]))
    width=sum(r['combined_observed_control_change_eV'] for r in rows)/abs(x[1]-x[0])
    return dict(schema='coherence-neutral-adiabatic-v1',cases=rows,
        two_host_slope_eV_per_eV=central,
        observed_control_slope_span_eV_per_eV=[central-width,central+width],
        slope_eligible=all(r['eligible_15percent'] for r in rows),
        occupied_Pd_d_centres_relative_EF_eV=x.tolist(),
        sources=dict(Pd=pd_source,Pd3Rh=rh_source),
        atomic_reference_calibration=calibration, endpoint_settings_and_energies=dict(Pd=pd,Pd3Rh=rh),
        definition='Neutral Ue_ad=E(host+DDfar)-E(host+He) after restoring GPAW PAW atomic reference offsets. Signed metal excess=Ue_ad-[2E(H_atom_doublet)-E(He_atom_singlet)]. Same neutral electron total within each host pair.',
        interpretation='Neutral fixed-ion adiabatic total and metal excess are distinct. Neither is an ionic-beam limit or an all-electron material calculation; Pd/Rh core remains frozen. Do not add either to the independent RPA/Yukawa proxy. Raw GPAW cross-species differences are not physical without restoring reference energies.',
        slope_scope='Two-host association at the declared sites, conditional on computed occupied Pd-projected centres; interval propagates observed controls only, not statistical confidence or full physical uncertainty.')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pd',type=Path,required=True);parser.add_argument('--rh',type=Path,required=True)
    parser.add_argument('--dos-rows',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--reference-calibration',type=Path,required=True)
    args=parser.parse_args()
    pd,ps=read_unit(args.pd,'Pd_hydro_+0pct');rh,rs=read_unit(args.rh,'Pd3Rh_L12')
    result=analyse(pd,rh,ps,rs,json.loads(args.dos_rows.read_bytes()),
                   json.loads(args.reference_calibration.read_bytes()))
    result['sources']['qualified_DOS_rows']=manifest(args.dos_rows)
    result['sources']['consumer']=manifest(Path(__file__))
    result['sources']['reference_calibration']=manifest(args.reference_calibration)
    write_json(args.output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
