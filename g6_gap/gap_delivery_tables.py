"""Publish P257 stopping/G bands and the requested 1keV LiT decision; NumPy only."""
import copy
import csv
import io
from pathlib import Path
import numpy as np
from screen_io import require,publish,write_json
from lid_return import return_curves
from stopping_threshold import _half_window,threshold_record
from substrate_framework.lithium_fuel_return import (
    lit_icru_electronic,lit_inherited_sqrt_electronic,lit_zbl_nuclear,
    stopping_per_formula,break_even_barns)


def write_csv(path,header,rows):
    text=io.StringIO(newline='')
    writer=csv.writer(text);writer.writerow(header);writer.writerows(rows)
    return publish(path,lambda stream:stream.write(text.getvalue().encode()))


def stopping_cases(stopping,analysis):
    e=np.asarray(stopping['energy_lab_keV'])
    require(np.array_equal(e,analysis['energy_lab_keV']),'Envelope energy support differs')
    central=np.asarray(stopping['stopping_eV_A'])
    error=np.asarray(analysis['combined_observed_L1_envelope_eV_A'])
    cases={}
    for name,values in [('q4',central),('low_drag',np.maximum(central-error,0)),('high_drag',central+error)]:
        record=copy.deepcopy(stopping)
        record['stopping_eV_A']=values.tolist()
        record['electronic_cross_section_per_formula']=stopping_per_formula(values,record['volume_A3']).tolist()
        cases[name]=record
    return cases


def stopping_table(cases,output):
    c=cases['q4'];e=np.asarray(c['energy_lab_keV']);keep=(e>=.1)&(e<=13)
    volume=c['volume_A3']
    columns=[e[keep],*[np.asarray(cases[k]['stopping_eV_A'])[keep] for k in ('low_drag','q4','high_drag')],
             lit_icru_electronic(e[keep])*10/volume,lit_inherited_sqrt_electronic(e[keep])*10/volume,
             lit_zbl_nuclear(e[keep])*10/volume]
    return write_csv(output/'stopping-band.csv',
        ['E_D_lab_keV','Se_low_eV_A','Se_q4_eV_A','Se_high_eV_A','Se_source_ICRU49_eV_A',
         'Se_inherited_sqrt_eV_A','Sn_ZBL_D_on_Li6T_eV_A'],zip(*columns))


def gain_table(results,output):
    central=results['q4']['curves'];model='source_ICRU49_native_replacement'
    rows=[];maxima={}
    for u in ('0.0','0.84','1.4'):
        e=np.asarray(central[model]['gains']['0.1']['energy_lab_keV']);keep=e>=.5
        baseline=np.asarray(central['source_ICRU49_baseline']['gains']['0.1']['G'][u])
        values={k:np.asarray(r['curves'][model]['gains']['0.1']['G'][u]) for k,r in results.items()}
        require(np.all(values['low_drag']+1e-16>=values['q4']) and np.all(values['q4']+1e-16>=values['high_drag']),
                'Gain envelope has reversed its physical stopping ordering')
        for i in np.flatnonzero(keep):
            rows.append([float(e[i]),float(u),float(baseline[i]),float(values['high_drag'][i]),
                         float(values['q4'][i]),float(values['low_drag'][i])])
        maxima[u]={k:dict(G=float(v[keep].max()),E_lab_keV=float(e[keep][np.argmax(v[keep])])) for k,v in values.items()}
    write_csv(output/'G-band.csv',['E_D_lab_keV','Ue_keV','G_ICRU_law_baseline',
        'G_high_drag_lower_return','G_q4','G_low_drag_upper_return'],rows)
    return maxima


def decision(cases):
    c=cases['q4'];energy=np.asarray(c['energy_lab_keV']);keep=(energy>=.1)&(energy<=13)
    nuclear=lit_zbl_nuclear(energy);one=int(np.flatnonzero(energy==1.)[0])
    windows={};at_one={}
    for name,r in cases.items():
        electronic=np.asarray(r['electronic_cross_section_per_formula'])
        window=_half_window(energy[keep],electronic[keep]/(2*nuclear[keep]))
        windows[name]=dict(electronic_equals_nuclear_bracket_keV=window['low_energy_window_upper_boundary_keV'],
                           nuclear_dominates_at_lower_edge=bool(electronic[keep][0]<=nuclear[keep][0]))
        total=electronic[one]+nuclear[one]
        at_one[name]=dict(electronic_cross_section_per_formula=float(electronic[one]),
            nuclear_ZBL_cross_section_per_formula=float(nuclear[one]),
            nuclear_fraction=float(nuclear[one]/total),break_even_barns=float(break_even_barns(total)))
    return dict(nuclear_dominated_windows=windows,at_1keV_lab=at_one,
        TiT1p7_reference_break_even_barns=200,TiT_reference_source='Axiom2220 comparator, not recalculated',
        scope='Electronic LiD BO stopping reused for LiT at the same volume; unchanged D-on-Li6 plus D-on-T ZBL. This is the declared model, not measured LiT loss or a physical gain verdict.')


def build_tables(stopping,analysis,bridge,output):
    output.mkdir(parents=True,exist_ok=False)
    cases=stopping_cases(stopping,analysis)
    results={}
    for name,record in cases.items():
        folder=output/name;folder.mkdir()
        results[name]=return_curves(record,folder,bridge)
    table=stopping_table(cases,output);maxima=gain_table(results,output)
    result=dict(schema='P257-gap-delivery-decision-v1',decision=decision(cases),
        half_reference_threshold=threshold_record(stopping),G_maxima=maxima,
        observed_settings_accuracy_pass=analysis['observed_15percent_absolute_floor_controls_pass'],
        observed_tolerance_ratio=analysis['max_required_domain_tolerance_ratio'],
        worst_energy_lab_keV=analysis['worst_energy_lab_keV'],stopping_table=table,
        literal_bridge=results['q4']['compiled_measured_Bragg_bridge'],
        Geant4_model_relative_overlap=results['q4']['native_join_relative_to_source_ICRU49'],
        envelope='Seven-control observed L1 envelope; nonnegative electronic lower bound, not statistical coverage or an infinite-basis bound.',
        returns='E_cut0.1keV retained. Eq7 bare/Ue0.84/1.4keV keeps the unshifted D-T S-factor and ZBL nuclear stopping. ICRU-model bridge above25.1keV is untuned; literal measured-compilation bridge fails30percent. All physical returns remain conditional.')
    write_json(output/'decision.json',result)
    return result
