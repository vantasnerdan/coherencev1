"""Consume completed bulk cases by byte identity, actual geometry and numeric checks."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from read_fingerprint import bound_file, read_fingerprint
from screen_io import manifest, require
from screen_length import fit_length, nearest_periodic_distance, ratio_and_force


def read_case(directory, expected_material, dos_supplement=None):
    directory = Path(directory)
    complete = json.loads((directory/'completed.json').read_bytes())
    producer = json.loads(bound_file(directory, complete['producer']).read_bytes())
    config = producer['config']
    require(config['material_id'] == expected_material,
            'Completed source is not the requested material')
    require(complete['case'] == config['case'], 'Completion and producer case differ')
    bound_file(directory, producer['input'])
    require(json.loads((directory/'input.json').read_bytes()) == config,
            'Producer/input mismatch')
    for reference in producer['producer_source']:
        bound_file(directory, reference)
    ground = json.loads(bound_file(directory, complete['ground']).read_bytes())
    require(ground['case'] == config['case'], 'Wrong ground-state identity')
    require(ground['ground_state_settings'] == config['ground_state'],
            'Ground-state sampling differs from declared settings')
    dos = dict(ground['dos'])
    if dos_supplement is not None:
        require(dos_supplement['gpw'] == ground['gpw'], 'Supplement belongs to another wavefunction')
        require(dos_supplement['native_completed'] == producer['source_ground']['completed'],
                'Supplement belongs to another completed native source')
        dos = dict(dos_supplement['dos'])
        dos['saved_DOS_postprocess_provenance'] = dos_supplement['source_data']
    require('n_EF_tetrahedron_states_per_eV_per_cell' in dos and 'by_metal_element' in dos,
            'Full saved-state DOS analysis is required; no Gaussian or element-moment substitute')
    require(np.isfinite(dos['n_EF_tetrahedron_states_per_eV_per_cell'])
            and dos['n_EF_tetrahedron_states_per_eV_per_cell'] >= 0,
            'Invalid tetrahedron DOS')
    with np.load(bound_file(directory, dos['spectrum']), allow_pickle=False) as spectrum:
        energy = spectrum['energies_eV']
        projection = spectrum['occupied_metal_d_states_per_eV']
        require(energy.shape == projection.shape and np.all(np.diff(energy) > 0)
                and np.isfinite(projection).all(), 'Invalid projected spectrum')
        integral = float(simpson(projection, x=energy))
        centre = float(simpson(energy*projection, x=energy)/integral)
        require(abs(integral/dos['occupied_d_states_per_cell']-1) < 1e-4,
                'Saved occupied d count fails independent quadrature')
        require(abs(centre-dos['occupied_d_centre_relative_EF_eV']) < .001,
                'Saved d centre fails independent quadrature')
        for symbol, recorded in dos['by_metal_element'].items():
            element = spectrum['occupied_'+symbol+'_d_states_per_eV']
            count = float(simpson(element, x=energy))
            centre = float(simpson(energy*element, x=energy)/count)
            require(count > 0 and abs(count/recorded['occupied_d_states_per_cell']-1) < 1e-4,
                    'Element-resolved occupied d count fails independent quadrature')
            require(abs(centre-recorded['occupied_d_centre_relative_EF_eV']) < .001,
                    'Element-resolved d centre fails independent quadrature')
    result = dict(material_id=expected_material, native_case=config['case'],
                  geometry=ground['geometry'], ground_state_settings=config['ground_state'],
                  initial_geometry=config['geometry'], dos=dos,
                  source=dict(completed=manifest(directory/'completed.json'),
                              producer=complete['producer'], ground=complete['ground'],
                              setup_bytes=producer['setups']),
                  model=config['model'], kind=config['kind'])
    if config['kind'] == 'rpa':
        response = read_fingerprint(directory, config['case'])
        points = response['points']
        length = fit_length([p['Q_head_Ainv'] for p in points],
                            [complex(*p['epsilon_macro_LF']) for p in points])
        symbols = ground['geometry']['symbols']
        positions = np.asarray(ground['geometry']['positions_A'])
        cell = ground['geometry']['cell_A']
        metal = [a for a, s in enumerate(symbols) if s != 'H']
        hydrogen = [a for a, s in enumerate(symbols) if s == 'H']
        metal_distance = nearest_periodic_distance(positions[metal], cell)
        if hydrogen:
            dd = nearest_periodic_distance(positions[hydrogen], cell)
            spacing_kind = 'actual ordered D sublattice'
        else:
            dd = metal_distance
            spacing_kind = 'registered ideal octahedral fcc probe sublattice; no D present'
        result['screening'] = dict(**length,
                                  **ratio_and_force(length['lambda_A'], dd, metal_distance),
                                  DD_spacing_convention=spacing_kind,
                                  points=points, source_response=complete['response'],
                                  eligible_for_Tier1=length['q_window_pass_10percent'],
                                  uncertainty_scope='Pd k/band/cutoff/eta controls plus this q-window; transferred model tolerance, not all-electron physical confidence')
    else:
        require(config['kind'] == 'dos', 'Wrong completed scientific role')
        require('response' not in complete, 'DOS-only case pretends to contain RPA')
    return result
