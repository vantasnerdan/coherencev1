"""Reuse a completed ground state only with exact case/settings/byte agreement."""
import json
import os
from pathlib import Path

from screen_io import manifest, require


def admit_saved_ground(config, source, output):
    source = Path(source).resolve()
    completed = json.loads((source/'completed.json').read_bytes())
    producer = json.loads((source/'producer.json').read_bytes())
    ground = json.loads((source/'ground.json').read_bytes())
    require(manifest(source/'ground.json') == completed['ground'],
            'Source ground record does not match its completion receipt')
    require(manifest(source/'producer.json') == completed['producer'],
            'Source producer does not match its completion receipt')
    for field, value in [('case', ground['case']), ('model', ground['model']),
                         ('geometry', ground['geometry']),
                         ('ground_state', ground['ground_state_settings']),
                         ('dos', ground['dos']['settings'])]:
        require(config[field] == value == producer['config'][field],
                f'Reused ground has a different {field}')
    references = [ground['gpw'], ground['dos']['spectrum'], completed['ground']]
    for reference in references:
        name = reference['path']
        require(Path(name).name == name, 'Source artifact is not a local filename')
        require(manifest(source/name) == reference,
                f'Reused ground artifact binding failed: {name}')
    for reference in references:
        os.link(source/reference['path'], output/reference['path'])
    return dict(source_case_directory=source.name,
                completed=manifest(source/'completed.json'),
                producer=manifest(source/'producer.json'),
                input=manifest(source/'input.json'),
                ground=completed['ground'], gpw=ground['gpw'])
