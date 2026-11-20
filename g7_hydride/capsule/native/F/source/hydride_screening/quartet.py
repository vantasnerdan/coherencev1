"""Actual neutral/bare probe quartets on one matched frozen relaxed host."""
import gc
import math
from pathlib import Path
import time

import numpy as np
from ase import Atom
from ase.units import Bohr
from substrate_framework.collision_bare_reference import PeriodicNuclearGreen
from substrate_framework.collision_observables import COULOMB_EV_A

from geometry import MASS_D
from native_state import require, scf, write_json


def probe_stem(role, charge, radius):
    return role+'-probe-q'+str(charge)+'-R'+format(radius, '.3f').replace('.', 'p')


def quartet(atoms, geometry, spec, radii, output, reuse=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    reuse = reuse or {}
    started = time.process_time()
    target_index = geometry['target_index']
    target = np.asarray(geometry['path']['target_A'])
    direction = np.asarray(geometry['path']['direction'])
    require(abs(float(np.linalg.norm(direction))-1) < 1e-12, 'Nonunit projectile direction')
    require(atoms[target_index].symbol == 'H', 'Target is not D/H')
    original = atoms.copy()
    original.calc = None
    states = {'occupied-host': scf(original, spec, 0, output, 'occupied-host',
                                   reuse=reuse.get('occupied-host'))}
    original.calc = None
    vacant = original.copy()
    del vacant[target_index]
    states['vacant-host'] = scf(vacant, spec, 0, output, 'vacant-host')
    vacant.calc = None
    gc.collect()
    require(abs(states['occupied-host']['electrons']-states['vacant-host']['electrons']-1) < 1e-7,
            'Target removal must remove one nucleus and one electron')
    rows = {}
    green = PeriodicNuclearGreen(original.cell.array)
    width = np.sqrt(5.)/Bohr
    for charge in (0, 1):
        key = 'q'+str(charge)
        states[key] = {}
        for role, host in [('occupied', original), ('vacant', vacant)]:
            probe = host.copy()
            probe += Atom('H', position=target-float(radii[0])*direction, mass=MASS_D)
            records = []
            for radius in radii:
                radius = float(radius)
                require(radius > 0, 'Nonpositive pair separation')
                probe.positions[-1] = target-radius*direction
                stem = probe_stem(role, charge, radius)
                record = scf(probe, spec, charge, output, stem,
                             extra=dict(R_A=radius, projectile_A=probe.positions[-1].tolist()),
                             reuse=reuse.get(stem))
                require(abs(record['electrons']-states[role+'-host']['electrons']-(1-charge)) < 1e-7,
                        'Wrong probe electron addition for declared charge')
                records.append(record)
            states[key][role+'-probes'] = records
            probe.calc = None
            del probe
            gc.collect()
        rows[key] = []
        for occupied, unoccupied in zip(states[key]['occupied-probes'], states[key]['vacant-probes']):
            radius = occupied['R_A']
            require(radius == unoccupied['R_A'], 'Unmatched quartet radii')
            raw = math.fsum([occupied['energy']['internal_energy_eV'],
                             -unoccupied['energy']['internal_energy_eV'],
                             -states['occupied-host']['energy']['internal_energy_eV'],
                             states['vacant-host']['energy']['internal_energy_eV']])
            correction = green.point_minus_gaussian_eV(-radius*direction, width)
            total = raw+correction
            rows[key].append(dict(R_A=radius, BO_quartet_soft_eV=raw,
                                 point_minus_gaussian_eV=correction,
                                 V_pair_total_eV=total, lowering_eV=COULOMB_EV_A/radius-total))
    result = dict(schema='hydride-quartet-v1', case=spec['case'], specification=spec,
                  geometry=geometry, states=states, rows=rows,
                  process_cpu_s=time.process_time()-started,
                  energy_zero='Matched frozen background: Eprobe,o-Eprobe,v-Ehost,o+Ehost,v; neutral target-only removal.',
                  charge_treatment='Periodic uniform compensating background for q1; same-cell matched probe q makes a q^2 constant cancel. Residual images require F cell control.',
                  scope='Finite-R static nonspin PBE/FD0.1 family; no experimental U0, prepared-state rate or exact isolated-limit claim.')
    write_json(output/'quartet.json', result)
    return result
