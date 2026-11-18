"""Portable saved Pd18 XML/atomic parser; source identities are in ledger.json."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


def binding(path):
    data = Path(path).read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def read_identity(xml):
    root = ET.parse(xml).getroot()
    atom = root.find('atom').attrib
    xc = root.find('xc_functional').attrib
    generator = root.find('generator')
    if atom['symbol'] != 'Pd' or float(atom['Z']) != 46 or float(atom['valence']) != 18 or float(atom['core']) != 28:
        raise ValueError('Expected actual Pd Z46/core28/valence18')
    if xc.get('name') != 'PBE' or generator.get('type') != 'scalar-relativistic':
        raise ValueError('Expected scalar-relativistic PBE dataset')
    states = [dict(s.attrib) for s in root.find('valence_states')]
    occupied = {}
    for state in states:
        occupation = float(state.get('f', 0.))
        if occupation > 0:
            if 'n' not in state:
                raise ValueError('An occupied state must identify its principal quantum number')
            occupied[int(state['n']), int(state['l'])] = occupation
    if occupied.get((4, 0)) != 2 or occupied.get((4, 1)) != 6 or abs(sum(occupied.values())-18) > 1e-10:
        raise ValueError('Expected occupied4s2/4p6 and total18 reference valence electrons')
    return dict(atom=atom, xc_functional=xc, generator_attributes=generator.attrib,
                generator_parameters=generator.text.strip(), valence_states=states)


def phase_diagnostics(path):
    with np.load(path, allow_pickle=False) as curves:
        energies = curves['energy_Hartree']
        if energies.ndim != 1 or not np.all(np.isfinite(energies)) or np.any(np.diff(energies) <= 0):
            raise ValueError('Expected finite ordered atomic energy grid')
        rows = []
        for ell in range(4):
            ae, ps = curves[f'ae_l{ell}'], curves[f'pseudo_l{ell}']
            if ae.shape != energies.shape or ps.shape != energies.shape or not np.all(np.isfinite(ae+ps)):
                raise ValueError('Incomplete atomic logarithmic derivative curves')
            difference = ae-ps
            offset = np.rint(difference[0])
            aligned = difference-offset
            wrapped = (difference+.5) % 1.-.5
            rows.append(dict(
                angular_momentum=ell, alignment_integer=float(offset),
                max_abs_aligned_arctan_logderivative_over_pi=float(np.max(np.abs(aligned))),
                max_abs_wrapped_arctan_logderivative_over_pi=float(np.max(np.abs(wrapped))),
                endpoint_winding_difference=float(np.rint(aligned[-1])),
                max_difference_energy_Hartree=float(energies[np.argmax(np.abs(aligned))])))
    return dict(radius_Bohr=3., source_API_quantity='unwrapped arctan(u_prime/u)/pi, not a free-space scattering phase',
                energy_range_Hartree=[float(energies[0]), float(energies[-1])],
                energy_step_Hartree=float(energies[1]-energies[0]), channels=rows,
                curves=binding(path),
                sampled_proxy_gate_pass=all(row['max_abs_aligned_arctan_logderivative_over_pi'] <= .02 and row['endpoint_winding_difference'] == 0 for row in rows),
                limitation='Sampled atomic diagnostics do not prove absence of all scalar-relativistic ghost states.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', required=True, type=Path)
    args = parser.parse_args()
    directory = args.directory
    result = dict(schema='pd18-saved-atomic-readout-v1',
                  identity=read_identity(directory/'Pd.18.PBE'), XML=binding(directory/'Pd.18.PBE'),
                  logarithmic_derivative_diagnostics=phase_diagnostics(directory/'logarithmic-derivatives.npz'),
                  generation_log=binding(directory/'generation.txt'),
                  qualification='Identity and sampled diagnostics only; solid numerical controls are reported separately. Atomic execution records remain private.')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
