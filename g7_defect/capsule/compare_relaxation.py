"""Index-preserving positions relative to the registered vacancy; no setup-equivalence claim."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from ase.geometry import find_mic


def binding(path):
    raw = path.read_bytes()
    return dict(path=path.name, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def aligned_change(paw_positions, paw_centre, ae_positions, ae_centre, cell):
    p, _ = find_mic(np.asarray(paw_positions)-paw_centre, cell, True)
    a, _ = find_mic(np.asarray(ae_positions)-ae_centre, cell, True)
    delta, _ = find_mic(p-a, cell, True)
    return delta


def compare(paw_path, ae_path):
    paw, ae = (json.loads(path.read_bytes()) for path in (paw_path, ae_path))
    require(paw['passed'] and paw['raw_fmax_eV_A'] <= .1, 'PAW geometry is not admitted at Axiom2518 threshold')
    require(ae['passed'] is False, 'This comparison must retain the unconverged AE trace status')
    require(paw['symbols'] == ae['symbols'] and paw['target_index'] == ae['target_index'], 'Different physical inventories or target labels')
    cell = np.asarray(paw['cell_A'])
    require(np.allclose(cell, ae['cell_A'], rtol=0, atol=1e-12), 'Different cells cannot define this one-cell comparison')
    require(np.allclose(paw['initial']['vacancy_initial_A'], ae['initial']['vacancy_initial_A'], rtol=0, atol=1e-12), 'Different original vacancy registration')
    p, a = np.asarray(paw['positions_A']), np.asarray(ae['positions_A'])
    pc, ac = (np.asarray(value['vacancy']['registered_centre_A']) for value in (paw, ae))
    change = aligned_change(p, pc, a, ac, cell)
    require(np.isfinite(change).all(), 'Nonfinite position difference')
    pd = np.asarray(paw['symbols']) == 'Pd'; target = paw['target_index']
    # A global translation must not masquerade as a setup-dependent displacement.
    shift = np.array([.31, -.27, .19])
    shifted = aligned_change(p+shift, pc+shift, a, ac, cell)
    registration_error = float(np.max(np.abs(shifted-change)))
    require(registration_error < 1e-11, 'Vacancy registration depends on the global origin')
    require(change.shape == (len(paw['symbols']), 3), 'Lost atom in the position comparison')
    return dict(schema='pd-defect-position-provenance-v1', sources=dict(paw=binding(paw_path), ae_trace=binding(ae_path)),
                paw_admission_fmax_eV_A=paw['raw_fmax_eV_A'], paw_force_gate_eV_A=.1,
                ae_trace_fmax_eV_A=ae['raw_fmax_eV_A'], ae_trace_converged=False,
                cell_A=cell.tolist(), registration='Same seeded atom labels; each configuration relative to its own registered vacancy centre; periodic minimum images via ASE.',
                displacement_PAW_minus_AE_A=change.tolist(),
                Pd_rms_displacement_A=float(np.sqrt(np.mean(np.sum(change[pd]**2, axis=1)))),
                Pd_max_displacement_A=float(np.linalg.norm(change[pd], axis=1).max()),
                target_displacement_A=change[target].tolist(), target_displacement_norm_A=float(np.linalg.norm(change[target])),
                target_offset_PAW_A=paw['vacancy']['target_offset_A'], target_offset_AE_trace_A=ae['vacancy']['target_offset_A'],
                global_translation_sensitivity_A=registration_error,
                scope='Provenance against an unconverged AE optimization trace, not two converged setup geometries, a relaxation-energy comparison, or a physical screening error bound.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('paw', 'ae_trace', 'output'):
        parser.add_argument('--'+name.replace('_', '-'), type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.paw, args.ae_trace)
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    with args.output.open('x') as stream:
        stream.write(payload)
    print(payload, end='')
