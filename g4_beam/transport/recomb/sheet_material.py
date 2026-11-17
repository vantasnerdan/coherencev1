"""Build the physical sheet model only from all three settled material inputs.

Local lower-coverage free energies are re-evaluated from Native electronic
references and spectra. The accepted full-minus-local1ML shift is an explicit
approximation; configurational entropy belongs to SheetThermo, not these rows.
"""
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from sheet_thermo import fit_anchored
from thermo_harmonic import harmonic_parts


def read_bound(base, reference):
    raw = (base/reference['path']).read_bytes()
    if len(raw) != reference['bytes'] or sha256(raw).hexdigest() != reference['sha256']:
        raise ValueError('Material input byte identity changed: '+reference['path'])
    return json.loads(raw)


def temperature_rows(rows):
    result = {float(row['T_K']): row for row in rows}
    if len(result) != len(rows) or not result:
        raise ValueError('Distinct actual temperature rows required')
    return result


def build_material(ledger_path, expected_sha256):
    ledger_path = Path(ledger_path)
    raw = ledger_path.read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Material ledger identity changed')
    ledger = json.loads(raw)
    cases = {float(case['coverage_ML']): case for case in ledger['cases']}
    if len(cases) != 3 or len(ledger['cases']) != 3 or set(cases) != {.25, .5, 1.}:
        raise ValueError('Incomplete material ledger: actual0.25/0.5/1ML inputs all required')
    base = ledger_path.parent
    full = read_bound(base, cases[1.]['source'])
    full_rows = temperature_rows(full['rows'])
    declared_full = temperature_rows(cases[1.]['mean_binding_rows'])
    if set(full_rows) != set(declared_full):
        raise ValueError('Full1ML ledger temperatures differ from the actual source')
    gas_point = read_bound(base, ledger['gas_reference'])
    gas_energy = float(gas_point['energy_extrapolated_eV'])
    lower_inputs = {}
    for coverage in (.25, .5):
        case = cases[coverage]
        hessian = read_bound(base, case['source'])
        modes = np.asarray(hessian['spectrum']['signed_quanta_eV'], dtype=float)
        if modes.shape != (3,) or np.any(~np.isfinite(modes)) or np.any(modes <= 0):
            raise ValueError('Three positive actual local-D modes required')
        if hessian['center_mobile_fmax_eV_A'] > .05 or case['full_mobile_fmax_eV_A'] > .05:
            raise ValueError('Material geometry fails the declared0.05eV/A criterion')
        clean_ref = case['clean_reference']
        clean = read_bound(base, clean_ref)
        multiplier = clean_ref.get('multiplier', 1)
        if multiplier not in (1, 2, 4):
            raise ValueError('Declared clean-cell reference or qualified periodic folding required')
        clean_energy = multiplier*float(clean[clean_ref['energy_field']])
        electronic = -(float(hessian['center_energy_eV'])-clean_energy-.5*gas_energy)
        if abs(electronic-case['electronic_binding_eV_per_D']) > 1e-10:
            raise ValueError('Declared binding differs from actual matched electronic references')
        declared = temperature_rows(case['mean_binding_rows'])
        if set(declared) != set(full_rows):
            raise ValueError('Lower-coverage temperatures differ from the full1ML source')
        lower_inputs[coverage] = dict(modes=modes, electronic=electronic, declared=declared)
    models, rows = {}, []
    for temperature, source in full_rows.items():
        full_binding = float(source['full_mobile_thermochemistry']['Eads_eV_per_D'])
        local_full_binding = float(source['frozen_Pd_sensitivity']['Eads_eV_per_D'])
        if (abs(declared_full[temperature]['full_mobile_binding_eV_per_D']-full_binding) > 1e-12
                or abs(declared_full[temperature]['local_D_only_binding_eV_per_D']-local_full_binding) > 1e-12):
            raise ValueError('Full1ML ledger binding differs from its actual source')
        correction = full_binding-local_full_binding
        calibrated, details = [], []
        for coverage, data in lower_inputs.items():
            local = (data['electronic']-harmonic_parts(temperature, data['modes'])['F_eV']
                     +.5*source['gas_RRHO']['correction_G_eV'])
            if abs(local-data['declared'][temperature]['local_D_only_binding_eV_per_D']) > 1e-10:
                raise ValueError('Local mean free energy differs from the actual spectrum/gas reference')
            mean = local+correction
            calibrated.append((coverage, mean))
            details.append(dict(coverage_ML=coverage, local_mean_binding_eV_per_D=local,
                full_minus_local1ML_calibration_eV_per_D=correction,
                calibrated_mean_binding_eV_per_D=mean))
        model, fit = fit_anchored(temperature, full_binding, calibrated,
                                 source_id='threecoverage-ledger:'+expected_sha256)
        models[temperature] = model
        rows.append(dict(T_K=temperature, full_mean_binding_eV_per_D=full_binding,
            lower_coverage_calibration=details, fit=fit,
            differential_binding_at_0ML_eV_per_D=model.differential_binding(0.),
            differential_binding_at_1ML_eV_per_D=model.differential_binding(1.)))
    return models, dict(schema='coherence-physical-sheet-material-v1', rows=rows,
        ledger=dict(path=ledger_path.name, bytes=len(raw), sha256=expected_sha256),
        excess_free_energy_eV_per_site='-(B1+w)*theta+w*theta**2',
        differential_binding_eV_per_D='B1+w-2*w*theta',
        configurational_entropy='Ideal theta*ln(theta)+(1-theta)*ln(1-theta), added once',
        calibration_assumption='Measured full-minus-local1ML correction transferred to0.25/0.5ML',
        capture_probability_computed=False)
