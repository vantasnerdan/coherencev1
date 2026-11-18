"""Portable saved ten-valence speed adapter; source identities are in ledger.json."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import elementary_charge, physical_constants
from g4_beam.screen.dynamic._upstream.electronic_speed import fermi_rms_velocity_au, speed_factor
from g4_beam.merit.api import ENERGIES

MASS_D_KG = physical_constants['deuteron mass'][0]
VELOCITY_AU_M_S = physical_constants['atomic unit of velocity'][0]
LAB_ENERGIES_KEV = [.5, 2., 5., 13.]


def load_models(path):
    raw = Path(path).read_bytes()
    source = json.loads(raw)
    models = {}
    controls = {}
    for case in source['cases']:
        name = case['case']
        material = case['material']
        expected = 10 if material == 'Pd' else 11 if material == 'PdD' else None
        expected_symbols = ['Pd'] if material == 'Pd' else ['Pd', 'H']
        if expected is None or case['geometry']['symbols'] != expected_symbols or case['expected_electrons'] != expected or abs(case['reconstructed_electrons']-expected) > 1e-5 or not case['no_new_SCF']:
            raise ValueError('Actual moment readout is not its declared ten-valence target')
        tensor = np.array(case['plasma_squared_au_real'])+1j*np.array(case['plasma_squared_au_imag'])
        velocity_au = fermi_rms_velocity_au(tensor, case['D_T_states_per_Hartree_cell'], case['volume_Bohr3'])
        if not np.isclose(velocity_au, case['band_RMS_velocity_au'], rtol=1e-12, atol=0):
            raise ValueError('Stored RMS velocity is inconsistent with its actual plasma/DOS moments')
        model = dict(material=material, velocity_au=velocity_au,
                     velocity_m_s=velocity_au*VELOCITY_AU_M_S,
                     U0_eV=case['U10_static_Yukawa_proxy_eV'], kmesh=case['kmesh'],
                     D_T_states_per_eV_cell=case['D_T_states_per_eV_cell'])
        if name.endswith('-control'):
            controls[name] = model
        else:
            if material in models:
                raise ValueError('Ambiguous duplicate material moment readouts')
            models[material] = model
    if set(models) != {'Pd', 'PdD'}:
        raise ValueError('Actual Pd and PdD moment readouts are both required')
    return models, controls, dict(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def screening_eV(model, laboratory_energy_keV):
    energy = np.asarray(laboratory_energy_keV, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0):
        raise ValueError('Finite nonnegative laboratory energy is required')
    speed = np.sqrt(2*energy*1000*elementary_charge/MASS_D_KG)
    return model['U0_eV']*speed_factor(speed/model['velocity_m_s'])


def comparison_rows(models, controls, lab_energies):
    rows = []
    for energy in lab_energies:
        pd = float(screening_eV(models['Pd'], energy))
        pdd = float(screening_eV(models['PdD'], energy))
        row = dict(E_lab_keV=float(energy), E_cm_keV=float(energy)/2,
                   U10_Pd_eV=pd, U10_PdD_eV=pdd, loading_PdD_minus_Pd_eV=pdd-pd,
                   fitted_loading_step_eV=[550., 790.],
                   signed_fraction_of_fitted_step=[(pdd-pd)/550., (pdd-pd)/790.])
        if 'Pd-k24-control' in controls:
            pd24 = float(screening_eV(controls['Pd-k24-control'], energy))
            row.update(U10_Pd_k24_velocity_control_eV=pd24,
                       loading_with_Pd_k24_velocity_eV=pdd-pd24,
                       Pd_velocity_grid_change_eV=pd24-pd)
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--moments', type=Path, required=True)
    args = parser.parse_args()
    models, controls, source = load_models(args.moments)
    result = dict(schema='pass4-tenvalence-Huke-speed-estimate-v1', source_moments=source,
                  models=models, velocity_controls=controls,
                  four_laboratory_energies=comparison_rows(models, controls, LAB_ENERGIES_KEV),
                  paper_grid=comparison_rows(models, controls, 2*ENERGIES),
                  constants=dict(deuteron_mass_kg=MASS_D_KG, atomic_velocity_m_s=VELOCITY_AU_M_S,
                                 elementary_charge_C=elementary_charge),
                  screening_scope='r->0 limit of accepted pass4 low-q ten-valence C/lambda Yukawa proxy; NOT full large-q U18/contact dielectric readout',
                  dynamic_scope='One isotropic Huke Eq10 factor using actual finite-width band-RMS velocity, not full frequency-dependent/local-field or real-time screening',
                  rate_mapping='Pd-minus-PdD contrast at the same incident lab energy, common PdD phenomenological calibration; never add absolute microscopic U to a fitted Ue',
                  uncertainty='Pd k12/k24 velocity sensitivity is measured; no PdD velocity-grid control exists here. Static proxies retain prior conditional numerical allowance; no global physical error bound.',
                  fitted_step_source='Historical pass4 +0.55 to +0.79keV loading step under the three charge readings; separate from literal Eq7 effective-fit definitions',
                  no_new_solid=True)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
