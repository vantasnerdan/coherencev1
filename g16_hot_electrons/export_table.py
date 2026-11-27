"""Write G15-consumable temperature, excitation, and stiffness rows."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


SCHEMA = 'g16-pdd-thermal-material-v1'


def flatten(readout, sha):
    if readout['schema'] != SCHEMA or readout['reference_temperature_K'] != 300:
        raise ValueError('Not the G16 response on its 300 K reference')
    out = []
    for r in readout['rows']:
        for k in (12, 16):
            v = r[f'mesh{k}']
            harmonic = v['harmonic']
            contrast = v['contrast']
            zero = v['zero']
            fit100 = v['quartic'].get('100', {}).get('fit')
            fit111 = v['quartic'].get('111', {}).get('fit')
            step = contrast.get('delta_K_step_sensitivity_eV_A2')
            independent = r['mesh_check']
            out.append(dict(temperature_K=r['temperature_K'], kmesh_axis_count=k, reference_temperature_K=300,
                            volume_A3_per_D_site=zero['volume_A3'],
                            electronic_excitation_eV_per_D_site=contrast['electronic_excitation_eV_per_D_site'],
                            electronic_free_energy_change_eV_per_D_site=contrast['electronic_free_energy_change_eV_per_D_site'],
                            K_Gamma_100_eV_A2=harmonic['K_force_eV_A2'],
                            K_change_eV_A2=contrast['K_change_eV_A2'], fractional_K_change=contrast['fractional_K_change'],
                            fractional_K_secant_per_eV=contrast['fractional_K_secant_per_eV'],
                            delta_K_mesh_sensitivity_eV_A2=independent['delta_K_mesh_sensitivity_eV_A2'],
                            delta_K_force_energy_difference_eV_A2=independent['delta_K_force_energy_difference_eV_A2'],
                            delta_K_step_sensitivity_eV_A2=step,
                            delta_K_sign_larger_than_sensitivities=(independent['delta_K_sign_larger_than_sensitivities']
                                                                    and (step is None or abs(contrast['K_change_eV_A2']) > step)),
                            D_centered_force_norm_eV_A=math.hypot(*zero['forces_eV_A'][1]),
                            D_centered_force_change_norm_eV_A=math.hypot(*contrast['centered_D_force_change_eV_A']),
                            K4_100_eV_A4=None if fit100 is None else fit100['K4_force_fit_eV_A4'],
                            K4_111_eV_A4=None if fit111 is None else fit111['K4_force_fit_eV_A4'],
                            pressure_GPa=zero['pressure_GPa'], pressure_change_GPa=contrast['pressure_change_GPa'],
                            readout_sha256=sha))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('readout', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    raw = args.readout.read_bytes()
    rows = flatten(json.loads(raw), hashlib.sha256(raw).hexdigest())
    with args.output.open('x', newline='') as f:
        csv_out = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n')
        csv_out.writeheader()
        csv_out.writerows(rows)


if __name__ == '__main__':
    main()
