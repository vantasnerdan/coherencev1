"""One late100 local-spectrum row; never updates any completed trajectory or clock."""
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'pd_x1_quantum'/'engine'))
from anharmonic import spectrum


def main():
    start=time.perf_counter()
    source=json.loads((ROOT/'updated/ni_oscillator_924_input.json').read_text())
    old=json.loads((ROOT/source['frozen_quantum_result']).read_text())
    material=old['material_inputs']
    A=source['late_K4_100_eV_A4'];D=material['K4_111_eV_A4'];B=(3*D-A)/2
    actual=spectrum(material['K2_eV_A2'],A,B,material['mass_eVc2'],10,material['J_sum_eV_A2'])
    keys=('zero_to_one_meV','parallel_one_to_two_lines','perpendicular_one_to_two_lines',
        'parity_basis_per_axis','largest_block')
    result={'source':source,'temperature_scope':'T0 local eigenenergy spectra; not a new bath/dynamics case',
        'scope':'Local excitation energies only. All static, pumped-clock and trajectory results retain their original frozen-input epoch.',
        'frozen':{key:old['spectrum'][key] for key in keys},
        'late100_old111':{key:actual[key] for key in keys},
        'late100_old111_K4_eV_A4':{'100':A,'111':D},
        'zero_to_one_shift_meV':actual['zero_to_one_meV']-old['spectrum']['zero_to_one_meV'],
        'zero_to_one_relative_shift':actual['zero_to_one_meV']/old['spectrum']['zero_to_one_meV']-1,
        'elapsed_s':time.perf_counter()-start}
    path=ROOT/'updated/ni_epoch/oscillator_924.json'
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
