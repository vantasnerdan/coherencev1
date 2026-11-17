"""Read actual original isotherm ranges; do not assert an electronic/free-energy bridge."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np


def main():
    source = Path('../../sources/coherencev1-e9c652c18e246aa5674bf0cb85dcbd1e8666261a/g4_beam/transport/transport.py').resolve()
    spec = importlib.util.spec_from_file_location('pinned_coherence_transport', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    x = np.array([1e-6, 1e-4, 0.001, 0.012, 0.1, 0.6, 0.7, 0.9, 0.95])
    rows = []
    for temperature in (300.0, 350.0, 450.0):
        iso = module.Isotherm(module.PD, temperature)
        mu = iso.mu(x)
        rows.append(dict(T_K=temperature, D_per_Pd=x.tolist(),
                         mu_relative_half_D2_1Pa_eV=mu.tolist(),
                         dilute_mu0_eV=float(iso.mu0()),
                         mu_minus_mu0_eV=(mu-iso.mu0()).tolist(),
                         baseline_trap_occupation_fraction=(1/(1+np.exp(-(mu-iso.mu0()+module.PD['Eb'])/iso.kT))).tolist()))
    temperatures = [dict(beam_energy_keV=e, cooling=c,
                         T_K=float(module.foil_temperature(e, module.PD, c)))
                    for e in (0.5, 6.5, 13.0) for c in ('electrolyte', 'dry')]
    result = dict(source=dict(bytes=source.stat().st_size,
                              sha256=hashlib.sha256(source.read_bytes()).hexdigest()),
                  original_model_reservoirs=rows,
                  original_foil_temperature_examples=temperatures,
                  limitation='mu-mu0 may be compared to electronic delta_mu ONLY after an explicit finite-T reference bridge, ionic free energies and hydride-background applicability. These ranges describe original model assumptions, not DFT-computed gas thermodynamics.')
    Path('baseline-reservoir-range.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
