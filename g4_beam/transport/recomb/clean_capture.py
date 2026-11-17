"""Conditional detailed-balance factor from a sampled electronic route ridge."""
import math
from coverage_inputs import make_model
from elbow_path_results import bound_json


def capture_from_elbow(reference, *, source, rows, assumed_capture):
    readout = bound_json(reference)
    if (readout.get('schema') != 'coherence-connected-elbow-readout-v1'
            or not readout.get('complete_route')):
        raise ValueError('A complete authenticated connected-route readout is required')
    energies = [float(point['energy_eV']) for point in readout['curve']]
    baseline = float(readout['clean_plus_gas_energy_eV'])
    sampled = max(energies)-baseline
    if (not all(math.isfinite(value) for value in energies+[baseline])
            or not math.isclose(sampled, readout['sampled_peak_excess_eV'], rel_tol=0, abs_tol=1e-12)):
        raise ValueError('Reported ridge disagrees with its actual sampled energy curve')
    barrier = max(0., sampled)
    factors, by_temperature = [], {}
    for row in rows:
        model = make_model(source, row, s_effective=1.)
        factor = math.exp(-barrier/model.kT)
        by_temperature[row['T_K']] = assumed_capture*factor
        factors.append(dict(temperature_K=row['T_K'], kT_eV=model.kT,
            relative_rate_factor=factor, conditional_capture=assumed_capture*factor))
    metadata = dict(source_readout=reference, sampled_ridge_excess_eV=sampled,
        nonnegative_electronic_ridge_eV=barrier, assumed_base_capture=assumed_capture,
        factors=factors, same_factor_in_adsorption_and_desorption=True,
        capture_probability_computed=False,
        interpretation='Conditional Boltzmann factor of the sampled restricted-path electronic ridge; not a saddle free energy or a computed capture probability',
        force_interpolated_energy_admitted_as_DFT=False)
    return lambda temperature: by_temperature[temperature], metadata
