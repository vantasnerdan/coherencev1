"""Consume bound saved DFT primitives without admitting incomplete rate physics.

The pinned transport stores Eads as POSITIVE binding eV/D relative to half D2,
uses it in a Langmuir equilibrium, and front_flux is impingement/sticking based.
An electronic barrier cannot replace s_clean or that finite-temperature reservoir.
"""
import hashlib
import json
import math
from pathlib import Path

KB_EV_K = 8.617333e-5


def bound_json(binding, root=None):
    path = Path(binding['path'])
    if root is not None and not path.is_absolute():
        path = Path(root)/path
    raw = path.read_bytes()
    if len(raw) != binding['bytes'] or hashlib.sha256(raw).hexdigest() != binding['sha256']:
        raise ValueError('Saved primitive byte identity mismatch')
    return json.loads(raw)


def primitive(binding, root=None):
    value = bound_json(binding, root)
    energy = value['energy_extrapolated_eV']
    if not math.isfinite(energy):
        raise ValueError('Nonfinite electronic energy')
    geometry = value['geometry']
    if len(value['forces_eV_A']) != len(geometry['symbols']):
        raise ValueError('Primitive force/geometry dimension mismatch')
    return value


def adsorption(clean, adsorbed, atom, molecule, n):
    if n < 1 or not isinstance(n, int):
        raise ValueError('Explicit integer adsorbate stoichiometry required')
    if atom['geometry']['symbols'] != ['H'] or not atom['spin_polarized']:
        raise ValueError('Atomic reservoir must be an actually spin-polarized isolated H')
    if molecule['geometry']['symbols'] != ['H', 'H']:
        raise ValueError('Molecular reservoir must be actually computed H2')
    if clean['geometry']['cell_A'] != adsorbed['geometry']['cell_A']:
        raise ValueError('Unmatched slab cells')
    es = clean['energy_extrapolated_eV']
    ea = adsorbed['energy_extrapolated_eV']
    eh = atom['energy_extrapolated_eV']
    eh2 = molecule['energy_extrapolated_eV']
    return {'atomic_H_adsorption_eV_per_D': (ea-es-n*eh)/n,
            'half_H2_adsorption_eV_per_D': (ea-es-0.5*n*eh2)/n,
            'positive_electronic_binding_eV_per_D': -(ea-es-0.5*n*eh2)/n,
            'atomic_reservoir_spin_polarized': True}


def barrier_relation(initial_eV, saddle_eV, final_eV):
    forward = saddle_eV-initial_eV
    reverse = saddle_eV-final_eV
    if min(forward, reverse) < -1e-8:
        raise ValueError('Chosen saddle is below an endpoint')
    return {'desorption_barrier_eV_per_D2_event': forward,
            'reverse_dissociation_barrier_eV_per_D2_event': reverse,
            'reaction_energy_eV_per_D2_event': final_eV-initial_eV,
            'detailed_balance_energy_residual_eV': (forward-reverse)-(final_eV-initial_eV)}


def rate_range(replacement, temperature_K, pair_occupation):
    if not replacement['eligibility']['barrier']:
        raise ValueError('No converged MEP/TS barrier: rate refused')
    if temperature_K <= 0 or not 0 <= pair_occupation <= 1:
        raise ValueError('Invalid temperature or pair occupation')
    pref = replacement['prefactor']
    if pref['source_kind'] != 'externally_assumed_range':
        raise ValueError('This interface requires declared external prefactor range')
    barrier = replacement['barriers']['desorption_barrier_eV_per_D2_event']
    # Per explicitly defined reactive pair. Mapping to flux requires pair density.
    factor = pair_occupation*math.exp(-barrier/(KB_EV_K*temperature_K))
    return {'D2_events_per_pair_s': [x*factor for x in pref['range_s_inv']],
            'DFT_unique_rate': False, 'temperature_K': temperature_K,
            'ZPE_entropy': 'omitted; electronic-barrier-only conditional estimate'}


def admit_original_transport(replacement, finite_temperature_correction_eV_per_D):
    if not replacement['eligibility']['adsorption']:
        raise ValueError('Adsorption controls or matched reservoirs incomplete')
    if not replacement['eligibility']['original_transport_Eads']:
        raise ValueError('Coverage-average electronic binding is not a differential Langmuir free energy')
    if finite_temperature_correction_eV_per_D is None:
        raise ValueError('Half-D2-at-1Pa finite-temperature reservoir correction required')
    return {'Eads': replacement['adsorption']['positive_electronic_binding_eV_per_D']
            + finite_temperature_correction_eV_per_D,
            's_clean': 'unchanged: cannot infer sticking from energetics'}
