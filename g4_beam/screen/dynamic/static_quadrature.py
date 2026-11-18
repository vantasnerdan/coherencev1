"""Assemble complete declared midpoint quadratures from actual saved banks."""
from pathlib import Path

import numpy as np
from g4_beam.screen.dynamic._upstream.microscopic_screening import q_residual_eV, reference_eV
from g4_beam.screen.dynamic._upstream.screening_quadrature import midpoint_classes, validate_site_group
from g4_beam.screen.dynamic.largeq_bank import binding, load_ground, require, response_points
from g4_beam.screen.dynamic.largeq_results import COULOMB_EV_A, RADII_A


REFERENCES = ((2., 4.), (1., 2.), (2., 2.), (3., 4.))


class SavedBanks:
    """Cache small validated response matrices, never repeated large GPW reads."""
    def __init__(self, base):
        self.base = Path(base)
        self.cache = {}

    def get(self, locator, material, valence):
        path = (self.base / locator).resolve()
        key = (path, material, valence)
        if key not in self.cache:
            bank = load_ground(path, material, valence)
            points = {point['index']: point for point in response_points(path, bank)}
            self.cache[key] = bank, points
        return self.cache[key]


def _point_class(q, n, members):
    numerator = np.rint(np.asarray(q) * (2 * n)).astype(int)
    require(np.max(np.abs(np.asarray(q) * (2 * n) - numerator)) < 1e-10,
            'Saved q is not on the declared midpoint grid')
    numerator = (numerator + n) % (2 * n) - n
    key = tuple(numerator)
    require(key in members, 'Saved q is outside the declared midpoint grid')
    return members[key]


def _complete_points(declaration, case, banks):
    n = declaration['n']
    classes = midpoint_classes(n)
    members = {tuple(member): index for index, row in enumerate(classes)
               for member in row['members_numerators']}
    selected, seen, anchor = [], set(), None
    for reference in declaration['points']:
        bank, points = banks.get(reference['bank'], case['material'], case['pd_valence'])
        index = reference['index']
        require(index in points, 'Declared response point is not producer-published')
        point = points[index]
        settings = point['settings']
        require(all(settings[name] == expected for name, expected in case['response'].items()),
                'Saved response does not match the declared cutoff/bands/broadening')
        require(bank['config']['ground_state'] == case['ground_state'],
                'Saved ground parameters differ from the declared numerical model')
        validate_site_group(bank['cell'], bank['ground']['geometry']['positions_A'],
                            bank['ground']['geometry']['symbols'], bank['source_A'],
                            bank['config']['ground_state']['kmesh'])
        if anchor is None:
            anchor = bank
        else:
            require(np.allclose(bank['cell'], anchor['cell'], rtol=0, atol=1e-9) and
                    np.allclose(bank['source_A'], anchor['source_A'], rtol=0, atol=1e-8),
                    'One quadrature mixes inequivalent cells or source sites')
        group = _point_class(settings['q_fractional'], n, members)
        require(group not in seen, 'Duplicate midpoint symmetry class; weights would be wrong')
        seen.add(group)
        selected.append((bank, point, classes[group]['weight'], reference))
    require(seen == set(range(len(classes))),
            'Incomplete midpoint quadrature cannot produce a total screening potential')
    require(abs(sum(row[2] for row in selected) - 1) < 1e-12,
            'Normalized quadrature weight invariant failed')
    return selected


def evaluate_case(case, banks, continuum_references):
    quadratures = []
    anchor = None
    for declaration in case['quadratures']:
        selected = _complete_points(declaration, case, banks)
        observed = selected[0][0]
        if anchor is None:
            anchor = observed
        else:
            require(np.allclose(observed['cell'], anchor['cell'], rtol=0, atol=1e-9) and
                    np.allclose(observed['source_A'], anchor['source_A'], rtol=0, atol=1e-8) and
                    np.allclose(observed['ground']['geometry']['positions_A'],
                                anchor['ground']['geometry']['positions_A'], rtol=0, atol=1e-8),
                    'One model mixes different physical geometries across q grids')
        profiles = []
        for (kappa, width), continuum in zip(REFERENCES, continuum_references):
            residual = np.zeros(len(RADII_A))
            for bank, point, weight, _ in selected:
                residual += q_residual_eV(
                    point['Q_Ainv'], point['inverse'], bank['source_A'], RADII_A,
                    bank['volume_A3'], weight, kappa, width, COULOMB_EV_A)
            total = residual + continuum  # Add the continuum reference ONCE.
            profiles.append(dict(kappa_Ainv=kappa, gaussian_width_Ainv=width,
                                 residual_eV=residual.tolist(), U_eV=total.tolist()))
        primary = np.asarray(profiles[0]['U_eV'])
        sensitivity = np.max(np.abs(np.asarray([p['U_eV'] for p in profiles]) - primary), axis=0)
        quadratures.append(dict(n=declaration['n'], full_midpoints=declaration['n']**3,
            computed_classes=len(selected), total_weight=sum(row[2] for row in selected),
            primary_U_eV=primary.tolist(), reference_profiles=profiles,
            reference_parameter_max_difference_eV=sensitivity.tolist(),
            actual_points=[dict(reference=row[3], weight=row[2],
                settings=row[1]['settings'], checks=row[1]['checks']) for row in selected]))
    return dict(id=case['id'], material=case['material'], pd_valence=case['pd_valence'],
                ground_state=case['ground_state'], response=case['response'],
                geometry=anchor['ground']['geometry'], source_A=anchor['source_A'].tolist(),
                cell_volume_A3=anchor['volume_A3'],
                quadratures=quadratures)


def evaluate_registry(registry, base):
    banks = SavedBanks(base)
    continuum = [reference_eV(RADII_A, kappa, width, COULOMB_EV_A)
                 for kappa, width in REFERENCES]
    ids = [case['id'] for case in registry['cases']]
    require(len(ids) == len(set(ids)) and ids, 'Declared model identities must be unique/nonempty')
    cases = [evaluate_case(case, banks, continuum) for case in registry['cases']]
    source_banks = []
    for (path, material, valence), (bank, points) in banks.cache.items():
        source_banks.append(dict(locator=str(path.relative_to(Path(base).resolve())),
            material=material, pd_valence=valence, files=bank['bindings'] +
            [record for point in points.values() for record in point['bindings']]))
    return dict(schema='saved-microscopic-static-profiles-v1',
        purpose=registry['purpose'], radii_A=RADII_A.tolist(), coulomb_eV_A=COULOMB_EV_A,
        convention='Spherical source-site test-charge interaction; full local fields; no self-energy half.',
        primary_reference=dict(kappa_Ainv=2., gaussian_width_Ainv=4.),
        cases=cases, source_banks=source_banks,
        qualification='Finite declared q/G/band/k/basis profiles, not an automatic physical error bound.',
        runtime_bindings=[binding(Path(__file__)),
            binding(Path(__file__).parent / '_upstream/microscopic_screening.py'),
            binding(Path(__file__).parent / '_upstream/screening_quadrature.py')])
