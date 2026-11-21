"""Source-bound, falsifiable contact mapping proof (stdlib and NumPy only).

For a normalized exp(-alpha*r*r) nuclear Gaussian the pair width is
sqrt(alpha/2).  At coincidence its cross energy is -2*M_H.  Thus the
same-alpha Z=2 replacement changes only the self reference: 2*M_H-4*M_H.
For the frozen finite-R mapping, erfc(wR)/R = 1/R-2*w/sqrt(pi)+O(R**2).
Every nonzero lattice term is regular.  Consequently
lim(Vtotal-C/R) = Q_He-2*M_H*Hartree-2*C*w/sqrt(pi)+regular_images.
The unit constants are independent; their small residual is NOT set to zero.
This program proves the mapping, not SCF convergence or an isolated limit.
"""
import argparse
import ast
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

C = 14.399645468667815
BOHR_A = .5291772105638411
WIDTH = math.sqrt(5.) / BOHR_A
MAPPING_SHA256 = '05ccba36601567b71a12af42c309cd978e1637b454454740c0f2909725598f09'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_binding(path):
    data = Path(path).read_bytes()
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def _expression(text):
    return ast.dump(ast.parse(text, mode='eval').body, include_attributes=False)


def verify_operator_source(path):
    """Check executable AST expressions, not a keyword/token proof."""
    tree = ast.parse(Path(path).read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'HydrogenAllElectronSetup')
    init = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == '__init__')
    require([ast.literal_eval(v) for v in init.args.defaults] == [10., 300.], 'Changed Hae alpha defaults')
    assignments = {}
    for node in init.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name) and target.value.id == 'self':
                    assignments[target.attr] = ast.dump(node.value, include_attributes=False)
    for key, expression in dict(E='0.0', Z='1', Nv='1', Nc='0', ni='0',
                                M='-(alpha1 / 2 / np.pi)**0.5',
                                Delta0='-1 / (4 * np.pi)**0.5',
                                ghat_l='[Spline.from_data(0, rc, 4 * alpha1**1.5 / np.pi**0.5 * np.exp(-alpha1 * r2_g))]',
                                vbar='Spline.from_data(0, rc, v_g)').items():
        require(assignments.get(key) == _expression(expression), 'Changed Hae source expression: ' + key)
    for statement in ('v_g = erf(alpha1**0.5 * r_g) - erf(alpha2**0.5 * r_g)',
                      'v_g[1:] *= (4 * np.pi)**0.5 / r_g[1:]',
                      'v_g[0] = 4 * (alpha1**0.5 - alpha2**0.5)'):
        wanted = ast.dump(ast.parse(statement).body[0], include_attributes=False)
        require(any(ast.dump(n, include_attributes=False) == wanted for n in init.body),
                'Changed Hae local-potential source expression')
    # Full source binding is retained in addition to the checked expressions.
    return source_binding(path)


def verify_ase_units(path, actual_hartree, codata_version):
    """Re-derive the actual ASE unit independently from its bound CODATA source."""
    tree = ast.parse(Path(path).read_text())
    assignments = {t.id: n.value for n in tree.body if isinstance(n, ast.Assign)
                   for t in n.targets if isinstance(t, ast.Name)}
    require(codata_version == ast.literal_eval(assignments['__codata_version__']) == '2014',
            'Unexpected actual ASE CODATA epoch')
    table = assignments['CODATA']
    require(isinstance(table, ast.Dict), 'Changed ASE CODATA source structure')
    entries = dict(zip((ast.literal_eval(k) for k in table.keys), table.values))
    selected = entries['2014']
    fields = dict(zip((ast.literal_eval(k) for k in selected.keys), selected.values))
    constants = dict(_c=299792458., _hplanck=6.626070040e-34, _e=1.6021766208e-19, _me=9.10938356e-31)
    for key, expected in constants.items():
        require(ast.literal_eval(fields[key]) == expected, 'Changed source ASE constant ' + key)
    require(ast.dump(fields['_mu0'], include_attributes=False) == _expression('4.0e-7*pi'),
            'Changed source ASE permeability')
    for text in ("u['_eps0'] = 1 / u['_mu0'] / u['_c']**2",
                 "u['_hbar'] = u['_hplanck'] / (2*pi)",
                 "u['Hartree'] = u['_me'] * u['_e']**3 / 16 / pi**2 / u['_eps0']**2 / u['_hbar']**2"):
        wanted = ast.dump(ast.parse(text).body[0], include_attributes=False)
        require(any(ast.dump(n, include_attributes=False) == wanted for n in ast.walk(tree)),
                'Changed ASE Hartree derivation')
    eps0 = 1 / (4e-7*math.pi) / constants['_c']**2
    hbar = constants['_hplanck'] / (2*math.pi)
    expected = constants['_me']*constants['_e']**3/16/math.pi**2/eps0**2/hbar**2
    require(math.isfinite(actual_hartree) and abs(actual_hartree-expected) <= 1e-12,
            'Native ASE Hartree is not the independently bound source unit')
    return dict(source=source_binding(path), CODATA_version=codata_version,
                independently_derived_Hartree_eV=expected, actual_Hartree_eV=actual_hartree,
                residual_eV=actual_hartree-expected)


def regular_mapping(cell_A):
    cell = np.asarray(cell_A, dtype=float)
    require(cell.shape == (3, 3) and np.isfinite(cell).all(), 'Invalid contact cell')
    volume = float(np.linalg.det(cell))
    require(volume > 0, 'Nonpositive cell volume')
    shifts = np.array([(i, j, k) for i in range(-4, 5) for j in range(-4, 5)
                       for k in range(-4, 5) if (i, j, k) != (0, 0, 0)]) @ cell
    distances = np.linalg.norm(shifts, axis=1)
    require(bool(np.all(distances > 0)), 'Degenerate periodic image')
    image_sum = math.fsum(math.erfc(WIDTH * float(r)) / float(r) for r in distances)
    background = -math.pi / (WIDTH**2 * volume)
    return dict(C_eV_A=C, BOHR_A=BOHR_A, gaussian_width_inverse_A=WIDTH,
                real_image_shell=4, volume_A3=volume,
                B_source_eV=2 * C * WIDTH / math.sqrt(math.pi),
                nonzero_image_eV=C * image_sum, zero_G_background_eV=C * background,
                regular_images_eV=C * (image_sum + background))


def verify_mapping(cell_A, H_M_Hartree, He_M_Hartree, ASE_Hartree_eV,
                   ae_source, mapping_source):
    ae_binding = verify_operator_source(ae_source)
    mapping_binding = source_binding(mapping_source)
    require(mapping_binding['sha256'] == MAPPING_SHA256, 'Finite-R operator is not the accepted frozen source')
    mh = -math.sqrt(10. / (2 * math.pi))
    require(abs(H_M_Hartree - mh) <= 1e-12 and abs(He_M_Hartree - 4 * mh) <= 1e-12,
            'Wrong H/He quadratic Gaussian self reference')
    require(math.isfinite(ASE_Hartree_eV) and 27. < ASE_Hartree_eV < 28., 'Invalid ASE Hartree constant')
    mapping = regular_mapping(cell_A)
    bn = -2 * mh * ASE_Hartree_eV
    # Exact rational algebra of the independent Gaussian self-reference terms.
    self_coeff = Fraction(2) - Fraction(2)**2
    require(self_coeff == -2, 'Gaussian quadratic self-reference algebra failed')
    # Source-bound erfc Taylor coefficients: (erfc(wR)-1)/R.
    coefficients = [Fraction(-2), Fraction(2, 3), Fraction(-1, 5)]
    require(coefficients[0] == -2, 'Wrong signed singular-subtraction limit')
    cell = np.asarray(cell_A, dtype=float)
    shifts = np.array([(i, j, k) for i in range(-4, 5) for j in range(-4, 5)
                       for k in range(-4, 5) if (i, j, k) != (0, 0, 0)]) @ cell
    direction = np.ones(3) / math.sqrt(3.)
    limit = -mapping['B_source_eV'] + mapping['regular_images_eV']
    samples = []
    for radius in (1e-2, 1e-3, 1e-4, 1e-5):
        distances = np.linalg.norm(shifts - radius * direction, axis=1)
        # Stable independent evaluation of source expression minus C/R.
        regularized = C * (-math.erf(WIDTH * radius) / radius +
                          math.fsum(math.erfc(WIDTH * float(r)) / float(r) for r in distances)
                          - math.pi / (WIDTH**2 * mapping['volume_A3']))
        samples.append(dict(R_A=radius, regularized_eV=regularized,
                            limit_error_eV=regularized - limit))
    require(abs(samples[-1]['limit_error_eV']) <= 1e-7, 'Numerical continuum limit failed')
    require(abs(samples[-1]['limit_error_eV']) < abs(samples[0]['limit_error_eV']) / 1000,
            'Numerical continuum sequence does not approach the derived limit')
    wrong_linear_M = (2 * mh - 2 * mh) * ASE_Hartree_eV
    require(abs(bn - wrong_linear_M) > 1., 'Wrong-M control is not sensitive')
    return dict(schema='hydride-contact-mapping-proof-v1',
                executed=True, scope='Source algebra and numerical Gaussian mapping only; no Native or isolation validation',
                operator_source=ae_binding, finite_R_mapping_source=mapping_binding,
                derivation=dict(gaussian_pair_width='sqrt(alpha1/2)/BOHR_A',
                                self_reference_coefficient=str(self_coeff),
                                erfc_regular_limit='-2*width/sqrt(pi)',
                                erfc_regular_series_coefficients=[str(v) for v in coefficients],
                                W0='Q_He + B_native - B_source + regular_images'),
                mapping=mapping, B_native_eV=bn,
                unit_constant_residual_eV=bn - mapping['B_source_eV'],
                numerical_continuum_samples=samples,
                sensitivity=dict(wrong_linear_He_M_predicted_shift_eV=wrong_linear_M,
                                 wrong_linear_He_M_error_eV=wrong_linear_M - bn,
                                 one_eV_atomic_offset_identity_error_eV=-1.,
                                 atomic_offset_permitted=False))


def verify_operator_identity(radial, source_paths, guard_paths, convention, charge, electrons):
    """Check actual radial operators and arbitrary-common-density functional algebra.

    This is not an executed coincident-atom SCF. The public guard stays enabled.
    """
    require(radial['schema'] == 'hydride-contact-radial-v1', 'Wrong actual radial schema')
    require(set(source_paths) == {'ae', 'setup', 'basis_data', 'hamiltonian', 'density'} and
            set(guard_paths) == {'legacy_calculator', 'utilities'}, 'Missing assembly/guard source')
    trees = {k: ast.parse(Path(p).read_text()) for k, p in source_paths.items()}
    def has(key, text):
        target = ast.dump(ast.parse(text).body[0], include_attributes=False)
        return any(ast.dump(n, include_attributes=False) == target for n in ast.walk(trees[key]))
    require(has('hamiltonian', 'e_coulomb += setup.M + np.dot(D_p, (setup.M_p + np.dot(setup.M_pp, D_p)))'),
            'Hamiltonian self-reference assembly changed')
    require(has('hamiltonian', 'e_zero = self.finegd.integrate(self.vbar_g, dens.nt_g, global_integral=False)'),
            'Hamiltonian linear local-potential functional changed')
    require(has('hamiltonian', 'e_coulomb = 0.5 * self.finegd.integrate(self.vHt_g, dens.rhot_g, global_integral=False)'),
            'Hamiltonian Hartree functional changed')
    require(has('density', 'Q_L[0] += setup.Delta0') and
            has('density', 'self.ghat.add(self.rhot_g, self.Q_aL)'), 'Density compensation assembly changed')
    guard_trees = {k: ast.parse(Path(p).read_text()) for k, p in guard_paths.items()}
    require(any(isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and
                n.func.id == 'check_atoms_too_close' for n in ast.walk(guard_trees['legacy_calculator'])),
            'Missing public coincident-atom guard call')
    require(any(isinstance(n, ast.FunctionDef) and n.name == 'check_atoms_too_close' and
                any(isinstance(c, ast.Raise) for c in ast.walk(n))
                for n in ast.walk(guard_trees['utilities'])), 'Missing public coincident-atom guard implementation')
    r = np.asarray(radial['r_Bohr'], dtype=float)
    require(r.shape == (801,) and np.isfinite(r).all() and r[0] == 0. and r[-1] == 2. and
            bool(np.all(np.diff(r) > 0)), 'Wrong actual radial sampling')
    h, he = radial['H'], radial['He']
    empty = dict(Delta_pL=[0, 1], Delta_iiL=[0, 0, 1], dO_ii=[0, 0], dC_ii=[0, 0],
                 M_pp=[0, 0], K_p=[0], M_p=[0], MB_p=[0], X_p=[0])
    for atom, z in ((h, 1), (he, 2)):
        for key, expected in dict(Z=z, Nv=z, Nc=0, ni=0, E=0., alpha1=10., alpha2=300.,
                                  occupation=[float(z)], empty_shapes=empty, Kc=0., MB=0.,
                                  Nct=0., pt_count=0, type='all-electron',
                                  ExxC=0., xc_correction=None).items():
            require(atom[key] == expected, 'Wrong actual constructed operator: ' + key)
        require(abs(atom['Delta0'] + z / math.sqrt(4 * math.pi)) <= 1e-12, 'Wrong actual Delta0')
        require(abs(atom['M'] - z*z*convention['H_M_Hartree']) <= 1e-12, 'Wrong actual quadratic M')
    arrays, residuals = {}, {}
    for name in ('ghat0_radial', 'compensation_radial', 'vbar_radial', 'nct_radial', 'tauct_radial'):
        arrays[name] = tuple(np.asarray(a[name], dtype=float) for a in (h, he))
        require(all(a.shape == r.shape and np.isfinite(a).all() for a in arrays[name]), 'Invalid actual radial ' + name)
    def check(label, a, b):
        error = float(np.max(np.abs(a-b)))
        require(error <= 1e-9, 'Actual operator identity failed: ' + label)
        residuals[label] = error
    check('ghat_He_equals_H', arrays['ghat0_radial'][1], arrays['ghat0_radial'][0])
    for i, a in enumerate((h, he)):
        check('Delta0_times_ghat_' + str(i), arrays['compensation_radial'][i],
              a['Delta0'] * arrays['ghat0_radial'][i])
        for name in ('nct_radial', 'tauct_radial'):
            check(name + '_' + str(i), arrays[name][i], np.zeros_like(r))
    check('compensation_He_equals_two_H', arrays['compensation_radial'][1], 2*arrays['compensation_radial'][0])
    check('vbar_He_equals_two_H', arrays['vbar_radial'][1], 2*arrays['vbar_radial'][0])
    # Exact symbolic linear combination of *functional* terms, valid for any
    # common orbitals/density and same-q background, not a fitted scalar offset.
    hh = {'T[n]': Fraction(1), 'EH[n+2g]': Fraction(1), 'XC[n]': Fraction(1),
          'external[n,2v]': Fraction(1), 'M_H': Fraction(2), 'E_H': Fraction(2)}
    united = dict(hh, M_H=Fraction(4))
    difference = {key: hh[key]-united[key] for key in hh}
    require(difference['M_H'] == -2 and all(v == 0 for k, v in difference.items() if k != 'M_H'),
            'Arbitrary-density functional subtraction failed')
    hartree = convention['ASE_Hartree_eV']
    expected = (2*h['M'] - he['M']) * hartree
    independent = -2*convention['H_M_Hartree'] * hartree
    require(abs(expected-independent) <= 1e-9, 'Actual operator self-reference difference failed')
    weights, dr = 4*math.pi*r*r, np.diff(r)
    def integral(a):
        return float(np.sum((a[:-1]+a[1:])*dr/2))
    gaussian_norm = integral(r*r*arrays['ghat0_radial'][0])
    require(abs(gaussian_norm-1.) <= 1e-3, 'Constructed ghat is not normalized')
    trials = []
    for exponent in (.5, 1., 2.):
        density = np.exp(-exponent*r*r)
        density *= electrons/integral(weights*density)
        exact = integral(weights*density*(2*arrays['vbar_radial'][0]-arrays['vbar_radial'][1]))/math.sqrt(4*math.pi)
        wrong = integral(weights*density*(2*arrays['vbar_radial'][0]-1.01*arrays['vbar_radial'][1]))/math.sqrt(4*math.pi)
        require(abs(exact*hartree) <= 1e-9 and abs(wrong*hartree) > 1e-4, 'Actual-vbar sensitivity failed')
        trials.append(dict(exponent_Bohr_m2=exponent, electrons=electrons, residual_eV=exact*hartree,
                           wrong_He_vbar_scale=1.01, wrong_vbar_residual_eV=wrong*hartree))
    wrong_m = (2*h['M']-2*h['M'])*hartree-independent
    require(abs(wrong_m) > 1., 'Wrong-M sensitivity is not meaningful')
    require(type(charge) is int and charge in (0, 1) and electrons > 0, 'Wrong actual identity branch')
    wrong_ne = electrons+charge-(1-charge)
    return dict(schema='hydride-contact-functional-identity-replay-v1', charge=charge, electrons=electrons,
                proof_kind='source-bound actual-operator functional identity; not coincident Native SCF',
                coincident_Native_performed=False, public_guard_obstruction='Coincident atoms rejected; guard unchanged',
                sources={k: source_binding(p) for k, p in source_paths.items()},
                guard_sources={k: source_binding(p) for k, p in guard_paths.items()},
                actual_radial_max_abs_residuals=residuals,
                symbolic_functional_difference={k: str(v) for k, v in difference.items()},
                source_expected_difference_eV=independent, constructed_operator_difference_eV=expected,
                residual_eV=expected-independent, passed=True, density_functional_trials=trials,
                sensitivity=dict(wrong_linear_He_M_residual_eV=wrong_m, one_eV_He_offset_residual_eV=-1.,
                                 wrong_branch_electrons=wrong_ne, wrong_charge_electron_mismatch=wrong_ne-electrons,
                                 wrong_charge_would_be_rejected=True),
                scope='Common-density operator identity for actual charge/Ne; He SCF and h qualification separate')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, help='Contact calculation with frozen source and raw states')
    args = parser.parse_args()
    from contact_replay import load_contact
    print(json.dumps(load_contact(args.directory)['mapping_proof'], indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
