"""Reduced ionic-response cause, projected through the unchanged target integrator.

Run from any directory: python g6_couple/response.py. No electronic response is
computed here. The point-charge and linear-response models are conditional.
"""
from pathlib import Path
import json
import sys
from contextlib import contextmanager
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from g4_beam.forward import model as M
from g4_beam.forward.consumer import forward_yields
from g4_beam.merit.api import reference
from g4_beam.merit.report import scaled_profile

HERE = Path(__file__).resolve().parent
C = 14.399645478  # eV Angstrom
QE = 1.602176634e-19
AMU = 1.66053906660e-27
MD = 2.01410 * AMU
EPS0 = 8.8541878128e-12
LAB = np.array([.5, 1., 2., 3., 5., 8., 13.])
PARAM = ROOT / 'g3_lattice/rerun/part2/updated/parameters.json'


def ionic_shift(lab_keV, a=4.084, screening_A=.541, frequency_THz=10.063741521544044,
                gamma_ps=2 / 4.136036943628514, domain_nm=np.inf):
    """Angular Doppler average, fixed one-charge-per-ion oscillator strength.

    Domain truncation only removes q < 2pi/L; it never multiplies by N.
    Full charge and a filled FCC sublattice are explicit optimistic assumptions.
    """
    density = 4 / a**3
    kappa = 1 / screening_A
    qmax = (6 * np.pi**2 * density)**(1 / 3)
    qmin = 0. if np.isinf(domain_nm) else 2 * np.pi / (10 * domain_nm)
    if qmin >= qmax:
        return 0.
    omega2 = (2 * np.pi * frequency_THz * 1e12)**2
    plasma2 = density * 1e30 * QE**2 / (EPS0 * MD)
    gamma = gamma_ps * 1e12
    speed = np.sqrt(2 * np.asarray(lab_keV) * 1000 * QE / MD)

    def one(v):
        def integrand(q):
            weight = q*q / (q*q + kappa*kappa)
            long2 = omega2 + plasma2 * weight
            if v == 0:
                angular = 1 / long2
            else:
                p = np.sqrt(long2 - gamma*gamma/4 + 0j)
                x = q * 1e10 * v
                angular = (np.arctanh(x / (p + .5j * gamma)) / (p * x)).real
            return 2*C/np.pi * weight * plasma2 * weight * angular
        return quad(integrand, qmin, qmax, epsabs=1e-10, limit=150)[0]
    result = np.array([one(v) for v in speed.ravel()]).reshape(speed.shape)
    return float(result) if result.ndim == 0 else result


def spectator_bound(amplitude_cap_A, a=4.084, screening_A=.541):
    """Triangle bound on DIFFERENCED spectator contact energy, not constant V.

    Every site's displacement must be pointwise <= A. This is not inferred
    from an RMS. The bound ignores shell cancellations and is not a signed Ue.
    """
    ijk = np.indices((17, 17, 17)).reshape(3, -1).T - 8
    basis = .5*np.array([[1, -1, 1], [1, 1, 0], [0, 0, 1]])
    r = np.linalg.norm(a * (ijk @ basis.T), axis=1)
    r = r[r > 1e-12]
    delta = 2 * amplitude_cap_A
    if delta >= r.min():
        raise ValueError('Displacement cap reaches the nearest reference site')
    potential = lambda x: C*np.exp(-x/screening_A)/x
    return float(np.sum(potential(r-delta)-potential(r)))


@contextmanager
def cross_section_adapter(transform):
    """Observe/transform only the reaction kernel; stopping/profile stay native."""
    original = M.reaction_sigma
    def adapted(ec, u=0., ch='p', reaction='proton'):
        return transform(original, ec, u, ch, reaction)
    M.reaction_sigma = adapted
    try:
        yield
    finally:
        M.reaction_sigma = original


def project(profile, *, extra_eV=0., dynamic=False, d2=.15, current=1., base_Ue_eV=40.):
    target = reference()
    # A compact interpolation of the smooth ionic susceptibility is evaluated
    # at EACH slowing-down node and each molecular constituent's local energy.
    grid = np.geomspace(.1, 16., 180)
    response = ionic_shift(grid) if dynamic else None
    def transform(original, ec, u, ch, reaction):
        shift = extra_eV
        if dynamic:
            shift = shift + np.interp(2*np.asarray(ec), grid, response)
        return original(ec, u + shift/1000, ch, reaction)
    with cross_section_adapter(transform):
        return forward_yields(LAB/2, profile, Ue_keV=base_Ue_eV/1000,
                              measurement={'D1_fraction':1-d2, 'D2_fraction':d2,
                                           'current_A':current},
                              stopping_fn=target.stopping_fn, reaction='proton')


def main():
    target = reference()
    cap = .0432
    bound = spectator_bound(cap)
    profiles = {'baseline':target.profiles,
                'near_surface_D_x50':scaled_profile(target.profiles, 50., depth=10e-9,
                                                     scale_ads=False)}
    scenarios = [('baseline', .15, 1.), ('near_surface_D_x50', .15, 1.),
                 ('D2_current_0.85', .85, 1.), ('half_current_frozen', .15, .5)]
    rows = []
    for name, d2, current in scenarios:
        profile = profiles.get(name, target.profiles)
        base = project(profile, d2=d2, current=current)
        frozen = project(profile, extra_eV=bound, d2=d2, current=current)
        dynamic = project(profile, dynamic=True, d2=d2, current=current)
        for i, energy in enumerate(LAB):
            rows.append(dict(intervention=name, E_lab_keV=float(energy),
                E_cm_keV=float(energy/2), base_Ue_eV=40.,
                baseline_physical_yield_per_C=float(base['physical_yield_per_C'][i]),
                signed_linear_locked_over_unlocked=1.,
                frozen_positive_envelope_factor=float(frozen['physical_yield_per_C'][i]/base['physical_yield_per_C'][i]),
                full_ionic_addition_factor=float(dynamic['physical_yield_per_C'][i]/base['physical_yield_per_C'][i]),
                frozen_normalized_R_factor=float(frozen['R'][i]/base['R'][i]),
                ionic_normalized_R_factor=float(dynamic['R'][i]/base['R'][i])))
    requirements = []
    for desired in (500., 2000.):
        acap = brentq(lambda a:spectator_bound(a)-desired, 0., 4.084/np.sqrt(2)/2*(1-1e-8))
        requirements.append(dict(desired_eV=desired, required_displacement_cap_A=acap,
            amplitude_over_nn=acap/(4.084/np.sqrt(2)),
            full_ionic_weight_multiplier_at_0p5keV=desired/ionic_shift(.5)))
    collision_reference = ionic_shift(2.)
    levers = []
    for desired in (500., 2000.):
        contrast = desired/collision_reference
        levers.append(dict(additional_screening_target_eV=desired,
            reference_E_lab_keV=2., reference_full_ionic_eV=collision_reference,
            additional_oscillator_weight_over_unit_reference=contrast,
            locked_over_unlocked_weight_if_unlocked_is_unit=1+contrast,
            locked_over_unlocked_effective_charge_if_weight_is_Z_squared=float(np.sqrt(1+contrast)),
            contact_potential_change_V=desired,
            fixed_shape_linear_polarization_contrast_over_spectator_envelope=desired/bound,
            interpretation='Targets only, with poles/kernel held fixed. Charge and weight changes require a new response calculation and sum-rule accounting. The spectator envelope is not a bound on electronic polarization.'))
    out = dict(model='Frozen screened spectator envelope plus fixed-strength ionic RPA; not electronic DFT',
        baseline='40 eV electronic Ue; unchanged transported PdD profiles and adsorbed sheet; proton threshold term retained',
        static_full_ionic_addition_eV=ionic_shift(0.),
        ionic_addition_eV=dict(zip(map(str, LAB), map(float, ionic_shift(LAB)))),
        finite_domain_ionic_eV_at_0p5keV={str(L):ionic_shift(.5,domain_nm=L) for L in (1.38,3.85,10.,100.)},
        frozen_pointwise_cap_A=cap, frozen_absolute_envelope_eV=bound,
        pure_Gamma_relative_spectator_change_eV=0.,
        susceptibility_locked_minus_unlocked_eV=0.,
        susceptibility_reason='Equal density, charge, curvature and loss give identical linear susceptibility. Coherent population is not additional oscillator strength.',
        requirements=requirements, projections=rows,
        response_lever_targets_at_2keV_lab=levers,
        conditional_domain_result='At fixed charge/density/poles, increasing L saturates at the full ionic addition; no finite domain attains 500-2000 eV in this model. Electronic/Floquet/nonlinear changes are a distinct open mechanism.',
        fourth_cause=dict(cause='lattice_state', fitted_parameters=None,
            near_surface_D_x50='Scale only top-10nm volumetric D; retain outer sheet; state-mediated component follows loading/coherence as well as density.',
            D2_fraction='Constituent energies and two D per molecular charge; no extra current-specific state factor at fixed pump and temperature.',
            half_current='Exactly unchanged physical yield/C for frozen target and independently maintained pump; beam-maintained state can relax and requires transport/pump history.',
            loading_temperature_signature='Correlated onset/loss with loading, phase-correlation length, pump threshold and temperature; hysteresis/relaxation tracks the lattice state rather than instantaneous current alone.',
            pair_at_rest='Pure Gamma displacement cancels. Pair rate requires joint relative coordinates plus a specified confining quantum potential.',
            DT='Recompute reduced mass, S factor and target isotope dynamics; the DD threshold interference is not transferable.',
            LiD='Recompute isotope modes, ionic density and electronic screening; no fitted LiD lattice state exists here.'),
        missing_electronic_response='Collision-scope paired frozen DFT snapshot with in-phase neighboring D versus undisplaced geometry; no substitute or zero is assigned here.')
    (HERE/'response.json').write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ('static_full_ionic_addition_eV', 'ionic_addition_eV',
          'frozen_absolute_envelope_eV','requirements')}, indent=2), flush=True)


if __name__ == '__main__':
    main()
