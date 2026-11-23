"""Distribution-free bounds for the fixed-profile bound-spectator reported yield.

Analytic inequalities apply to every isotropic spectral measure, not a K grid.
Only NEW outputs are written; the pinned source consumer owns transport/readout.
"""
import argparse
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import resource
import time

import numpy as np
import sympy as sp

from yield_kernel import DEFAULT_BASE, FixedProfileShellKernel, source_common

EXC_MAX = .100  # keV
DOMAIN_MAX = 16  # keV; checked against all actual local energies


def binding(path):
    path = Path(path)
    raw = path.read_bytes()
    return dict(path=str(path.resolve()), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def symbolic_receipt(dd, screening):
    """Check the square/flux identities and exact coefficient lower bounds."""
    a, b, mu, exc, rem = sp.symbols('a b mu exc rem', real=True)
    energy = (b*b-a*a)/2-a*b*mu+exc-rem
    square = b*b+exc-rem-(a-b)**2/2-a*b*(1+mu)
    residual = sp.expand(energy-square)
    E, U, L, BG = sp.symbols('E U L BG', positive=True)
    S = sp.Function('S')(E)
    sigma = S/sp.sqrt(E*(E+U))*sp.exp(-BG/sp.sqrt(E+U))
    reduced = sp.sqrt(2/L)*S/sp.sqrt(E+U)*sp.exp(-BG/sp.sqrt(E+U))
    flux_residual = sp.simplify(sigma*sp.sqrt(2*E/L)-reduced)
    log_derivative = sp.diff(sp.log(reduced), E)
    expected = sp.diff(S, E)/S+(BG-sp.sqrt(E+U))/(2*(E+U)**sp.Rational(3, 2))
    derivative_residual = sp.simplify(log_derivative-expected)
    # At mu=-1, Kcap is the upper root; substitute a=b+sqrt(2b²+2exc).
    root = b+sp.sqrt(2*b*b+2*exc)
    cap_residual = sp.simplify(energy.subs({a: root, mu: -1, rem: 0}))
    assert all(x == 0 for x in (residual, flux_residual, derivative_residual, cap_residual))
    coeff = {}
    D = Fraction(DOMAIN_MAX)
    for ch, values in dd.BH.items():
        aa = [Fraction(str(v)) for v in values]
        # Drop positive nonconstant terms; bound negative terms at D.
        lower_s = aa[0]+sum(min(c, 0)*D**i for i, c in enumerate(aa[1:], 1))
        lower_ds = aa[1]+sum(i*min(c, 0)*D**(i-1) for i, c in enumerate(aa[2:], 2))
        assert lower_s > 0 and lower_ds > 0
        coeff[ch] = dict(coefficients=[str(v) for v in aa],
                         S_lower_keV_mb=str(lower_s), S_prime_lower_mb=str(lower_ds))
    bg_gap = Fraction(str(dd.BG))**2-D-Fraction(str(screening))
    assert bg_gap > 0
    return dict(energy_residual=str(residual), flux_residual=str(flux_residual),
                log_derivative_residual=str(derivative_residual), Kcap_root_residual=str(cap_residual),
                domain_keV=[0, DOMAIN_MAX], BH_lower_bounds=coeff,
                BG_squared_minus_D_minus_U_keV=str(bg_gap),
                threshold_product_without_sqrt_2_over_L={
                    ch: float(values[0]/np.sqrt(screening)*np.exp(-dd.BG/np.sqrt(screening)))
                    for ch, values in dd.BH.items()})


def flux(common, energy, lab, channel):
    """Actual Equation 7 amplitude and relative/beam flux, positive E only."""
    return common.DD.sigma(energy, common.U, channel)*np.sqrt(2*energy/lab)


def envelope(common, material, nominal, removal, arrays):
    lab, path = common.response(material, nominal)
    ref_lab, ref_path = common.response(material, 6.5)
    assert np.all(lab > 0) and np.all(ref_lab > 0)
    assert np.all(path >= 0) and np.all(ref_path >= 0)
    assert max(lab.max(), ref_lab.max())+EXC_MAX < DOMAIN_MAX
    lmax = float(lab.max())
    cap = float((np.sqrt(lmax)+np.sqrt(2*lmax+2*EXC_MAX))**2)
    rem_max = lmax+EXC_MAX if removal else 0.
    upper_rate = flux(common, lab+EXC_MAX, lab, 'p')
    lower_energy = (ref_lab-cap)/2-rem_max
    positive = lower_energy > 0
    lower_rate = np.zeros_like(ref_lab)
    lower_rate[positive] = .5*flux(common, lower_energy[positive], ref_lab[positive], 'n')
    upper, lower = float(path@upper_rate), float(ref_path@lower_rate)
    assert upper > 0 and lower > 0
    # Invoke the real readout: no hand-reimplemented stitch/current normalization.
    def rates(e, ch):
        if e == nominal and ch == 'p':
            return upper
        if e == 6.5 and ch == 'n':
            return lower
        raise ValueError((e, ch))
    reported = {str(k): float(common.readout(material, [nominal], rates, current_power=k)[0])
                for k in (0, 1)}
    label = f'{material}_{nominal}_{int(removal)}'
    for key, value in dict(lab=lab, path=path, ref_lab=ref_lab, ref_path=ref_path,
                           upper_rate=upper_rate, lower_rate=lower_rate,
                           lower_energy=lower_energy).items():
        arrays[label+'_'+key] = value
    obs, error = common.observed(material, [nominal])
    return dict(material=material, CM_keV=nominal, nonnegative_common_removal=removal,
                local_Lmax_keV=lmax, Kcap_keV=cap, proton_active_removal_cap_keV=rem_max,
                numerator_upper_per_C=upper, neutron_reference_lower_per_C=lower,
                positive_reference_nodes=int(positive.sum()), reference_nodes=len(ref_lab),
                reported_upper=reported, observed=float(obs[0]), sigma_lnR=float(error[0]),
                observed_over_upper={k: float(obs[0]/r) for k, r in reported.items()},
                ln_observed_over_upper_over_sigma={k: float(np.log(obs[0]/r)/error[0])
                                                  for k, r in reported.items()})


def hot_counterexample(base, rows):
    """A free on-shell hot target violates the same reported envelope."""
    kernel = FixedProfileShellKernel(base, angular_order=64)
    energies = [.25, .5]
    # 10 keV lies above both proton-active compact caps.
    reported = {str(k): kernel.reported('PdD', energies, [10000.], [1.], current_power=k).tolist()
                for k in (0, 1)}
    limits = {r['CM_keV']: r for r in rows
              if r['material'] == 'PdD' and r['nonnegative_common_removal']}
    factors = {k: [value/limits[e]['reported_upper'][k] for e, value in zip(energies, values)]
               for k, values in reported.items()}
    assert all(v > 1 for values in factors.values() for v in values)
    return dict(model='free on-shell isotropic shell; not bound spectator', target_K_keV=10.,
                CM_keV=energies, angular_order=64, reported=reported,
                reported_over_broader_bound=factors,
                identity_failure='Efree=(L+K)/2-sqrt(L*K)*mu can exceed L+excmax; K is supplied initial energy')


def saved_consumer(common, rows, nodes):
    """Consume persisted nodes through the actual readout, closing the save bridge."""
    residuals = []
    with np.load(nodes) as saved:
        for row in rows:
            e = row['CM_keV']
            label = f"{row['material']}_{e}_{int(row['nonnegative_common_removal'])}"
            upper = float(saved[label+'_path']@saved[label+'_upper_rate'])
            lower = float(saved[label+'_ref_path']@saved[label+'_lower_rate'])
            def rates(energy, channel):
                if energy == e and channel == 'p':
                    return upper
                if energy == 6.5 and channel == 'n':
                    return lower
                raise ValueError((energy, channel))
            for k in (0, 1):
                replay = float(common.readout(row['material'], [e], rates, current_power=k)[0])
                residuals.append(abs(replay/row['reported_upper'][str(k)]-1))
    maximum = max(residuals)
    assert maximum <= 8*np.finfo(float).eps
    return dict(reported_values_consumed=len(residuals), max_relative_discrepancy=maximum)


def derivation():
    return [
        'Let L,K be nonnegative lab kinetic energies in keV, mu uniform on [-1,1], '
        '0<=exc<=0.100 and r>=0 a common state-dependent removal. The same normalized '
        'spectral measure is used at every nominal energy, species and transport depth.',
        'E=(L-K)/2-sqrt(L*K)*mu+exc-r = L+exc-r-(sqrt(K)-sqrt(L))^2/2'
        '-sqrt(L*K)*(1+mu) <= L+exc. Rates are zero for E<=0.',
        'For E>0, F(E,L)=sigma_Eq7(E,0.040)*sqrt(2E/L)='
        'sqrt(2/L)*S_BH(E)/sqrt(E+0.040)*exp(-BG/sqrt(E+0.040)). '
        'It has finite positive right limit at zero. The zero-forbidden-side extension is nondecreasing. '
        "dlogF/dE=S'/S+(BG-sqrt(E+U))/(2*(E+U)^(3/2))>0: "
        'exact rational BH coefficient bounds give S,Sprime>0 on [0,16], and BG^2>16+U.',
        'At mu=-1, E>0 implies sqrt(K)<sqrt(L)+sqrt(2L+2exc); hence any '
        'proton-active state has K<Kcap=(sqrt(Lmax)+sqrt(2Lmax+2excmax))^2. '
        'For arbitrary r>=0 it also has r<Lmax+excmax.',
        'Numerator for any proton-active state is <=A=sum_i w_i F(L_i+excmax,L_i). '
        'On the reference angular half mu in [-1,0], E>=(Lref-Kcap)/2-rmax, '
        'with rmax=0 or Lmax+excmax. Keep only positive reference nodes; '
        'B=sum_j wref_j F((Lref_j-Kcap)/2-rmax,Lref_j)/2>0.',
        'If q is the measure of proton-active states, N<=q*A and D>=q*B; '
        'the remaining states add zero N and nonnegative D. Thus N<=(A/B)*D '
        'for every normalized measure. With D>0 this bounds N/D independently of n(K). '
        'If D=0, N=0 and the reported ratio is undefined, not a finite prediction. '
        'Actual common.readout(A,B) supplies the unchanged stitch, focusing and both current powers.',
        'Premises: the finite pinned transport quadrature/profile, positive weights, isotropy '
        'conditional on spectral state, common energy/depth-independent spectral measure, '
        '0<=exc<=100eV, common nonnegative angle/beam-independent removal (possibly K-dependent recoil), '
        'equal masses/infinite spectator, and unchanged on-shell 40eV DD amplitude. '
        'This is not an exact many-body/off-shell nuclear theorem or a continuum transport enclosure. '
        'No site/density/profile/S-factor parameters are fitted.'
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--base-root', type=Path, default=DEFAULT_BASE)
    parser.add_argument('--prior-report', type=Path)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (180, 180))
    start = time.process_time()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    common = source_common(args.base_root)
    symbolic = symbolic_receipt(common.DD, common.U)
    arrays = {}
    rows = [envelope(common, mat, e, rem, arrays)
            for mat in common.MATERIALS for e in (.25, .5) for rem in (False, True)]
    sensitivity = hot_counterexample(args.base_root, rows)
    nodes = out/'envelope-nodes.npz'
    with nodes.open('xb') as stream:
        np.savez_compressed(stream, **arrays)
    saved_receipt = saved_consumer(common, rows, nodes)
    error_rows = list(csv.DictReader(line for line in
        (args.base_root/'g4_beam/forward/errbars.csv').read_text().splitlines()
        if not line.startswith('#')))
    for row in rows:
        source = next(r for r in error_rows if r['material'] == row['material']
                      and float(r['E_cm']) == row['CM_keV'])
        row['observed_linear_half_error'] = float(source['half'])
    root = Path('<PRIVATE_PATH>')
    paths = [Path(__file__), Path(__file__).with_name('yield_kernel.py'),
             Path(__file__).with_name('bound_forward.py'),
             root/'anharmonic-preregistration-001.json', root/'well-model-002.json',
             root/'well-model-003.json']
    sources = FixedProfileShellKernel(args.base_root).source_bindings()
    # Imported dynamic modules and any extra data files are also byte-bound.
    import sys
    for module in tuple(sys.modules.values()):
        file = getattr(module, '__file__', None)
        if file and Path(file).is_relative_to(args.base_root):
            paths.append(Path(file))
    sources['additional_inputs'] = [binding(p) for p in sorted(set(paths))]
    findings = [
        dict(id='energy-flux-monotonicity', disposition='established as stated',
             evidence='Exact symbolic residuals and rational coefficient inequalities; symbolic receipt.'),
        dict(id='distribution-free-envelope', disposition='established as stated',
             evidence='Analytic measure inequality plus saved fixed-quadrature sums and actual common.readout.'),
        dict(id='nonnegative-removal-extension', disposition='established as stated',
             evidence='Proton-active r<Lmax+excmax and reference angular-half bound; same removal across beams.'),
        dict(id='hot-target-sensitivity', disposition='established as stated', evidence=sensitivity),
        dict(id='zero-reference-degeneracy', disposition='scope limitation retained',
             evidence='All mass at sufficiently large K has N=D=0. Universal statement is N<=C*D; '
                      'reported ratio requires D>0. No zero denominator is replaced.'),
        dict(id='empirical-comparison', disposition='computed',
             evidence='Zero-removal PdD upper bounds at both energies are below observation for both '
                      'current powers. Broader arbitrary nonnegative-removal envelope still excludes '
                      'the0.25keV point by millions, but does not exclude the0.5keV observation for '
                      'current_power1 (upper0.01985 versus0.01668); current_power0 gap is only1.87 '
                      'sigma_lnR. These standardized log gaps are diagnostics, not discovery significances.'),
        dict(id='accuracy-and-physics', disposition='limitations retained',
             evidence='Floating-point fixed-quadrature sums, not outward-rounded continuum transport '
                      'enclosures; no new transport discretization or experimental uncertainty reduction. '
                      'Excitations above100eV and dynamical/energy-dependent off-shell physics are excluded.')]
    prior_cpu = 0.
    if args.prior_report:
        prior_report = json.loads(args.prior_report.read_bytes())
        prior_cpu = prior_report['runtime']['total_process_cpu_seconds']
        findings.append(dict(id='sensitivity-pilot-scope-correction', disposition='corrected',
                             evidence='First actual calculation used a5keV free shell and already broke '
                                      'both reported bounds, but its comment incorrectly said5keV exceeded '
                                      'both caps (the0.5keV cap is6.520keV). Current10keV shell exceeds '
                                      'both. Original complete findings/results are retained unchanged.',
                             original=binding(args.prior_report),
                             original_findings=prior_report['findings']))
    report = dict(schema='g8-spectator-envelope-v1', route_verdict='established as stated',
                  evidence_scope='Exact model inequalities with audited fixed-source numerical readout',
                  sources=sources, outputs=[binding(nodes)], symbolic=symbolic,
                  derivation=derivation(), results=rows, sensitivity=sensitivity, findings=findings,
                  saved_consumer=saved_receipt,
                  prior_numerical_cpu_seconds=prior_cpu,
                  cumulative_conservative_charge_seconds=np.ceil(prior_cpu)+np.ceil(
                      resource.getrusage(resource.RUSAGE_SELF).ru_utime+
                      resource.getrusage(resource.RUSAGE_SELF).ru_stime)+2,
                  runtime=dict(cpu_seconds=time.process_time()-start,
                               total_process_cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+
                               resource.getrusage(resource.RUSAGE_SELF).ru_stime,
                               cpu_limit_seconds=180, affinity=sorted(__import__('os').sched_getaffinity(0)),
                               maximum_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                               numpy_version=np.__version__, sympy_version=sp.__version__))
    final = out/'spectator-envelope-final.json'
    with final.open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(final=binding(final), runtime=report['runtime'],
                          PdD=[r for r in rows if r['material'] == 'PdD']), indent=2), flush=True)


if __name__ == '__main__':
    main()
