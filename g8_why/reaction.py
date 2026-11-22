"""Response-weighted positive reaction-law and retained threshold-term comparison."""
import json
import hashlib
import numpy as np
from scipy.optimize import least_squares
from common import (DD, ENERGIES, MATERIALS, OUT, U, response, yield_C, readout,
                    observed, paper_rows, baseline_rows, baseline_proof, write_csv)

# Piecewise-positive intrinsic S/S_BH over the full sampled reaction-energy range.
EDGES = np.array([.05, .10, .18, .28, .40, .60, .90, 1.30, 2., 3., 4.5, 5.5, 7.5])
NB = len(EDGES)-1


def basis(material, e, ch):
    lab, w = response(material, float(e))
    ec = np.maximum(lab/2, .05)
    j = np.clip(np.searchsorted(EDGES, ec, side='right')-1, 0, NB-1)
    return np.bincount(j, weights=w*DD.sigma(ec, U, ch), minlength=NB)


def resonance_parts(material, e, ch):
    lab, w = response(material, float(e))
    ec = np.maximum(lab/2, .05)
    base = DD.sigma(ec, U, ch)
    if ch == 'n':
        return np.array([w@base, 0., 0.])
    sb = DD.S_BH(ec, ch)
    sr = np.sqrt(2*986.)*197326.9804**3*40.e-6/(7.*937806.5**1.5*ec**2)*10.
    penetration = base/sb
    interference = 2*np.sqrt(sb*sr/3)*np.cos(np.deg2rad(115.))
    # a=sqrt(Gamma_p/40meV), keeping radius, phase, ER~0 convention fixed.
    return np.array([w@base, w@(penetration*interference), w@(penetration*sr)])


def run():
    baseline = baseline_rows()
    write_csv(OUT/'baseline.csv', baseline)
    cache = {}
    energy_sets = {}
    data = {}
    for m in MATERIALS:
        es = np.array([float(r['E_cm']) for r in paper_rows()
                       if r['material'] == m and float(r['E_cm']) < 6.5])
        energy_sets[m] = es
        data[m] = observed(m, es)
        for e in np.unique(np.r_[es, ENERGIES, 4.5, 6.5]):
            for ch in ('p', 'n'):
                cache[m, float(e), ch] = basis(m, e, ch)
    results, rows = {}, []

    def fit_basis(materials, power):
        def evaluate(x, m, es):
            mult = np.exp(np.r_[x, 0.])  # gauge: last energy bin =1, not measured S
            y = lambda e, ch: float(cache[m, e, ch]@mult)
            return readout(m, es, y, current_power=power)
        def residual(x):
            return np.concatenate([(np.log(evaluate(x, m, energy_sets[m]))-np.log(data[m][0]))/data[m][1]
                                   for m in materials])
        # Start from required response-weighted shape, not zero-gradient tiny tails.
        centers = np.sqrt(EDGES[:-1]*EDGES[1:])
        means = []
        for m in materials:
            obs, _ = observed(m, ENERGIES)
            base = readout(m, ENERGIES, lambda e, ch: yield_C(m, e, ch), current_power=power)
            means.append(np.interp(np.log(centers), np.log(ENERGIES), np.log(obs/base)))
        x0 = np.mean(means, axis=0)[:-1]
        fits = [least_squares(residual, start, bounds=(-80., 80.), max_nfev=400,
                             ftol=2e-10, xtol=2e-10, gtol=2e-8)
                for start in (x0, np.zeros(NB-1))]
        fit = min(fits, key=lambda f: f.fun@f.fun)
        sv = np.linalg.svd(fit.jac, compute_uv=False)
        rank = int(np.sum(sv > sv[0]*1e-7))
        return dict(chi2=float(fit.fun@fit.fun), n=sum(len(energy_sets[m]) for m in materials),
                    fitted_parameters=NB-1, jacobian_rank_relative_1e7=rank,
                    singular_values=sv.tolist(), positive_multipliers=np.exp(np.r_[fit.x, 0.]).tolist(),
                    success=bool(fit.success), optimality=float(fit.optimality)), evaluate, fit.x

    for power in (0, 1):
        for tag, materials in [('common', MATERIALS)]+[(m, (m,)) for m in MATERIALS]:
            result, evaluate, x = fit_basis(materials, power)
            results[f'basis_{tag}_k{power}'] = result
            mult = np.exp(np.r_[x, 0.])
            for m in materials:
                predicted = evaluate(x, m, ENERGIES)
                base = readout(m, ENERGIES, lambda e, ch: yield_C(m, e, ch), current_power=power)
                for i, e in enumerate(ENERGIES):
                    for ch in ('p', 'n'):
                        physical = float(cache[m, float(e), ch]@mult)/yield_C(m, e, ch)
                        rows.append(dict(model=f'basis_{tag}', current_power=power, material=m,
                            E_cm_keV=float(e), channel=ch, physical_same_E_factor=physical,
                            predicted_own65_R=float(predicted[i]), normalized_factor=float(predicted[i]/base[i])))

    rc = {(m, float(e), ch): resonance_parts(m, e, ch)
          for m in MATERIALS for e in np.unique(np.r_[energy_sets[m], ENERGIES, 4.5, 6.5]) for ch in ('p', 'n')}
    for power in (0, 1):
        for rejoin in (False, True):
            for tag, materials in [('common', MATERIALS)]+[(m, (m,)) for m in MATERIALS]:
                def values(loga, m, es):
                    a = float(np.exp(loga[0]))
                    y = lambda e, ch: float(rc[m, e, ch]@np.array([1., a, a*a]))
                    return readout(m, es, y, current_power=power, rejoin=rejoin)
                def residual(loga):
                    return np.concatenate([(np.log(values(loga, m, energy_sets[m]))-np.log(data[m][0]))/data[m][1]
                                           for m in materials])
                fit = least_squares(residual, [5.], bounds=(-30., 40.), max_nfev=120)
                a = float(np.exp(fit.x[0]))
                key = f'resonance_{tag}_k{power}_rejoin{rejoin}'
                results[key] = dict(chi2=float(fit.fun@fit.fun), n=sum(len(energy_sets[m]) for m in materials),
                    gamma_p_eV_conditional=.040*a*a, strength_relative_40meV=a*a,
                    phase_deg=115., radius_fm=7., ER_keV=0., total_width_identifiable=False,
                    success=bool(fit.success))
                if not rejoin:
                    pred = values(fit.x, materials[0], ENERGIES)
                    for m in materials:
                        pred = values(fit.x, m, ENERGIES)
                        base = readout(m, ENERGIES, lambda e, ch: yield_C(m, e, ch), current_power=power)
                        for i, e in enumerate(ENERGIES):
                            y = float(rc[m, float(e), 'p']@np.array([1., a, a*a]))
                            rows.append(dict(model=f'resonance_{tag}', current_power=power, material=m,
                                E_cm_keV=float(e), channel='p', physical_same_E_factor=y/yield_C(m,e),
                                predicted_own65_R=float(pred[i]), normalized_factor=float(pred[i]/base[i])))
    # Existing 40meV term is not its earlier high-fitted-reference plateau metric.
    retained = []
    for m in MATERIALS:
        y = lambda e, ch: float(rc[m, e, ch]@np.ones(3))
        shape = readout(m, ENERGIES, y)
        for i, e in enumerate(ENERGIES):
            retained.append(dict(material=m, E_cm_keV=float(e),
                proton_same_E_factor=y(e, 'p')/yield_C(m, e),
                total_DD_same_E_factor=(y(e,'p')+y(e,'n'))/(yield_C(m,e,'p')+yield_C(m,e,'n')),
                predicted_own65_R=float(shape[i])))
    write_csv(OUT/'reaction_floats.csv', rows)
    write_csv(OUT/'reaction_retained.csv', retained)
    shape_bounds = []
    for m in MATERIALS:
        lo, hi = rc[m, .25, 'p'], rc[m, 2.5, 'p']
        bound = (1+1/np.sqrt(3))/(1-1/np.sqrt(3))*max(lo[0]/hi[0], lo[2]/hi[2])
        obs, _ = observed(m, [.25, 2.5])
        required = obs[0]/obs[1]*DD.I_cup(.25)/DD.I_cup(2.5)
        shape_bounds.append(dict(material=m, any_strength_any_phase_p_shape_upper=float(bound),
            observed_p_shape_k1=float(required), required_over_upper=float(required/bound)))
    write_csv(OUT/'reaction_shape_bound.csv', shape_bounds)
    # Decision-targeted replay: basis partition and resonance at a=1 reproduce kernels.
    proof = baseline_proof()
    partition = max(abs(sum(cache[m, float(e), ch])/yield_C(m,e,ch)-1)
                    for m in MATERIALS for e in ENERGIES for ch in ('p','n'))
    res_parity = max(abs((rc[m,float(e),'p']@np.ones(3))/
                    yield_C(m,e,'p',lambda ec,ch: DD.reaction_S(ec,ch)/DD.S_BH(ec,ch))-1)
                    for m in MATERIALS for e in ENERGIES)
    assert partition < 3e-13 and res_parity < 3e-13
    proof.update(positive_basis_partition_max_relative=partition, retained_resonance_max_relative=res_parity,
                 reaction_script_sha256=hashlib.sha256((OUT/'reaction.py').read_bytes()).hexdigest())
    (OUT/'reaction_proof.json').write_text(json.dumps(proof, indent=2)+'\n')
    (OUT/'reaction_results.json').write_text(json.dumps(dict(energy_bin_edges_keV=EDGES.tolist(),
        convention='fixed smooth40 p/n4.5 calibration, retained s=1.0413535612834366; high-voltage focus; k0 primary, k1 alternative',
        fits=results, proton_only_shape_bounds=shape_bounds), indent=2)+'\n')
    print(json.dumps({'proof':proof, 'fits':results, 'baseline':baseline}, indent=2))


if __name__ == '__main__':
    run()
