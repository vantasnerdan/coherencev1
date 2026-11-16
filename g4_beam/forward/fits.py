"""Fit every (low-energy physics candidate) x (measurement hypothesis) to the digitized Fig 2 points.
Log-space chi^2 with digitized error bars (errbars.csv). Writes fits.csv and fits.json.
Run: python3 -I fits.py   (from this directory)"""
import json, os, sys, itertools
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model as M

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = M.load_data(os.path.join(HERE, 'errbars.csv'))
E_CUT = float(os.environ.get('E_CUT', 0.1))       # lab keV, lower end of slowing-down integral
TARGETS = {m: M.Target(m, E_cut=E_CUT) for m in M.MATS}

MEAS = {   # discrete measurement hypotheses
    'M0 as stated':            dict(k=0, b_hi=0, b_lo=0),
    'M1 fig x 1/I':            dict(k=1, b_hi=0, b_lo=0),
    'M2 cup focus >5.5kV':     dict(k=0, b_hi=1, b_lo=0),
    'M3 1/I + focus':          dict(k=1, b_hi=1, b_lo=0),
    'M4 cup fraction all V':   dict(k=0, b_hi=1, b_lo=1),
    'M5 1/I + fraction all V': dict(k=1, b_hi=1, b_lo=1),
}

# low-energy physics / floor candidates: (physics, neutral mode, extra parameter spec)
CANDS = {
    'none':          dict(phys='const'),
    'neutral d0':    dict(phys='const', neutral=0),
    'neutral d1':    dict(phys='const', neutral=1),
    'neutral5 d0':   dict(phys='const', neutral=0, En_fixed=5.0),
    'neutral5 d1':   dict(phys='const', neutral=1, En_fixed=5.0),
    'sigmoid_tt':    dict(phys='sigmoid_tt'),
    'sigmoid_xs':    dict(phys='sigmoid_xs'),
    'layer':         dict(phys='layer'),
    'czerski':       dict(phys='czerski'),
    'free floor':    dict(phys='const', floor=True),
}


def param_spec(c):
    """List of (name, init, lo, hi, log?) for candidate c."""
    sp = [('s', 1.0, 0.3, 3.0, True)]
    sp += [('Ue_%s' % m, u, 0.0, 15.0, False) for m, u in zip(M.MATS, (1.0, 1.8, 0.6, 1.2))]
    ph = c['phys']
    if ph in ('sigmoid_tt', 'sigmoid_xs'):
        sp += [('Ues_%s' % m, 3.0, 0.0, 30.0, False) for m in M.MATS] + [('E0', 0.6, 0.05, 3.0, False)]
    if ph == 'layer':
        sp += [('Ues_%s' % m, 4.0, 0.0, 60.0, False) for m in M.MATS] + [('d', 3.0, 0.2, 100.0, True)]
    if ph == 'czerski':
        sp += [('Gp', 40.0, 1e-3, 1e7, True)]
    if 'neutral' in c:
        sp += [('f_n', 0.05, 1e-5, 5.0, True)]
        if 'En_fixed' not in c:
            sp += [('E_n', 6.0, 0.5, 26.0, False)]
    if c.get('floor'):
        sp += [('F_%s' % m, 3e-3, 1e-7, 0.3, True) for m in M.MATS]
    return sp


def unpack(c, sp, v):
    d = {n: (np.exp(x) if lg else x) for (n, _, _, _, lg), x in zip(sp, v)}
    return d


def model_all(cname, mname, d):
    c = CANDS[cname]; out = {}
    for m in M.MATS:
        p = {'Ue': d['Ue_%s' % m]}
        if 'Ues_%s' % m in d: p['Ue_s'] = d['Ues_%s' % m]
        if 'E0' in d: p['E0'] = d['E0']
        if 'd' in d: p['d'] = d['d']
        if 'Gp' in d: p['Gp'] = d['Gp']
        meas = dict(MEAS[mname]); meas['s'] = d['s']
        if 'neutral' in c:
            meas.update(f_n=d['f_n'], delta=c['neutral'], E_n=c.get('En_fixed', d.get('E_n')))
        if c.get('floor'):
            meas['F'] = d['F_%s' % m]
        out[m] = M.relative_yield(TARGETS[m], DATA[m][0], c['phys'], p, meas)
    return out


def residuals(v, cname, mname, sp, emin=0.0):
    d = unpack(CANDS[cname], sp, v)
    R = model_all(cname, mname, d)
    r = []
    for m in M.MATS:
        E, y, s = DATA[m]; ok = E >= emin
        r.append(((np.log(y) - np.log(np.maximum(R[m], 1e-300)))/s)[ok])
    return np.concatenate(r)


def fit(cname, mname, emin=0.0, starts=4, seed=0, fix=None):
    c = CANDS[cname]; sp = param_spec(c)
    if fix:
        sp = [s for s in sp if s[0] not in fix]
    lo = np.array([np.log(s[2]) if s[4] else s[2] for s in sp])
    hi = np.array([np.log(s[3]) if s[4] else s[3] for s in sp])
    x0 = np.array([np.log(s[1]) if s[4] else s[1] for s in sp])
    rng = np.random.default_rng(seed); best = None

    def fun(v):
        d = unpack(c, sp, v)
        if fix: d.update(fix)
        R = model_all(cname, mname, d); r = []
        for m in M.MATS:
            E, y, s = DATA[m]; ok = E >= emin
            r.append(((np.log(y) - np.log(np.maximum(R[m], 1e-300)))/s)[ok])
        return np.concatenate(r)
    for i in range(starts):
        xi = x0 if i == 0 else np.clip(x0 + rng.normal(0, 1, len(x0))*(hi - lo)*0.15, lo + 1e-9, hi - 1e-9)
        try:
            r = least_squares(fun, xi, bounds=(lo, hi), x_scale='jac', max_nfev=3000)
        except Exception as e:
            continue
        if best is None or r.cost < best.cost:
            best = r
    d = unpack(c, sp, best.x)
    if fix: d.update(fix)
    chi2 = 2*best.cost; n = len(best.fun); kpar = len(sp)
    # 1-sigma from Jacobian (linearised), in natural units
    try:
        J = best.jac; cov = np.linalg.pinv(J.T @ J)
        err = np.sqrt(np.diag(cov))
        errs = {s[0]: (err[i]*d[s[0]] if s[4] else err[i]) for i, s in enumerate(sp)}
    except Exception:
        errs = {}
    return dict(cand=cname, meas=mname, chi2=chi2, n=n, k=kpar, dof=n - kpar,
                aic=chi2 + 2*kpar, bic=chi2 + kpar*np.log(n), params=d, errs=errs,
                per_mat={m: float(np.sum(fun(best.x)[_slices(emin)[m]]**2)) for m in M.MATS})


def _slices(emin):
    out = {}; i = 0
    for m in M.MATS:
        n = int(np.sum(DATA[m][0] >= emin)); out[m] = slice(i, i + n); i += n
    return out


def write_csv(results, path):
    keys = ['s'] + ['Ue_%s' % m for m in M.MATS] + ['Ues_%s' % m for m in M.MATS] + ['E0', 'd', 'Gp', 'f_n', 'E_n'] \
        + ['F_%s' % m for m in M.MATS]
    with open(path, 'w') as f:
        f.write('# Fits of the forward model to digitized Fig 2 (48 points, E_cm 0.25-6.0 keV; 6.5 is the reference).\n')
        f.write('# chi2 in ln R with digitized error bars; k = free parameters; AIC = chi2 + 2k. E_cut = %.3g keV lab (other cutoffs named in the candidate).\n' % E_CUT)
        f.write('# Ue in keV (Ue_* bulk, Ues_* surface/low-energy asymptote), d nm (mean projected depth), Gp meV (a = 8 fm),\n')
        f.write('# f_n neutral flux per unit current, E_n lab keV. chi2_<mat> = per-material contribution.\n')
        f.write('candidate,measurement,n,k,chi2,AIC,BIC,' + ','.join(keys) + ',' + ','.join('chi2_%s' % m for m in M.MATS) + '\n')
        for r in results:
            vals = [('%.4g' % r['params'][kk]) if kk in r['params'] else '' for kk in keys]
            f.write('%s,%s,%d,%d,%.2f,%.2f,%.2f,' % (r['cand'], r['meas'], r['n'], r['k'], r['chi2'], r['aic'], r['bic'])
                    + ','.join(vals) + ',' + ','.join('%.2f' % r['per_mat'][m] for m in M.MATS) + '\n')


if __name__ == '__main__':
    results = []
    for cname, mname in itertools.product(CANDS, MEAS):
        r = fit(cname, mname)
        results.append(r)
        print('%-12s %-24s chi2=%7.2f n=%d k=%2d AIC=%7.2f  ' % (cname, mname, r['chi2'], r['n'], r['k'], r['aic'])
              + ' '.join('%s=%.3g' % (kk, vv) for kk, vv in r['params'].items()), flush=True)
    # high-energy-only bulk fits (2.5-6.5 keV), constant Ue, s free, each measurement hypothesis
    for mname in MEAS:
        r = fit('none', mname, emin=2.5); r['cand'] = 'bulk 2.5-6.5'
        results.append(r)
        print('%-12s %-24s chi2=%7.2f n=%d ' % ('bulk2.5-6.5', mname, r['chi2'], r['n'])
              + ' '.join('%s=%.3g' % (kk, vv) for kk, vv in r['params'].items()), flush=True)
    # paper's sigmoid exactly as published (SI Table 2), no free parameters except s
    for ph in ('sigmoid_xs', 'sigmoid_tt'):
        fix = {}
        for m in M.MATS:
            us, ub, e0 = M.SI2[m]; fix['Ue_%s' % m] = ub; fix['Ues_%s' % m] = us; fix['E0'] = e0
        r = fit(ph, 'M0 as stated', fix=fix, starts=1); r['cand'] = ph + ' SI2 fixed'
        results.append(r)
        print('%-12s chi2=%7.2f n=%d s=%.3g' % (r['cand'], r['chi2'], r['n'], r['params']['s']))
    json.dump(results, open(os.path.join(HERE, 'fits.json'), 'w'), indent=1, default=float)
    write_csv(results, os.path.join(HERE, 'fits.csv'))
