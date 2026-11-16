"""Follow-up checks: (1) layer model with wide bounds (d down to 1e-6 nm, Ue_s up to 200 keV);
(2) Czerski fit vs low-energy cutoff E_cut; (3) sanity of the layer integral; (4) Jacobian errors of bulk Ue."""
import json, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model as M, fits as F
T = M.Target('PdD'); E0 = np.array([0.5, 2., 5., 13.])
print('layer sanity: Ue_s=Ue_b ->', np.allclose(T.Y_layer(E0, 1.4, 1.4, 3.), T.Y_const(E0, 1.4), rtol=1e-6),
      '; d huge ->', np.allclose(T.Y_layer(E0, 3., 1.4, 1e4), T.Y_const(E0, 3.), rtol=1e-6))
orig = F.param_spec
def wide(c):
    sp = orig(c)
    if c['phys'] == 'layer':
        sp = [(n, (40.0 if n.startswith('Ues') else i), (lo), (200.0 if n.startswith('Ues') else hi), lg) if n.startswith('Ues') else
              ((n, 1e-3, 1e-7, 100.0, True) if n == 'd' else (n, i, lo, hi, lg)) for (n, i, lo, hi, lg) in sp]
    return sp
F.param_spec = wide
out = []
for mn in ('M0 as stated', 'M4 cup fraction all V'):
    r = F.fit('layer', mn, starts=6); r['cand'] = 'layer wide bounds'; out.append(r)
    print('layer wide %-22s chi2=%.2f ' % (mn, r['chi2']) + ' '.join('%s=%.3g' % kv for kv in r['params'].items()))
F.param_spec = orig
for ec in (0.05, 0.2, 0.4):
    F.TARGETS = {m: M.Target(m, E_cut=ec) for m in M.MATS}
    for mn in ('M0 as stated', 'M4 cup fraction all V'):
        r = F.fit('czerski', mn, starts=4); r['cand'] = 'czerski Ecut=%g' % ec; out.append(r)
        print('czerski Ecut=%.2f %-22s chi2=%.2f ' % (ec, mn, r['chi2']) + ' '.join('%s=%.3g' % kv for kv in r['params'].items()))
F.TARGETS = {m: M.Target(m, E_cut=0.1) for m in M.MATS}
json.dump(out, open('extra_checks.json', 'w'), indent=1, default=float)
res = json.load(open('fits.json'))
print('\nJacobian 1-sigma on bulk Ue (keV):')
for r in res:
    if r['cand'] in ('bulk 2.5-6.5', 'neutral d1', 'sigmoid_tt', 'neutral5 d1', 'free floor') and r['meas'] in ('M0 as stated', 'M4 cup fraction all V', 'M5 1/I + fraction all V'):
        print('%-14s %-24s ' % (r['cand'], r['meas']) + ' '.join('%s=%.2f+-%.2f' % (m, r['params']['Ue_%s' % m], r['errs'].get('Ue_%s' % m, np.nan)) for m in M.MATS)
              + '  s=%.2f+-%.2f' % (r['params']['s'], r['errs'].get('s', np.nan))
              + ('  f_n=%.3g+-%.2g' % (r['params']['f_n'], r['errs']['f_n']) if 'f_n' in r['params'] else '')
              + ('  E_n=%.2f+-%.2f' % (r['params']['E_n'], r['errs']['E_n']) if 'E_n' in r['params'] else ''))
