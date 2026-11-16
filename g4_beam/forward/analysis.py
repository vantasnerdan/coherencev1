"""Plots and the absolute-yield check, from fits.json.
Outputs: fig_fits.png, fig_ue.png, abs_yield.csv, ranges.csv."""
import json, os, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model as M
import fits as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = json.load(open(os.path.join(HERE, 'fits.json')))
RES += json.load(open(os.path.join(HERE, 'extra_checks.json')))       # layer with wide bounds, Czerski vs E_cut
F.write_csv(RES, os.path.join(HERE, 'fits.csv'))
T = F.TARGETS
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
SER = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
LS = ['-', '--', '-.', ':', (0, (5, 1, 1, 1)), (0, (1, 1)), '-', '--']


def best_of(cands):
    rs = [r for r in RES if r['cand'] in cands]
    return min(rs, key=lambda r: r['aic'])


def model_curve(r, m, E, physics_only=False):
    """Evaluate a fitted candidate on an arbitrary energy grid for material m.
    physics_only: drop the measurement layer (no charge tilt, no neutrals, s = 1); a free floor is kept."""
    c = F.CANDS[r['cand'].replace(' SI2 fixed', '').replace(' wide bounds', '')]; d = r['params']
    p = {'Ue': d['Ue_%s' % m]}
    for kk, kn in (('Ues_%s' % m, 'Ue_s'), ('E0', 'E0'), ('d', 'd'), ('Gp', 'Gp')):
        if kk in d: p[kn] = d[kk]
    meas = dict(F.MEAS[r['meas']]); meas['s'] = d['s']
    if 'neutral' in c:
        meas.update(f_n=d['f_n'], delta=c['neutral'], E_n=c.get('En_fixed', d.get('E_n')))
    if c.get('floor'):
        meas['F'] = d['F_%s' % m]
    if physics_only:
        meas = {'F': meas['F']} if 'F' in meas else {}
    return M.relative_yield(T[m], E, c['phys'], p, meas)


# ---------------- ranges ----------------
with open(os.path.join(HERE, 'ranges.csv'), 'w') as f:
    f.write('# D ranges (nm): CSDA path length and mean projected range (first-moment transport); paper quotes 4.4 nm (E_cm 0.25) and 67.5 nm (6.5)\n')
    f.write('# integrals start at 10 eV lab (E_cut = 0.01 keV)\n')
    f.write('material,E_cm,E_lab,path_nm,projected_nm\n')
    for m in M.MATS:
        Tr = M.Target(m, E_cut=0.01)
        for e in (0.25, 0.5, 1.0, 2.0, 2.5, 4.0, 6.5):
            f.write('%s,%.2f,%.1f,%.1f,%.1f\n' % (m, e, 2*e, Tr.path_range(2*e), Tr.proj_range(2*e)))

# ---------------- figure 1: data and fits ----------------
show = [('none', ['none']), ('energetic neutral', ['neutral d0', 'neutral d1']),
        ('neutral at 5 keV lab', ['neutral5 d0', 'neutral5 d1']),
        ('sigmoid Ue(E) thick-target', ['sigmoid_tt']), ('surface layer Ue(z)', ['layer', 'layer wide bounds']),
        ('Czerski resonance', ['czerski']), ('free floor', ['free floor'])]
Eg = np.linspace(0.2, 6.5, 200)
fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True)
for ax, m in zip(axs.ravel(), M.MATS):
    E, y, s = F.DATA[m]
    for i, (lab, cs) in enumerate(show):
        r = best_of(cs)
        ax.plot(Eg, model_curve(r, m, Eg), color=SER[i], ls=LS[i], lw=1.8,
                label='%s [%s] AIC %.0f' % (lab, r['meas'].split()[0], r['aic']))
    ax.axvline(4.25, color=GRID, lw=3, zorder=0)
    rp = [r for r in RES if r['cand'] == 'sigmoid_xs SI2 fixed'][0]
    ax.plot(Eg, model_curve(rp, m, Eg), color=INK2, ls=(0, (2, 2)), lw=1.2, label='paper Eq 16, SI Table 2 (cross-section ratio)')
    ax.errorbar(E, y, yerr=y*s, fmt='s', ms=5, color=INK, ecolor=INK, elinewidth=1, capsize=2, zorder=5, label='Fig 2 (digitized)')
    ax.set_yscale('log'); ax.set_ylim(1e-4, 1.5); ax.set_xlim(0, 6.7)
    ax.set_title(m, color=INK, loc='left', fontsize=12)
    ax.grid(True, color=GRID, lw=0.6); ax.tick_params(colors=INK2)
    for sp in ax.spines.values(): sp.set_color(INK2)
for ax in axs[1]: ax.set_xlabel('E_cm (keV)', color=INK)
for ax in axs[:, 0]: ax.set_ylabel('R = Y(E)/Y(6.5 keV)', color=INK)
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, fontsize=8.5, frameon=False)
fig.suptitle('Fig 2 data with the best-AIC fit of each low-energy candidate (measurement hypothesis in brackets).\nThe step at 4-4.5 keV (grey band) is the fitted proton-to-neutron stitch factor s.', color=INK, fontsize=10.5)
fig.tight_layout(rect=(0, 0.1, 1, 0.97))
fig.savefig(os.path.join(HERE, 'fig_fits.png'), dpi=130)

# ---------------- figure 2: per-energy Ue implied ----------------
def tt_R(m, U, E):
    """Analyst's model: constant Ue thick target, beam-averaged, as-stated normalisation, no stitch."""
    return M.relative_yield(T[m], np.atleast_1d(E), 'const', {'Ue': U}, {})


def invert_tt(m, E, R):
    out = []
    for e, rr in zip(np.atleast_1d(E), np.atleast_1d(R)):
        g = lambda U: np.log(tt_R(m, U, e)[0]) - np.log(rr)
        try:
            out.append(brentq(g, 0.0, 40.0, xtol=1e-3))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def invert_xs(E, R):
    out = []
    for e, rr in zip(np.atleast_1d(E), np.atleast_1d(R)):
        g = lambda U: np.log(M.sigma(e, U)/M.sigma(6.5, U)) - np.log(rr)
        try:
            out.append(brentq(g, 0.0, 40.0, xtol=1e-3))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


Ei = np.array([0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True)
ue_rows = []
for ax, m in zip(axs.ravel(), M.MATS):
    E, y, s = F.DATA[m]
    for i, (lab, cs) in enumerate(show):
        r = best_of(cs)
        U = invert_tt(m, Ei, model_curve(r, m, Ei, physics_only=True))
        ax.plot(Ei, U, color=SER[i], ls=LS[i], lw=1.8, label=lab)
        ue_rows.append((m, lab, r['meas'], U))
    us, ub, e0 = M.SI2[m]
    ax.plot(Eg, M.Ue_sigmoid(Eg, us, ub, e0), color=INK2, ls=(0, (2, 2)), lw=1.2, label='paper Eq 16 (SI Table 2)')
    ax.plot(E, invert_tt(m, E, y), 's', ms=6, color=INK, label='data, thick-target inversion')
    ax.plot(E, invert_xs(E, y), 'o', ms=6, mfc='none', mec=INK2, label='data, cross-section-ratio inversion (paper method)')
    ax.set_ylim(0, 8); ax.set_xlim(0, 6.3)
    ax.set_title(m, color=INK, loc='left', fontsize=12)
    ax.grid(True, color=GRID, lw=0.6); ax.tick_params(colors=INK2)
    for sp in ax.spines.values(): sp.set_color(INK2)
for ax in axs[1]: ax.set_xlabel('E_cm (keV)', color=INK)
for ax in axs[:, 0]: ax.set_ylabel('implied per-energy Ue (keV)', color=INK)
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, fontsize=8.5, frameon=False)
fig.suptitle('Per-energy Ue implied by each candidate\'s physics (fitted measurement layer removed), by constant-Ue thick-target inversion.\nPoints: the raw data inverted the same way (filled) and by the paper\'s cross-section ratio (open).', color=INK, fontsize=10.5)
fig.tight_layout(rect=(0, 0.1, 1, 0.97))
fig.savefig(os.path.join(HERE, 'fig_ue.png'), dpi=130)
with open(os.path.join(HERE, 'ue_implied.csv'), 'w') as f:
    f.write('# Ue (keV) implied by each candidate after removing its fitted measurement layer (constant-Ue thick-target inversion).\n')
    f.write('material,candidate,measurement,' + ','.join('U%.2f' % e for e in Ei) + '\n')
    for m, lab, ms, U in ue_rows:
        f.write('%s,%s,%s,' % (m, lab, ms) + ','.join('%.2f' % u for u in U) + '\n')
    for m in M.MATS:
        E, y, s = F.DATA[m]
        f.write('%s,data tt-inversion,,' % m + ','.join('%.2f' % u for u in invert_tt(m, E, y)) + ',E=' + ' '.join('%.2f' % e for e in E) + '\n')
        f.write('%s,data xs-inversion,,' % m + ','.join('%.2f' % u for u in invert_xs(E, y)) + '\n')

# ---------------- absolute yields (SI Table 1) ----------------
E_CHARGE = 1.602177e-19
rows = []
bulk = {}
for r in RES:
    if r['cand'] == 'bulk 2.5-6.5':
        bulk[r['meas']] = r['params']
bestfull = best_of([c for c in F.CANDS])
for m in M.MATS:
    cpc = M.counts_per_C_table(m)
    for lab, U in [('bare', 0.0), ('paper Ue_b', M.SI2[m][1]),
                   ('tt bulk M0', bulk['M0 as stated']['Ue_%s' % m]),
                   ('tt bulk M4', bulk['M4 cup fraction all V']['Ue_%s' % m]),
                   ('tt bulk M5', bulk['M5 1/I + fraction all V']['Ue_%s' % m]),
                   ('best full fit (%s, %s)' % (bestfull['cand'], bestfull['meas'].split()[0]), bestfull['params']['Ue_%s' % m])]:
        C = T[m].cumulative(U, 'p')
        for e in (3.0, 2.5):
            Yb = M.beam_average(lambda x: T[m].Y_const(x, U, C=C), np.array([2*e]))[0]   # per unit charge (D/e)
            per_x = Yb/T[m].x
            i = list(M.SI1_E).index(e)
            xeps = cpc[i]*E_CHARGE/per_x
            rows.append((m, lab, U, e, cpc[i], per_x/E_CHARGE, xeps, xeps*T[m].n))
with open(os.path.join(HERE, 'abs_yield.csv'), 'w') as f:
    f.write('# Counts/C from SI Table 1 (S:B x 0.007 cph / (I x 3600 s)); model counts/C per unit (x eps) = Y_beam/x / e,\n')
    f.write('# x = D per host atom near the surface, eps = proton detection efficiency (solid angle/4pi x transmission).\n')
    f.write('material,screening,Ue_keV,E_cm,counts_per_C,model_counts_per_C_per_xeps,x_eps,nD_eps_cm3\n')
    for row in rows:
        f.write('%s,%s,%.3f,%.1f,%.3g,%.4g,%.3g,%.3g\n' % row)
for row in rows:
    print('%-4s %-28s Ue=%.2f E=%.1f cpc=%.3g  x*eps=%.3g  nD*eps=%.3g' % (row[0], row[1], row[2], row[3], row[4], row[6], row[7]))

# D/M implied relative to TiD if the detection efficiency is common (3.0 keV point)
print('\nx relative to TiD (common eps), 3.0 keV:')
for lab in ('bare', 'tt bulk M0', 'tt bulk M4', 'tt bulk M5'):
    xe = {row[0]: row[6] for row in rows if row[1] == lab and row[3] == 3.0}
    print('%-12s ' % lab + ' '.join('%s=%.3g' % (m, xe[m]/xe['TiD']) for m in M.MATS))
