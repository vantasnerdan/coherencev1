"""Digitize Fig 2 error bars from paper/fig-001.png, starting at the marker centres in work/fig2_digitized.csv.
Two end-finders per marker: (w) walk up/down the marker's column (+-1 px) through pixels of the marker colour,
allowing gaps <= 6 px; (c) error-bar caps = rows where >= 11 of the 17 pixels centred on the column match the
colour, outside the marker body. Matplotlib bars are symmetric in linear R, so for each marker the (up, down)
half-widths pair, chosen from {w, c} x {w, c}, closest to symmetric is kept. Three bars partly hidden behind
other markers were read by eye (VISUAL). TiD at 1.0 keV is hidden behind Ti; its R is the midpoint of its caps.
Output: errbars.csv (E_cm, material, R, half-width, sigma_lnR = half/R floored at 0.05, method)."""
import numpy as np, os
from PIL import Image
S = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
im = np.asarray(Image.open(os.path.join(S, 'paper/fig-001.png')).convert('RGB')).astype(int)
XS = (1012-103.5)/6.5; Y0 = 51.5; YS = (736.5-51.5)/4
X0 = {'Pd': 103.5, 'PdD': 103.5, 'Ti': 1073.5, 'TiD': 1073.5}
COL = {'Pd': (105, 105, 105), 'PdD': (139, 0, 0), 'Ti': (128, 128, 128), 'TiD': (0, 0, 128)}
VISUAL = {(1.5, 'Ti'): (5.6e-4, 'visual: top cap hidden behind TiD marker; walk ends'),
          (0.5, 'TiD'): (8.1e-4, 'visual: blue caps, bar hidden under Ti bar'),
          (2.5, 'Ti'): (1.23e-3, 'visual: upper half from walk, symmetric')}
r2y = lambda R: Y0 - np.log10(R)*YS
y2r = lambda y: 10**(-(y - Y0)/YS)

def match(px, m, tol):
    return np.abs(px - np.array(COL[m])).max(-1) <= tol

def walk(x, y, m):
    ok = match(im[:, x-1:x+2], m, 9 if m in ('Pd', 'Ti') else 25).any(1); ends = []
    for step in (-1, 1):
        yy = last = y; gap = 0
        while 52 < yy < 736:
            yy += step
            if ok[yy]: last = yy; gap = 0
            else:
                gap += 1
                if gap > 6: break
        ends.append(last)
    return y2r(ends[0]), y2r(ends[1])

def caps(x, y, m, hidden=False):
    full = match(im[:, x-8:x+9], m, 8 if m in ('Pd', 'Ti') else 25).sum(1) >= 11
    cs = np.where(full)[0]
    if hidden:
        up = [c for c in cs if y-80 < c < y-12]; dn = [c for c in cs if y+12 < c < y+120]
        return y2r(max(up)), y2r(min(dn))
    top = bot = y
    while full[top-1] and y-top < 14: top -= 1
    while full[bot+1] and bot-y < 14: bot += 1
    up = [c for c in cs if 52 < c < top-1]; dn = [c for c in cs if bot+1 < c < 736]
    return y2r(max(up) if up else top), y2r(min(dn) if dn else bot)

rows = [l for l in open(os.path.join(S, 'work/fig2_digitized.csv')) if not l.startswith('#')]
hdr = rows[0].strip().split(','); data = [r.strip().split(',') for r in rows[1:]]
pts = [(float(d[0]), m, float(d[hdr.index(m)])) for m in ['Pd', 'PdD', 'Ti', 'TiD'] for d in data if d[hdr.index(m)] != '']
pts.append((1.0, 'TiD', None))
out = []
for E, m, R in pts:
    x = int(round(X0[m] + E*XS))
    if E >= 6.5:
        out.append((E, m, R, 0.0, 'ref')); continue
    if R is None:
        hi, lo = caps(x, int(round(r2y(1.1e-3))), m, hidden=True)
        R = 0.5*(hi + lo); out.append((E, m, R, 0.5*(hi - lo), 'cc asym=1.00')); continue
    y = int(round(r2y(R)))
    hw, lw = walk(x, y, m); hc, lc = caps(x, y, m)
    best = None
    for ku, u in (('w', hw - R), ('c', hc - R)):
        for kd, d in (('w', R - lw), ('c', R - lc)):
            if u > 0 and d > 0:
                a = abs(np.log(u/d))
                if best is None or a < best[0]: best = (a, 0.5*(u + d), ku + kd + ' asym=%.2f' % np.exp(a))
    half, meth = best[1], best[2]
    if (E, m) in VISUAL: half, meth = VISUAL[(E, m)]
    out.append((E, m, R, half, meth))
out.sort(key=lambda o: (['Pd', 'PdD', 'Ti', 'TiD'].index(o[1]), o[0]))
with open(os.path.join(S, 'agents/measurement/errbars.csv'), 'w') as f:
    f.write('# Fig 2 error bars digitized from paper/fig-001.png (digitize_errorbars.py).\n')
    f.write('# Ends found by column walk (w) and cap detection (c); the up/down pair closest to symmetric is kept\n')
    f.write('# (matplotlib bars are symmetric in linear R). half = mean half-width; sigma_lnR = half/R floored at 0.05.\n')
    f.write('# TiD 1.0 keV marker is hidden behind Ti; its R is the linear midpoint of its caps.\n')
    f.write('E_cm,material,R,half,sigma_lnR,method\n')
    for E, m, R, h, meth in out:
        sl = max(h/R, 0.05) if E < 6.5 else 0.0
        f.write('%.2f,%s,%.4g,%.3g,%.3f,%s\n' % (E, m, R, h, sl, meth))
