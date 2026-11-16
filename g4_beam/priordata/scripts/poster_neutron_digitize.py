"""Digitize the neutron-rate plot of the Colborne et al. 2025 poster (OSTI 2569797, p.3, 'Initial Results').
Input: crop resized to 1460x1080 (neut.png). Axis calibration read by eye (estimate): y = 1e4 at 42 px, 1e1 at 985 px;
x = 2.5 keV at 184 px, 7.5 keV at 1373 px. Marker colours sampled from the image."""
import sys
import numpy as np
from PIL import Image
from scipy import ndimage

im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
cols = {'Ti': (0, 0, 0), 'TiDx': (42, 122, 128), 'Pd': (108, 168, 84), 'PdDx': (204, 170, 112)}
def y2v(py): return 10 ** (4 - 3 * (py - 42) / (985 - 42))
def x2e(px): return 2.5 + 5 * (px - 184) / (1373 - 184)
out = {}
for m, rgb in cols.items():
    mask = np.abs(im - np.array(rgb)).sum(axis=2) < (60 if m != 'Ti' else 90)
    mask[:200, :700] = False                       # legend
    mask[:, :140] = False; mask[1000:, :] = False    # axes labels
    mask = ndimage.binary_opening(mask, structure=np.ones((9, 9)))
    lab, n = ndimage.label(mask)
    for i in range(1, n + 1):
        ys, xs = np.nonzero(lab == i)
        e = round(x2e(xs.mean()) * 2) / 2
        out.setdefault(e, {})[m] = y2v(ys.mean())
print("E_cm_keV,Ti,TiDx,Pd,PdDx  # neutron counts per hour, estimate from plot")
for e in sorted(out):
    print(f"{e}," + ",".join(f"{out[e].get(m, float('nan')):.3g}" for m in cols))
