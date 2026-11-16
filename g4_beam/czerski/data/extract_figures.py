#!/usr/bin/env python3
"""Extract the colored experimental markers in arXiv 2409.02112v3.

Usage: python extract_figures.py /path/to/2409.02112v3.pdf > figure_points.csv
Requires pdftocairo. See README.md for calibration and limitations.
"""
import csv
import math
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

# Text tick centers, PDF points. y-values here are 10**exponent for log panels.
PANELS = {
    (5, 'a'): ('2a', 'Zr', 'yield', 116.04, 270.20, 0., 20., 367.52, 495.43, 3., -6., 'log'),
    (5, 'b'): ('2b', 'Zr', 'peak_channel', 333.56, 483.56, 0., 10., 395.77, 502.61, 200., 180., 'lin'),
    (6, 'a'): ('3a', 'Ti', 'yield', 110.77, 274.04, 0., 20., 381.16, 500.76, 3., -5., 'log'),
    (6, 'b'): ('3b', 'Ti', 'enhancement', 338.08, 480.89, 1., 6., 373.84, 501.36, 5., 0., 'log'),
    (12, 'a'): ('4a', 'Pd', 'yield', 127.41, 286.18, 0., 20., 283.37, 420.85, 3., -6., 'log'),
    (12, 'b'): ('4b', 'Pd', 'peak_channel', 347.33, 494.93, 0., 10., 297.23, 406.51, 305., 285., 'lin'),
    (23, 'a'): ('7a', 'Zr', 'yield', 119.44, 244.76, 2., 20., 294.99, 413.83, 3., -1., 'log'),
    (23, 'b'): ('7b', 'Zr', 'peak_channel', 351.13, 496.39, 0., 10., 298.29, 422.29, 204., 180., 'lin'),
}
BLUE = 'rgb(0%, 0%, 100%)'
RED = 'rgb(100%, 0%, 0%)'
GREENS = {'rgb(42.352295%, 65.490723%, 7.843018%)', 'rgb(0%, 57.647705%, 0%)'}


def markers(svg):
    root = ET.parse(svg).getroot()
    for element in root.iter():
        color = element.get('fill', '')
        if element.get('fill-rule') != 'evenodd' or color not in ({BLUE, RED} | GREENS):
            continue
        numbers = [float(n) for n in re.findall(r'-?\d+(?:\.\d+)?', element.get('d', ''))]
        if not numbers:
            continue
        xs, ys = numbers[::2], numbers[1::2]
        x, y = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        transform = element.get('transform')
        if transform:
            a, b, c, d, tx, ty = map(float, re.findall(r'-?\d+(?:\.\d+)?', transform))
            x, y = a*x + c*y + tx, b*x + d*y + ty
            w, h = abs(a*w), abs(d*h)
        if w < 7 and h < 7:
            yield color, x, y


def main(pdf):
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['figure_panel', 'target', 'beam_species', 'observable', 'energy_lab_keV_per_deuteron', 'energy_cm_keV', 'value', 'unit', 'x_pdf_pt', 'y_pdf_pt', 'uncertainty_state'])
    with tempfile.TemporaryDirectory() as tmp:
        for page in (5, 6, 12, 23):
            svg = Path(tmp) / f'p{page}.svg'
            subprocess.run(['pdftocairo', '-f', str(page), '-l', str(page), '-svg', pdf, str(svg)], check=True)
            seen = set()
            for color, x, y in markers(svg):
                if not 110 < x < 510:
                    continue
                side = 'a' if x < 300 else 'b'
                fig, target, obs, x1, x2, E1, E2, y1, y2, v1, v2, scale = PANELS[(page, side)]
                # Legend marks, not observations.
                if (page == 5 and side == 'a' and 176 < x < 181 and 420 < y < 447
                    or page == 5 and side == 'b' and 423 < x < 429 and 380 < y < 408
                    or page == 6 and side == 'a' and 166 < x < 171 and 448 < y < 473
                    or page == 6 and side == 'b' and 437 < x < 443 and 407 < y < 413
                    or page == 12 and side == 'a' and 211 < x < 216 and 356 < y < 384
                    or page == 12 and side == 'b' and 429 < x < 434 and 312 < y < 337
                    or page == 23 and side == 'a' and 167 < x < 172 and 374 < y < 379
                    or page == 23 and side == 'b' and 431 < x < 436 and 321 < y < 326):
                    continue
                # The common figure includes duplicate paths for a few exactly coincident markers.
                key = (fig, color, round(x, 2), round(y, 2))
                if key in seen:
                    continue
                seen.add(key)
                energy = E1 + (x-x1)*(E2-E1)/(x2-x1)
                value = v1 + (y-y1)*(v2-v1)/(y2-y1)
                if scale == 'log':
                    value = 10**value
                if fig == '3b':
                    beam = 'unspecified_in_panel'
                elif page == 23:
                    beam = 'D1+'
                else:
                    beam = 'D1+' if color == BLUE else ('D2+' if color in GREENS else 'D3+')
                if obs == 'peak_channel':
                    e_cm, e_lab = energy, 2*energy
                    unit = 'channel'
                else:
                    e_lab, e_cm = energy, energy/2
                    unit = 'counts/(microampere*kilosecond)' if obs == 'yield' else 'dimensionless'
                writer.writerow([fig, target, beam, obs, f'{e_lab:.4f}', f'{e_cm:.4f}', f'{value:.6g}', unit, f'{x:.3f}', f'{y:.3f}', 'statistical_error_bars_plotted_not_digitized'])


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
