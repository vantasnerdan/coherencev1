"""Digitize the 150 visible Figure 3 samples from Chen et al. publisher PDF page 4.

Input: page converted with `pdftocairo -svg -f 4 -l 4 chen.pdf chen-page4.svg`.
The vertical bars are centered on the visible markers. Read their center and
half width instead of measuring obscured pixels under overlapping markers.
"""
import argparse
import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from compile_external import FIELDS

SVG = '{http://www.w3.org/2000/svg}'
NUM = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?')
SERIES = {'rgb(25.881958%, 25.881958%, 25.881958%)': 'beam-only control',
          'rgb(58.039856%, 6.666565%, 0%)': 'loading sequence'}
# Vector PDF points for each panel's 0 and 120 minute ticks and phase I/II edges.
X_TICKS = [(74.239, 217.503, 80.212, 151.844), (241.865, 385.131, 247.228, 325.426),
           (409.493, 552.757, 415.468, 488.285)]
Y_ZERO, Y_150 = 207.955, 83.076


def get_rows(svg):
    scale_y = 150 / (Y_ZERO - Y_150)
    out = {}
    for el in ET.parse(svg).iter(f'{SVG}path'):
        series = SERIES.get(el.get('stroke'))
        path = el.get('d', '')
        matrix = el.get('transform', '')
        if series is None or not matrix or not re.fullmatch(r'M[^A-Z]+L[^A-Z]+', path):
            continue
        a, b, c, d, e, f = map(float, NUM.findall(matrix))
        x1, y1, x2, y2 = map(float, NUM.findall(path))
        if b or c or abs(x1-x2) > .001:
            continue
        x = e+a*x1
        y = f+d*(y1+y2)/2
        bar = scale_y*abs(d*(y1-y2))/2
        if not 49 < y < 212.2:
            continue
        for panel, (x0, x120, phase_i, phase_ii) in enumerate(X_TICKS, 1):
            elapsed = (x-x0)*120/(x120-x0)
            bin_i = round((elapsed-2.5)/5)
            time = 2.5+5*bin_i
            # PDF also contains cropped-out drawing paths off the visible half-bin grid.
            if not 0 <= bin_i < 25 or abs(elapsed-time) > .04:
                continue
            phase = 'background' if x < phase_i else ('I' if x < phase_ii or series=='beam-only control' else 'II')
            key = (panel, series, time)
            if key in out:
                raise ValueError(f'Duplicate visible sample {key}')
            out[key] = dict(row_id=f'chen2025:fig3:{panel}:{series}:{bin_i+1:02}', study_id='chen2025',
                            lab='UBC/TRIUMF', run_id=f'Fig 3 campaign {panel} {series}', reaction='D+D', channel='n+3He',
                            target='Pd', host='Pd', target_geometry='thick', foil_thickness_um=300,
                            target_preparation='Pd membrane; not the same run before and after the campaign anneal',
                            E_cm_keV=15 if phase!='background' else '', E_lab_keV=30 if phase!='background' else '',
                            energy_frame_source='Nominal D+ sheath lab energy, CM=lab/2' if phase!='background' else '',
                            elapsed_plot_min=time, treatment_phase=phase, observable=f'{(Y_ZERO-y)*scale_y:.2f}',
                            observable_kind='neutron_rate_background' if phase=='background' else 'neutron_rate_time_bin',
                            observable_unit='neutrons/s', sigma_minus=f'{bar:.2f}', sigma_plus=f'{bar:.2f}',
                            evidence_class='figure_digitized_measurement', citation='https://doi.org/10.1038/s41586-025-09042-7',
                            source_locator=f'Fig. 3 panel {panel}, PDF point ({x:.3f},{y:.3f})',
                            conditions_note='Error columns are plotted vertical half-widths; confidence level not specified in figure caption. The paired sequence has a 0.5 mA plasma setpoint in phases I and II, not a delivered D+ current trace. Repeated bins within a run are not replicates of the treatment.')
    if set(out) != {(i, s, 2.5+5*k) for i in range(1,4) for s in SERIES.values() for k in range(25)}:
        raise ValueError('Missing visible Fig. 3 sample')
    return [out[k] for k in sorted(out)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source_svg', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    rows = get_rows(args.source_svg)
    with args.output.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS + ['elapsed_plot_min', 'treatment_phase'])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
