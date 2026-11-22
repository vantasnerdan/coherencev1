"""Report conditional absolute model demands and actual SI1 detected counts/C."""
import csv
from common import OUT, DD, MATERIALS, write_csv


def rows():
    with (OUT/'baseline.csv').open() as stream:
        baseline = list(csv.DictReader(stream))
    out = []
    for row in baseline:
        m, e = row['material'], float(row['E_cm_keV'])
        channel = 'p' if e <= 4. else 'n'
        y = float(row[f'smooth40_{channel}_reactions_per_true_C'])
        count = None
        if e in DD.SI1_E:
            count = float(DD.counts_per_C_table(m)[list(DD.SI1_E).index(e)])
        for power in (0, 1):
            q = float(row[f'observed_normalized_excess_k{power}'])
            out.append(dict(material=m, E_cm_keV=e, channel=channel, current_power=power,
                required_model_reactions_per_true_C_if_smooth40_neutron_reference=q*y,
                smooth40_reactions_per_true_C=y,
                SI1_detected_counts_per_cup_C=count,
                empirical_reactions_per_true_C_times_unknown_efficiency=None if count is None else
                    count*float(DD.collected_fraction(e,1.,0.)),
                note='absolute model demand assumes unchanged smooth40 neutron6.5 anchor; empirical Y=true cup_fraction*detected_counts/eta_p'))
    return out


if __name__=='__main__':
    result=rows()
    write_csv(OUT/'readout_floats.csv',result)
    print('Saved40 rows: conditional model reactions/true C, actual SI1 detected counts/cup C; no invented efficiency.')
