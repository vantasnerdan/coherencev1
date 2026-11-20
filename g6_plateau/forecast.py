"""Fixed-profile, fixed-fit plateau counterfactuals at equal arrival energy.

Use the existing exact splitting and stopping quadrature. Only reaction kernels
and the explicitly named LiD replacement column differ from the reference.
Run from repo root: python -m g6_plateau.forecast
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from g4_beam.forward import model as beam
from g4_beam.forward.consumer import QE, csda_nodes, forward_yields
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.merit.api import ENERGIES, bestfits, reference
from g4_beam.transport import bca

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
A_LID = 4.0648  # Angstrom, conventional cell; Pandey & Stoneham (1985), Table 2
M_LI = 6.94  # atomic mass, natural isotopic mix
M_T = 3.0160493  # amu; keep the atomic-mass convention of the beam code
D_T_BG = 34.3827  # sqrt(keV), Bosch–Hale Table IV
D_T_NUM = [6.927e4, 7.454e8, 2.050e6, 5.2002e4]
D_T_DEN = [1., 6.38e1, -9.95e-1, 6.981e-5, 1.728e-4]
GAUSS_X, GAUSS_W = np.polynomial.hermite_e.hermegauss(5)
GAUSS_W /= GAUSS_W.sum()


def dt_stopping(energy, layer):
    """Exactly the retained LS electronic law; only fuel recoil mass changes."""
    _, _, _, density, _, _, _ = layer
    energy_eV = energy*1e3
    return transport_stopping(energy, layer) + density*1e-18*(
        bca.sn_zbl(energy_eV, 1, M_T) - bca.sn_zbl(energy_eV, 1, beam.M_D)
    )


def settings(cause, fit, *, hypothesis='known'):
    """Snapshot a fit without optimizing bulk U or altering the measured layer."""
    params = fit['params']
    u = params['Ue_PdD'] if hypothesis != 'no_screen' else 0.
    candidate = {'screeningrise': 'sigmoid_tt', 'surface': 'layer',
                 'neutrals': 'fast_neutrals', 'lattice': 'none', 'resonance': 'none'}[cause]
    extra = {}
    if candidate == 'sigmoid_tt':
        extra.update(Ue_s=params['Ues_PdD'] if hypothesis != 'no_screen' else 0., E0=params['E0'])
    elif candidate == 'layer':
        extra.update(Ue_s=params['Ues_PdD'] if hypothesis != 'no_screen' else 0., d_nm=params['d'])
    elif candidate == 'fast_neutrals':
        extra.update(f_n=params['f_n'], E_n_lab_keV=params['E_n'], delta=int(fit['cand'][-1]))
    bh, bl = fit.get('b_hi', 0), fit.get('b_lo', 0)
    current = lambda e: beam.I_cup(e)*1e-6/beam.collected_fraction(e, bh, bl)
    measurement = dict(current_A=current, neutral_reference_current_A=beam.I_INT*1e-6)
    return u, candidate, extra, measurement


def cross_section(lab, shift, isotope, reaction, channel=None):
    """Screened Eq. 7, in mb. No D–D resonance is transferred to D–T."""
    if isotope == 'DT':
        ec = np.asarray(lab)*M_T/(M_T+beam.M_D)
        s = np.polynomial.polynomial.polyval(ec, D_T_NUM)/np.polynomial.polynomial.polyval(ec, D_T_DEN)
        return s/np.sqrt(ec*(ec+shift))*np.exp(-D_T_BG/np.sqrt(ec+shift))
    if isotope != 'DD':
        raise ValueError(f'Unsupported isotope pair: {isotope}')
    ec = np.maximum(np.asarray(lab)/2, .05)  # exact forward-consumer convention
    if channel is None:
        return beam.reaction_sigma(ec, shift, 'p', reaction)+beam.reaction_sigma(ec, shift, 'n', reaction)
    return beam.reaction_sigma(ec, shift, channel, reaction)


def projected_yield(e, column, u, candidate, params, measurement, *, isotope='DD',
                    reaction='proton', channel=None, lattice_shift=None):
    """True reactions/C. e is the nominal DD CM axis, not DT CM or lab."""
    layers, nads, _ = column(e)
    stop = dt_stopping if isotope == 'DT' else transport_stopping
    share2 = measurement.get('D2_fraction', .15)
    spread = measurement.get('spread_FWHM', .05)/2.3548
    boundary = params['d_nm']*1e-9 if candidate == 'layer' else None

    def ion(e0):
        def xs(energy, depth):
            shift = u
            if boundary is not None:
                shift = np.where(depth < boundary, params['Ue_s'], u)
            elif candidate == 'sigmoid_tt':
                # Drive coordinate stays the original fitted *incident lab* energy/2.
                shift += (params['Ue_s']-u)/(1+np.exp((e0/2-params['E0'])/.5))
            if lattice_shift is not None:
                shift = shift+lattice_shift(energy)*1e-3  # eV -> keV
            return cross_section(energy, shift, isotope, reaction, channel)

        en, z, weights, _ = csda_nodes(e0, layers, stop, break_m=boundary)
        return nads*float(xs(e0, 0.))*1e-31 + float(weights@xs(en, z))

    ions = sum(w*((1-share2)*ion(2*e*(1+spread*x))+2*share2*ion(e*(1+spread*x)))
               for x, w in zip(GAUSS_X, GAUSS_W))
    neutrals = 0.
    if candidate == 'fast_neutrals':
        ratio = params['f_n']
        if params['delta']:
            ratio *= measurement['neutral_reference_current_A']/measurement['current_A'](e)
        neutrals = ratio*ion(params['E_n_lab_keV'])
    return (ions+neutrals)/QE


def lid_column(_):
    density = 4/(A_LID*1e-10)**3
    return [(0., 1e-5, density, density, 3, M_LI, density)], 0., {}


def uniform_bulk_fit():
    data = json.loads((ROOT/'g4_beam/forward/fits.json').read_text())
    return min((fit for fit in data if fit['cand']=='bulk 2.5-6.5'), key=lambda fit: fit['chi2'])


def forecast(lattice, collision, out):
    fits = bestfits()
    bulk = uniform_bulk_fit()
    near = collision['contrast_by_R'][0]
    sensitivity = {r['E_lab_Dplus_keV']: r['hypothetical_yield_ratio']
                   for r in collision['counterfactual']['yield_factors']}
    rows = []
    parity = []
    column = reference().profiles
    for cause, fit in [*fits.items(), ('lattice', bulk), ('resonance', bulk)]:
        u, candidate, params, meas = settings(cause, fit)
        _, _, zeros, _ = settings(cause, fit, hypothesis='no_screen')
        shift_eV = float(lattice['susceptibility_locked_minus_unlocked_eV'])
        extra = (lambda _energy: shift_eV) if cause == 'lattice' else None
        for ee in ENERGIES:
            e = float(ee)
            args = (e, column, u, candidate, params, meas)
            dd = projected_yield(*args, lattice_shift=extra)
            dt = projected_yield(*args, isotope='DT', lattice_shift=extra)
            dd_bh = projected_yield(*args, reaction='eq7', lattice_shift=extra)
            # Two independent LiD questions. Never imply a measured LiD U or sheet.
            li_transfer = projected_yield(e, lid_column, u, candidate, params, meas)
            li_unscreened = projected_yield(e, lid_column, 0., candidate, zeros, meas)
            snapshot = cause == 'lattice'
            rows.append(dict(cause=cause, E_lab_Dplus_keV=2*e, E_cm_nominal_DD_keV=e,
                             u_PdD_bulk_keV=u, source_fit_chi2=fit['chi2'], source_fit_meas=fit['meas'],
                             PdD_DD_events_per_true_C=dd, PdT_DT_events_per_true_C=dt,
                             DT_over_DD_same_lab=dt/dd,
                             LiD_over_PdD_transfer_U_no_LiD_sheet=li_transfer/dd,
                             LiD_over_PdD_zero_U_no_LiD_sheet=li_unscreened/dd,
                             DD_resonance_over_Bosch_Hale=dd/dd_bh,
                             lattice_collision_snapshot=('g6_couple/collision_projection.json' if snapshot else None),
                             lattice_snapshot_BO_lowering_change_eV_at_0p1A=(
                                 near['delta_BO_lowering_eV'] if snapshot else None),
                             lattice_snapshot_encounter_work_change_eV_at_0p1A=(
                                 near['delta_encounter_work_eV'] if snapshot else None),
                             lattice_snapshot_hypothetical_rigid_shift_DDp_yield_ratio=(
                                 sensitivity.get(2*e) if snapshot else None)))
            if cause in fits:
                ref = 0.
                for ch in ('p', 'n'):
                    result = forward_yields([e], column, Ue_keV=u, candidate=candidate, parameters=params,
                                            measurement={**meas, 'channel': ch}, stopping_fn=transport_stopping)
                    ref += float(np.asarray(result['physical_yield_per_C'])[0])
                parity.append(abs(dd/ref-1))
    worst = max(parity)
    if worst > 2e-12:
        raise ArithmeticError(f'Retained DD consumer disagrees ({worst:.4g} relative)')
    write_csv(out, rows)
    return {'DD_production_parity_max_relative': worst, 'rows': len(rows)}


def write_csv(out, rows):
    with out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def discriminators(lattice):
    with (ROOT/'g4_beam/causes/laboratory_ratios.csv').open() as f:
        rows = list(csv.DictReader(f))
    interventions = {'near_surface_D_factor': 'near_surface_D_x50',
                     'D2_fraction': 'D2_current_0.85', 'current': 'half_current_frozen'}
    projections = {(r['intervention'], r['E_cm_keV']): r for r in lattice['projections']}
    for row in rows:
        e = float(row['E_cm_keV'])
        row['E_lab_Dplus_keV'] = 2*e
        act = projections[interventions[row['knob']], e]
        base = projections['baseline', e]
        # Total response to changing beam or target vs. the incremental
        # locked/unlocked contrast at each fixed intervention.
        row.update(lattice_total_ratio_act_over_base_at_reference_Ue=(
                       act['baseline_physical_yield_per_C']/base['baseline_physical_yield_per_C']),
                   lattice_reference_Ue_keV=act['base_Ue_eV']*1e-3,
                   lattice_state_ratio_locked_over_unlocked_in_act=act['signed_linear_locked_over_unlocked'],
                   lattice_positive_frozen_contact_envelope_in_act=act['frozen_positive_envelope_factor'],
                   lattice_full_ionic_addition_in_act_not_contrast=act['full_ionic_addition_factor'],
                   lattice_regime='g6_couple frozen target/pump except stated profile or beam switch')
    return rows


def isotope_discriminator(forecasts):
    lab_keV = .5
    causes = ('screeningrise', 'surface', 'neutrals')
    with forecasts.open() as f:
        selected = {r['cause']: r for r in csv.DictReader(f)
                    if float(r['E_lab_Dplus_keV']) == lab_keV}
    ratios = {f'{c}_ratio': float(selected[c]['DT_over_DD_same_lab']) for c in causes}
    lattice = selected['lattice']
    return dict(knob='target_tritium_fraction', value=1.,
                regime='same-lab-energy counterfactual; replace D by T in frozen PdD profile and sheet',
                E_cm_keV=lab_keV/2, E_lab_Dplus_keV=lab_keV,
                separation_factor=max(ratios.values())/min(ratios.values()), **ratios,
                observable='total D-T / D-D fusion events per true charged C at the same lab energy',
                lattice_total_ratio_act_over_base_at_reference_Ue=float(lattice['DT_over_DD_same_lab']),
                lattice_reference_Ue_keV=float(lattice['u_PdD_bulk_keV']),
                lattice_regime='fixed bulk PdD shift; no transferred lattice-state electronic response')


def main():
    lattice = json.loads((ROOT/'g6_couple/response.json').read_text())
    collision = json.loads((ROOT/'g6_couple/collision_projection.json').read_text())
    forecasts = HERE/'forecasts.csv'
    status = forecast(lattice, collision, forecasts)
    write_csv(HERE/'discriminators.csv', [*discriminators(lattice), isotope_discriminator(forecasts)])
    print(json.dumps(status, indent=2))


if __name__ == '__main__':
    main()
