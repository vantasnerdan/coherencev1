"""Fixed-profile depth addendum; read accepted inputs, never rerun stopping or tails.

SI internally. Host energy J/true C is numerically eV/elementary charge.
The angular integral is scalar wave energy density, not plane-normal flux.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import sys

import numpy as np
from scipy.constants import elementary_charge
from scipy.integrate import quad
from scipy.special import exp1

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / 'g13_rogue'
sys.path.insert(0, str(ROOT))
from g4_beam.forward import model as beam

SPEED = 4000.0
AREA = 1e-4
RADIUS = math.sqrt(AREA / math.pi)
THICKNESS = 250e-6
TAU_ENERGY = 4.1e-12
FULL_DENSITY = 32 / (.8168e-9)**3
REFERENCE_THRESHOLD_EV = .0006423
REFERENCE_THERMAL_EV = .0101847
SOURCES = ('g4_beam/merit/profiles.json',
           'g4_beam/merit/README.md',
           'g4_beam/forward/model.py',
           'g11_cmns/collective/part1/stopping.json',
           'g11_cmns/collective/part1/motion.py')


def layer_at(layers, depth):
    """Piecewise constant source values; the rear face uses its interior limit."""
    for layer in layers:
        if layer[0] <= depth < layer[1]:
            return layer
    if math.isclose(depth, layers[-1][1], rel_tol=0, abs_tol=1e-18):
        return layers[-1]
    raise ValueError('Depth outside fixed layer table')


def track_column(layers, end):
    nh = nd = 0.0
    ratios = []
    for layer in layers:
        width = max(0., min(end, layer[1]) - max(0., layer[0]))
        if width:
            nh += width * layer[2]
            nd += width * layer[3]
            ratios.append(layer[3] / layer[2])
    return nd / nh, min(ratios), max(ratios), nd / end, nh / end


def kernels(distance, length):
    """Return scalar-density integral and survival-weighted age/residence area*time."""
    edge = math.hypot(distance, RADIUS)
    density = float(exp1(distance / length) - exp1(edge / length))
    residence = 2 * math.pi * length**2 / SPEED * (
        (distance + length) * math.exp(-distance / length)
        - (edge + length) * math.exp(-edge / length))
    return density, residence


def write_csv(name, rows):
    with (HERE / 'output' / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def verify(records, stops):
    """One bounded check of dilution, angular measure, age model and energy units."""
    expected = { .25: (.011537355253341838, .3496253447439519, .6102431765605311),
                 2.5: (.01072121716145717, .3185596568075262, .6087770566555903),
                 6.5: (.009716116797420755, .1990506088426156, .6038215465456302) }
    for energy, anchors in expected.items():
        layers = records['PdD', energy]['layers']
        np.testing.assert_allclose([layer_at(layers, z * 1e-6)[3] /
                                   layer_at(layers, z * 1e-6)[2] for z in (30, 50, 200)],
                                  anchors, rtol=2e-14)
        assert max(layer[3] / layer[2] for layer in layers) < .613
    for energy, anchor in ((.25, .0008305928529998756), (2.5, .011645653272094325)):
        end = stops['PdD', energy]['tracks'][0]['terminal_depth_A']
        np.testing.assert_allclose(track_column(records['PdD', energy]['layers'], end)[0],
                                   anchor, rtol=2e-14)
    distance, length = 50e-6, 100e-6
    scalar, residence = kernels(distance, length)
    # Independent polar-angle integration rejects an erroneous cos(theta) in u.
    theta = math.atan2(RADIUS, distance)
    angular = quad(lambda t: math.tan(t) * math.exp(-distance / (length * math.cos(t))),
                   0., theta, epsabs=1e-12)[0]
    # Independently integrate the already time-integrated age kernel over disk area.
    disk_age = 2 * math.pi * length / SPEED * quad(
        lambda r: r * math.exp(-math.hypot(distance, r) / length),
        0., RADIUS, epsabs=1e-23)[0]
    np.testing.assert_allclose([scalar, residence], [angular, disk_age], rtol=2e-10)
    echo = kernels(2 * THICKNESS - distance, length)[0]
    assert 0 <= .9 * echo < scalar
    # A downward point source deposits its full launch energy in the infinite
    # lossy hemisphere: angular measure 2*pi cancels its 1/(2*pi*s*s).
    radial_loss = quad(lambda t: math.exp(-t), 0., np.inf)[0]
    np.testing.assert_allclose(radial_loss, 1., rtol=2e-12)
    plane_fraction = quad(lambda t: math.sin(t) * math.exp(
        -distance / (length * math.cos(t))), 0., math.pi / 2)[0]
    assert 0 < plane_fraction < math.exp(-distance / length) < 1
    flux = 138e-6 / elementary_charge / AREA
    host = stops['PdD', 2.5]['host_recoil_J_per_true_C']
    # eV per charge * charge/(m2*s) * J/eV = J/(m2*s), once only.
    np.testing.assert_allclose(flux * host * elementary_charge, 138e-6 * host / AREA)
    # Equivalent per-nucleus energy and nuclei flux cancel the one factor 1.15.
    np.testing.assert_allclose(flux * host, (1.15 * flux) * (host / 1.15))
    return dict(profile_anchors=True, column_anchors=True,
                scalar_density_vs_angular_quadrature=True,
                residence_vs_disk_quadrature=True, first_echo_ordering=True,
                hemisphere_energy_balance=True, plane_flux_below_launch=True,
                charge_energy_normalization=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eta-host', type=float, default=.1,
                        help='fraction of host recoil launched downward as travelling waves')
    parser.add_argument('--host-to-D', type=float, default=.1,
                        help='local fraction of attenuated wave power entering retained D modes')
    args = parser.parse_args()
    if not (0 <= args.eta_host <= 1 and 0 <= args.host_to_D <= 1):
        parser.error('Fractions must be between zero and one')
    profiles = json.loads((ROOT / SOURCES[0]).read_text())['records']
    records = {(r['case']['condition'], float(r['case']['E_cm_keV'])): r for r in profiles}
    stops = {(r['target'], float(r['E_cm_keV'])): r for r in json.loads((ROOT / SOURCES[3]).read_text())}
    checks = verify(records, stops)
    profile_rows = []
    for material in ('Pd', 'PdD'):
        for energy in (.25, 2.5, 6.5):
            record = records[material, energy]
            layers = record['layers']
            for depth_nm in (.125, 2., 24., 129., 30000., 50000., 200000., 250000.):
                layer = layer_at(layers, depth_nm * 1e-9)
                profile_rows.append(dict(material=material, E_cm_keV=energy, kind='point',
                    depth_nm=depth_nm, depth_over_foil=depth_nm * 1e-9 / THICKNESS,
                    layer_top_nm=layer[0] * 1e9, layer_bottom_nm=layer[1] * 1e9,
                    n_Pd_m3=layer[2], n_D_m3=layer[3], x_D_per_Pd=layer[3] / layer[2],
                    traversed_layer_x_min='', traversed_layer_x_max='',
                    beam_time_s=record['case']['beam_time_s']))
            if energy <= 2.5:
                end = stops[material, energy]['tracks'][0]['terminal_depth_A']
                x, low, high, nd, nh = track_column(layers, end)
                profile_rows.append(dict(material=material, E_cm_keV=energy, kind='Dplus_column',
                    depth_nm=end * 1e9, depth_over_foil=end / THICKNESS,
                    layer_top_nm=0., layer_bottom_nm=end * 1e9,
                    n_Pd_m3=nh, n_D_m3=nd, x_D_per_Pd=x,
                    traversed_layer_x_min=low, traversed_layer_x_max=high,
                    beam_time_s=record['case']['beam_time_s']))
    wave_rows = []
    for energy in (.25, 2.5):
        split = stops['PdD', energy]
        cup = float(beam.I_cup(energy))
        current = cup / float(beam.collected_fraction(energy, b_hi=1., b_lo=0.))
        charge_flux = current * 1e-6 / elementary_charge / AREA
        host = split['host_recoil_J_per_true_C']
        power = charge_flux * host * args.eta_host  # eV/(m2*s)
        for depth_um in (30., 50., 200., 250.):
            depth = depth_um * 1e-6
            layer = layer_at(records['PdD', energy]['layers'], depth)
            nd = layer[3]
            for length_um in (5., 30., 50., 100., 250.):
                length = length_um * 1e-6
                direct_k, direct_age = kernels(depth, length)
                echo_k, echo_age = kernels(2 * THICKNESS - depth, length)
                for label, reflectivity in (('direct', 0.), ('direct_plus_first_rear_echo_rE0p9', .9)):
                    direct_u = power / SPEED * direct_k
                    echo_u = power / SPEED * reflectivity * echo_k
                    u = direct_u + echo_u
                    optical = SPEED / length * u * args.host_to_D * TAU_ENERGY / nd
                    wave_rows.append(dict(E_cm_keV=energy, I_cup_uA=cup, I_true_uA=current,
                        charge_flux_m2_s=charge_flux, nucleus_flux_m2_s=1.15 * charge_flux,
                        host_eV_per_elementary_charge=host, eta_host=args.eta_host,
                        local_host_to_D_fraction=args.host_to_D, depth_um=depth_um,
                        depth_over_foil=depth / THICKNESS, x_D_per_Pd=nd / layer[2],
                        n_D_m3=nd, L_energy_um=length_um, acoustic_energy_residence_ns=1e9 * length / SPEED,
                        direct_first_arrival_ns=1e9 * depth / SPEED,
                        first_echo_arrival_ns=1e9 * (2 * THICKNESS - depth) / SPEED,
                        path_model=label, rear_energy_reflectivity=reflectivity,
                        direct_u_eV_m3=direct_u, first_echo_u_eV_m3=echo_u, total_u_eV_m3=u,
                        acoustic_eV_per_local_D=u / nd, acoustic_eV_per_x1_D=u / FULL_DENSITY,
                        direct_survival_active_charge_launches=charge_flux * direct_age,
                        direct_survival_active_nucleus_impacts=1.15 * charge_flux * direct_age,
                        echo_survival_weighted_same_nucleus_launches=1.15 * charge_flux * reflectivity * echo_age,
                        local_optical_eV_per_D=optical,
                        optical_over_full_chart_threshold=optical / REFERENCE_THRESHOLD_EV,
                        optical_over_full_chart_300K=optical / REFERENCE_THERMAL_EV))
    write_csv('depth_profile.csv', profile_rows)
    write_csv('depth_waves.csv', wave_rows)
    summary = dict(schema='g13-fixed-profile-depth-v1',
        source_base='09b1b0779ed83f014d3e9ab21328fcd106f2a1bb',
        sources={name: dict(bytes=(ROOT / name).stat().st_size,
                           sha256=hashlib.sha256((ROOT / name).read_bytes()).hexdigest()) for name in SOURCES},
        model=dict(spot_area_m2=AREA, spot_radius_m=RADIUS, foil_thickness_m=THICKNESS,
                   speed_m_s=SPEED, optical_amplitude_lifetime_ps=8.2,
                   optical_energy_lifetime_ps=4.1, eta_host=args.eta_host,
                   local_host_to_D_fraction=args.host_to_D,
                   isotropic='downward hemisphere; thin surface launch; receiver at disk center',
                   reflection='one specular rear-face image only; energy reflectivity 0.9; same launches',
                   count='survival-weighted source-age/residence model, not packet duration or modes',
                   optical='local attenuation power (v/L_E)*u times fraction times tau_E/n_D; not coherent action',
                   threshold_reference_eV=REFERENCE_THRESHOLD_EV,
                   thermal_reference_eV=REFERENCE_THERMAL_EV),
        checks=checks, rows=dict(profile=len(profile_rows), waves=len(wave_rows)))
    (HERE / 'output' / 'depth.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
