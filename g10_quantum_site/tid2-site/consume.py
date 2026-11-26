"""One scoped G8 numerical consumer and G9 signed-table mapping, no envelope runs."""
import hashlib
import importlib.util
import itertools
import json
import sys

import numpy as np

from solve import OUT, dump, table
from well import HERE, ROOT, MASSES, Well


def source_evidence():
    zero = json.loads((HERE/'inputs/zero.json').read_text())
    position = np.asarray(zero['positions_A']); cell = np.asarray(zero['cell_A'])
    n = np.ones(3)/np.sqrt(3)
    titanium = []
    for index, element in enumerate(zero['physical_species']):
        if element != 'Ti':
            continue
        for image in itertools.product((-1, 0, 1), repeat=3):
            delta = position[index]+np.asarray(image)@cell-position[8]
            titanium.append((float(np.linalg.norm(delta)), float(delta@n), delta.tolist()))
    nearest_distance = min(row[0] for row in titanium)
    nearest = [row for row in titanium if abs(row[0]-nearest_distance) < 1e-8]
    assert len(nearest) == 4
    assert sum(row[1] < 0 for row in nearest) == 1 and sum(row[1] > 0 for row in nearest) == 3
    cuts = json.loads((HERE/'inputs/cuts.json').read_text())
    source_rows = []
    for direction, groups in cuts['directions'].items():
        ray = np.array([int(v) for v in direction]); ray = ray/np.linalg.norm(ray)
        for group in groups:
            for point in group['signed_points']:
                u = point['displacement_A']; case = direction+('m' if u < 0 else 'p')+f'{round(abs(u)*1000):03d}'
                path = HERE/f'inputs/{case}.json'; data = path.read_bytes(); native = json.loads(data)
                assert hashlib.sha256(data).hexdigest() == point['source']['sha256']
                actual_displacement = np.array(native['positions_A'])[8]-position[8]
                assert np.allclose(actual_displacement, u*ray, atol=1e-12)
                energy = native['energy_free_eV']-zero['energy_free_eV']
                force = float(np.asarray(native['forces_eV_A'])[8]@ray)
                assert abs(energy-point['delta_energy_free_eV']) < 1e-12
                assert abs(force-point['target_force_parallel_eV_A']) < 1e-12
                source_rows.append(dict(case=case, path=str(path.relative_to(ROOT)), sha256=hashlib.sha256(data).hexdigest(),
                                        energy_eV=energy, force_along_positive_ray_eV_A=force))
    anchors = Well().anchors()
    symmetry = {key: dict(max_signed_energy_residual_eV=max(abs(g['twofold_symmetry_diagnostic']['plus_minus_energy_free_eV']) for g in cuts['directions'][key]),
                          max_rotated_force_residual_eV_A=max(g['twofold_symmetry_diagnostic']['maximum_rotated_force_mismatch_eV_A'] for g in cuts['directions'][key])) for key in ('100', '110')}
    output = dict(accepted_base='463910dcf13ce69d74b138a21d53183129a779de', source_atom_index=8,
                  axes='conventional cubic x100 y010 z001', source_atom_A=position[8].tolist(),
                  nearest_Ti_distance_A=nearest_distance, nearest_Ti_vectors_and_signed_projections=nearest,
                  signed_111='+ moves away from the single Ti at -[111], toward the other three; - moves toward the single Ti',
                  symmetry=symmetry, points=source_rows,
                  max_anchor_energy_residual_eV=max(abs(row['model_eV']-row['source_eV']) for row in anchors),
                  max_anchor_slope_residual_eV_A=max(abs(row['model_slope_eV_A']-row['source_slope_eV_A']) for row in anchors),
                  mesh=json.loads((HERE/'inputs/setting.json').read_text()),
                  scope='native bytes, coordinate sign and interpolation checked; no electronic runs replayed')
    dump(OUT/'source_evidence.json', output)
    return output


def consumers():
    sys.path.insert(0, str(ROOT/'g8_why/collision/source'))
    from bound_forward import BoundSpectatorKernel
    all_rows = []
    for isotope in MASSES:
        for temperature in (77, 300):
            path = OUT/f'{isotope}_{temperature}K_g8_joint.json'
            joint = json.loads(path.read_text()); states = joint['states']
            total = sum(sum(state['probabilities']) for state in states)
            kinetic = sum(sum(p*k for p, k in zip(state['probabilities'], state['kinetic_eV'])) for state in states)
            assert abs(total-1) < 1e-12
            assert all(np.isfinite(state['kinetic_eV']).all() and min(state['probabilities']) >= 0 for state in states)
            all_rows.append(dict(isotope=isotope, temperature_K=temperature, joint_normalization=total,
                                 mean_kinetic_eV=kinetic, maximum_excitation_eV=max(state['excitation_eV'] for state in states),
                                 maximum_kinetic_node_eV=max(max(state['kinetic_eV']) for state in states),
                                 state_file=str(path.relative_to(ROOT))))
    states = json.loads((OUT/'D_300K_g8_joint.json').read_text())['states']
    kernel = BoundSpectatorKernel(ROOT, angular_order=32)
    value = kernel.fold('TiD', [.25], states, removal_eV=0.)
    assert np.isfinite(value).all() and value[0] > 0
    # The numerical interface is exercised, but the pinned TiD transport profile
    # is NOT TiD2 fluorite or an assigned foil. No reported yield/chart comparison.
    g8 = dict(reader='g8_why/collision/source/bound_forward.py:BoundSpectatorKernel.fold',
              state_measure='D_300K_g8_joint.json', parser_model_exercised=True,
              interface_profile='existing TiD ONLY, not a new TiD2 profile', nominal_cm_keV=.25,
              finite_positive_fold=True, result_scope='numerical schema/interface smoke only; not a TiD2 rate',
              isotope_scope='H/T retain their own Fourier/kinetic masses; equal-mass D-D consumer executed only for D',
              rows=all_rows,
              spectator_envelope='The retained TiD2 conditional excitation nodes are far below the existing G8 spectator envelope excitation cap of 100 eV; the envelope is distribution-free in momentum, so this site table supplies a shape inside that assumption, not a new envelope or Ti foil yield.')
    # G9 has no density-CSV API: consume the real module's coordinate constants,
    # import the table, and exercise exactly the existing signed ballistic map.
    spec = importlib.util.spec_from_file_location('tid2_g9_edge', ROOT/'g9_edge/edge.py')
    edge = importlib.util.module_from_spec(spec); spec.loader.exec_module(edge)
    plus = np.genfromtxt(OUT/'D_300K_g9_plus111.csv', delimiter=',', names=True)
    minus = np.genfromtxt(OUT/'D_300K_g9_minus111.csv', delimiter=',', names=True)
    assert np.allclose(plus['s_A'], -minus['s_A'][::-1])
    assert np.allclose(plus['marginal_A_inverse'], minus['marginal_A_inverse'][::-1])
    Z0 = float(edge.MATCH); v = float(np.sqrt(4*1000*.25/edge.MD))
    s = plus['s_A']; z = Z0+s; shift = s/v
    table(OUT/'g9_mapping.csv', ['s_A', 'marginal_A_inverse', 'Z0_A', 'z_target_A', 'speed_A_fs', 'arrival_shift_fs', 'opposite_normal_z_target_A', 'opposite_normal_arrival_shift_fs'],
          [s, plus['marginal_A_inverse'], s*0+Z0, z, s*0+v, shift, Z0-s, -shift])
    g9 = dict(source='g9_edge/edge.py:History, MATCH and MD',
              convention='registered n=+[111]/sqrt(3); s=displacement dot n; z_target=Z0+s; arrival_shift=s/v',
              opposite_normal='n -> -n gives s -> -s, reflected density, z_target=Z0-s and arrival_shift=-s/v for the same crystal displacement',
              normalized_plane_marginal=float(np.trapezoid(plus['marginal_A_inverse'], s)),
              speed_A_fs=v, registered_example_Z0_A=Z0,
              scope='actual table import and applicable G9 ballistic coordinate semantics; no crystallographic density API, slab target or Ti foil surface assigned',
              distinct_observables='marginal_A_inverse is an integral over the entire transverse plane; on_axis_line_density_A_minus3 is only rho(s*n), not normalized as a 1D probability')
    bindings = []
    for name in ['g8_why/collision/source/bound_forward.py', 'g8_why/collision/source/yield_kernel.py', 'g8_why/common.py', 'g8_why/collision/evidence/envelope.json', 'g9_edge/edge.py', 'g13_rogue/run.py']:
        data = (ROOT/name).read_bytes()
        bindings.append(dict(path=name, sha256=hashlib.sha256(data).hexdigest(), bytes=len(data)))
    result = dict(g8=g8, g9=g9, external_accepted_base_bindings=bindings)
    dump(OUT/'consumer_evidence.json', result)
    return result


if __name__ == '__main__':
    source_evidence()
    consumers()
