"""Finite oscillator diagonalization, normalized position/Fourier observables."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.constants import hbar
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.ndimage import map_coordinates, spline_filter
from scipy.special import roots_hermite

from well import HERE, ROOT, KB, MASSES, Well, kinetic

OUT = HERE / 'output'
BLOCKS = [((0, 0, 0), (1, 1, 1), 1), ((1, 0, 0), (0, 1, 1), 3)]


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def table(path, headers, columns):
    np.savetxt(path, np.column_stack(columns), delimiter=',', header=','.join(headers), comments='', fmt='%.13g')


def hermite(x, n, length, momentum=False):
    q = np.asarray(x)*(length if momentum else 1/length)
    a = np.empty((len(q), n), complex if momentum else float)
    a[:, 0] = np.pi**(-.25)*np.exp(-q*q/2)
    if n > 1:
        a[:, 1] = np.sqrt(2)*q*a[:, 0]
    for j in range(2, n):
        a[:, j] = np.sqrt(2/j)*q*a[:, j-1]-np.sqrt((j-1)/j)*a[:, j-2]
    a *= np.sqrt(length) if momentum else 1/np.sqrt(length)
    if momentum:
        a *= (-1j)**np.arange(n)
    return a


def values(v, maps):
    return np.einsum('ai,bj,ck,ijks->abcs', *maps, v, optimize=True)


def project_potential(V, maps, p, q):
    pairs = [(maps[:, i::2, None]*maps[:, None, j::2]).reshape(len(maps), -1) for i, j in zip(p, q)]
    raw = np.einsum('abc,ai,bj,ck->ijk', V, *pairs, optimize=True)
    m = maps.shape[1]//2
    return raw.reshape((m,)*6).transpose(0, 2, 4, 1, 3, 5).reshape(m**3, m**3)


def solve(isotope, n=20, order=40, variant='base'):
    mass = MASSES[isotope]; c = kinetic(mass)
    length = (2*c/5.58)**.25
    nodes, weights = roots_hermite(order)
    grid = nodes*length
    maps = hermite(grid, n, length)*np.sqrt(weights*np.exp(nodes*nodes)*length)[:, None]
    xyz = np.stack(np.broadcast_arrays(grid[:, None, None], grid[None, :, None], grid[None, None, :]), axis=-1)
    V = Well(variant)(xyz)
    ladder = np.diag(np.sqrt(np.arange(1, n+2)), 1)
    deriv = (ladder-ladder.T)/(np.sqrt(2)*length)
    T = (-c*deriv@deriv)[:n, :n]
    m = n//2; I = np.eye(m)
    spectra = []
    for p, q, mult in BLOCKS:
        diagonal = []
        for par in (p, q):
            ts = [T[j::2, j::2] for j in par]
            H = project_potential(V, maps, par, par)
            H += np.kron(np.kron(ts[0], I), I)+np.kron(np.kron(I, ts[1]), I)+np.kron(np.kron(I, I), ts[2])
            diagonal.append(H)
        cross = project_potential(V, maps, p, q)
        H = np.block([[diagonal[0], cross], [cross.T, diagonal[1]]])
        energies, vectors = eigh(H, check_finite=False, overwrite_a=True)
        spectra.append((energies, vectors))
    return dict(n=n, order=order, length=length, mass=mass, variant=variant, spectra=spectra)


def save_spectrum(path, spec):
    np.savez_compressed(path, n=spec['n'], order=spec['order'], length=spec['length'], mass=spec['mass'],
                        variant=spec['variant'], E0=spec['spectra'][0][0], C0=spec['spectra'][0][1],
                        E1=spec['spectra'][1][0], C1=spec['spectra'][1][1])


def load_spectrum(path):
    z = np.load(path)
    return dict(n=int(z['n']), order=int(z['order']), length=float(z['length']), mass=float(z['mass']),
                variant=str(z['variant']), spectra=[(z['E0'], z['C0']), (z['E1'], z['C1'])])


def thermal(spec, temperature):
    minimum = min(e[0] for e, _ in spec['spectra'])
    Z = sum(b[2]*np.exp(-(e-minimum)/(KB*temperature)).sum() for b, (e, _) in zip(BLOCKS, spec['spectra']))
    for b, (e, v) in zip(BLOCKS, spec['spectra']):
        weights = np.exp(-(e-minimum)/(KB*temperature))/Z
        selected = weights > 1e-13
        yield b, e[selected]-minimum, v[:, selected], weights[selected]


def density(spec, temperature, grid, momentum=False):
    maps = hermite(grid, spec['n'], spec['length'], momentum)
    rho = np.zeros((len(grid),)*3)
    retained = 0.
    for (p, q, mult), energies, v, w in thermal(spec, temperature):
        representative = np.zeros_like(rho)
        for start in range(0, len(w), 4):
            wave = None
            for i, par in enumerate((p, q)):
                block = v[i*(spec['n']//2)**3:(i+1)*(spec['n']//2)**3, start:start+4]
                val = values(block.reshape((spec['n']//2,)*3+(-1,)), [maps[:, j::2] for j in par])
                wave = val if wave is None else wave+val
            representative += np.einsum('ijks,s->ijk', abs(wave)**2, w[start:start+4])
        rho += representative
        if mult == 3:
            rho += representative.transpose(1, 2, 0)+representative.transpose(2, 0, 1)
        retained += mult*w.sum()
    return rho, float(1-retained)


def interpolator(grid, rho):
    coefficients = spline_filter(rho, order=3)
    def sample(points):
        points = np.asarray(points)
        coordinates = ((points.reshape(-1, 3)-grid[0])/(grid[1]-grid[0])).T
        return map_coordinates(coefficients, coordinates, order=3, prefilter=False, mode='constant', cval=0.).reshape(points.shape[:-1])
    return sample


def radial(grid, rho, radii, angular_order=20):
    mu, weights = np.polynomial.legendre.leggauss(angular_order)
    phi = np.arange(2*angular_order)*np.pi/angular_order
    directions = np.stack(np.broadcast_arrays(np.sqrt(1-mu[:, None]**2)*np.cos(phi),
                           np.sqrt(1-mu[:, None]**2)*np.sin(phi), mu[:, None]+phi*0), axis=-1).reshape(-1, 3)
    angular_weights = np.repeat(weights, 2*angular_order)/(4*angular_order)
    sample = interpolator(grid, rho)
    d = sample(radii[:, None, None]*directions[None, :, :])@angular_weights
    measure = 4*np.pi*radii*radii*d
    tail = -cumulative_trapezoid(measure[::-1], radii[::-1], initial=0)[::-1]
    return dict(coordinate=radii, density=d, measure=measure, tail=tail)


def cutoff(x, y):
    peak = int(np.argmax(y)); ids = np.flatnonzero((np.arange(len(x)) > peak) & (y <= y[peak]*1e-6))
    if not len(ids):
        return None
    i = int(ids[0])
    return float(x[i-1]+(x[i]-x[i-1])*(y[peak]*1e-6-y[i-1])/(y[i]-y[i-1]))


def moment(grid, rho, power):
    dx = grid[1]-grid[0]
    return float(np.einsum('i,ijk->', grid**power, rho)*dx**3)


def observables(spec, temperature, fine=False):
    x = np.linspace(-1.2, 1.2, 141 if fine else 101)
    k = np.linspace(-52, 52, 177 if fine else 129)
    rho, omitted = density(spec, temperature, x)
    mom, _ = density(spec, temperature, k, True)
    radial_x = radial(x, rho, np.linspace(0, 1.2, 301), 28 if fine else 20)
    radial_k = radial(k, mom, np.linspace(0, 52, 301), 28 if fine else 20)
    summary = dict(thermal_weight_omitted=omitted, position_normalization=moment(x, rho, 0), momentum_normalization=moment(k, mom, 0),
                   mean_x_A=moment(x, rho, 1), mean_x2_A2=moment(x, rho, 2), mean_k2_axis_A_inverse2=moment(k, mom, 2),
                   mean_kinetic_eV=3*kinetic(spec['mass'])*moment(k, mom, 2),
                   position_density_1e6_A=cutoff(radial_x['coordinate'], radial_x['density']),
                   position_measure_1e6_A=cutoff(radial_x['coordinate'], radial_x['measure']),
                   momentum_density_1e6_A_inverse=cutoff(radial_k['coordinate'], radial_k['density']),
                   momentum_measure_1e6_A_inverse=cutoff(radial_k['coordinate'], radial_k['measure']),
                   probability_r_gt_080=float(np.interp(.8, radial_x['coordinate'], radial_x['tail'])),
                   radial_position_normalization=float(radial_x['tail'][0]), radial_momentum_normalization=float(radial_k['tail'][0]),
                   momentum_inversion_defect_over_peak=float(abs(mom-mom[::-1, ::-1, ::-1]).max()/mom.max()))
    for name, dist in [('position', radial_x), ('momentum', radial_k)]:
        for kind in ('density', 'measure'):
            at = cutoff(dist['coordinate'], dist[kind])
            summary[f'{name}_{kind}_tail_at_1e6'] = None if at is None else float(np.interp(at, dist['coordinate'], dist['tail']))
    sample = interpolator(x, rho)
    s = np.linspace(-1.2, 1.2, 301)
    n = np.ones(3)/np.sqrt(3); a = np.array([1., -1., 0.])/np.sqrt(2); b = np.cross(n, a)
    transverse = np.linspace(-1.2, 1.2, 161 if fine else 121)
    aa, bb = np.meshgrid(transverse, transverse, indexing='ij')
    projected = []
    for value in s:
        plane = sample(value*n+aa[..., None]*a+bb[..., None]*b)
        projected.append(np.trapezoid(np.trapezoid(plane, transverse, axis=1), transverse))
    projected = np.asarray(projected)
    line = sample(s[:, None]*n)
    summary['projected_111_normalization'] = float(np.trapezoid(projected, s))
    summary['projected_111_positive_probability'] = float(np.trapezoid(projected[s >= 0], s[s >= 0]))
    summary['projected_111_mean_A'] = float(np.trapezoid(s*projected, s))
    summary['projected_111_third_moment_A3'] = float(np.trapezoid(s**3*projected, s))
    summary['axis_MSD_A2'] = [moment(x, rho, 2), float(np.einsum('j,ijk->', x*x, rho)*(x[1]-x[0])**3), float(np.einsum('k,ijk->', x*x, rho)*(x[1]-x[0])**3)]
    summary['axis_anisotropy_max_over_min'] = max(summary['axis_MSD_A2'])/min(summary['axis_MSD_A2'])
    return summary, radial_x, radial_k, (s, projected, line)


def joint_states(spec, temperature, bins=80):
    order = max(40, spec['n']+12)
    q, weights = roots_hermite(order)
    grid = q/spec['length']
    maps = hermite(grid, spec['n'], spec['length'], True)
    quadrature = weights*np.exp(q*q)/spec['length']
    measure = quadrature[:, None, None]*quadrature[None, :, None]*quadrature[None, None, :]
    K = kinetic(spec['mass'])/spec['length']**2*(q[:, None, None]**2+q[None, :, None]**2+q[None, None, :]**2)
    # Quadratic bins retain fine low-energy resolution and the full node support.
    edges = np.linspace(0, np.sqrt(K.max()+1e-12), bins+1)**2
    binid = np.minimum(np.searchsorted(edges, K.ravel(), side='right')-1, bins-1)
    states = []; raw_total = 0.; raw_kinetic = 0.
    for (p, qpar, mult), en, v, w in thermal(spec, temperature):
        for j in range(len(w)):
            wave = None
            for blockid, par in enumerate((p, qpar)):
                block = v[blockid*(spec['n']//2)**3:(blockid+1)*(spec['n']//2)**3, j:j+1]
                val = values(block.reshape((spec['n']//2,)*3+(1,)), [maps[:, k::2] for k in par])[..., 0]
                wave = val if wave is None else wave+val
            prob = (abs(wave)**2*measure).ravel()
            probability = np.bincount(binid, weights=prob, minlength=bins)
            energy = np.bincount(binid, weights=prob*K.ravel(), minlength=bins)
            active = probability > 1e-16
            joint = mult*w[j]*probability[active]
            raw_total += joint.sum(); raw_kinetic += mult*w[j]*energy[active].sum()
            states.append(dict(excitation_eV=float(en[j]), kinetic_eV=(energy[active]/probability[active]).tolist(), probabilities=joint.tolist()))
    for state in states:
        state['probabilities'] = (np.array(state['probabilities'])/raw_total).tolist()
    return dict(schema='g8-bound-spectral-states', states=states,
                metadata=dict(angular_assumption='isotropized conditional kinetic measure of each real 3D eigenstate; no momentum-excitation independence assumption',
                              p_equals_hbar_k=True, mass_u=spec['mass'], temperature_K=temperature, quadrature_order=order,
                              probability_before_normalization=float(raw_total), mean_kinetic_eV=float(raw_kinetic/raw_total),
                              kinetic_compression='80 quadratic bins; each retains probability and conditional kinetic mean; all Gauss-Hermite support retained',
                              scope='state measure only; no TiD2 foil yield or removal energy'))


def levels(spec):
    rows = []
    e0 = min(e[0] for e, _ in spec['spectra'])
    for block, (energies, _) in zip(BLOCKS, spec['spectra']):
        for i, energy in enumerate(energies[:24]):
            rows.append(dict(energy_eV=float(energy), excitation_eV=float(energy-e0), permutation_multiplicity=block[2], block=str(block[:2]), block_index=i))
    rows.sort(key=lambda row: row['energy_eV'])
    # Degeneracy clusters are observed numerically, not assigned inversion labels.
    clusters = []
    for row in rows:
        if not clusters or row['energy_eV']-clusters[-1]['energy_eV'] > 1e-6:
            clusters.append(dict(energy_eV=row['energy_eV'], excitation_eV=row['excitation_eV'], multiplicity=row['permutation_multiplicity']))
        else:
            clusters[-1]['multiplicity'] += row['permutation_multiplicity']
    for i, row in enumerate(clusters[:-1]):
        row['next_cluster_gap_eV'] = clusters[i+1]['energy_eV']-row['energy_eV']
    return clusters[:14]


def transitions():
    well = Well(); output = []
    for isotope, mass in MASSES.items():
        pair = []
        for dx in (.002, .001):
            x = np.arange(-2+dx, 2, dx)
            c = kinetic(mass)
            potential = well(np.stack((x, x*0, x*0), axis=-1))
            en, vectors = eigh_tridiagonal(2*c/dx**2+potential, np.full(len(x)-1, -c/dx**2), select='i', select_range=(0, 7))
            pair.append(dict(dx_A=dx, levels_eV=en.tolist(), first_gap_eV=float(en[1]-en[0]), second_gap_eV=float(en[2]-en[1]),
                             shift_per_quantum_eV=float(en[2]-2*en[1]+en[0]), probability_abs_u_gt_080=(vectors[abs(x) > .8]**2).sum(axis=0).tolist()))
        output.append(dict(isotope=isotope, coordinate='[100] signed one-dimensional frozen-neighbour coordinate, not 3D gaps',
                           mass_u=mass, result=pair[1], spacing_control=pair[0], first_gap_control_delta_eV=pair[1]['first_gap_eV']-pair[0]['first_gap_eV'],
                           second_gap_control_delta_eV=pair[1]['second_gap_eV']-pair[0]['second_gap_eV']))
    dump(OUT/'transitions_100.json', output)
    return output


def produce():
    OUT.mkdir(exist_ok=True)
    rows = []
    for isotope in MASSES:
        spec = solve(isotope)
        save_spectrum(OUT/f'{isotope}_spectrum.npz', spec)
        for temperature in (77, 300):
            name = f'{isotope}_{temperature}K'
            summary, position, momentum, projection = observables(spec, temperature)
            summary.update(isotope=isotope, mass_u=MASSES[isotope], temperature_K=temperature, n=spec['n'], quadrature_order=spec['order'], levels=levels(spec))
            for label, data in [('position', position), ('momentum', momentum)]:
                unit = ['r_A', 'density_A_minus3', 'radial_measure_A_inverse', 'tail_probability'] if label == 'position' else ['k_A_inverse', 'density_A3', 'radial_measure_A', 'tail_probability']
                columns = list(data.values())
                if label == 'momentum':
                    # p table uses atomic units: p/(atomic_mass * Angstrom/fs), with explicit SI too.
                    p_si = hbar*data['coordinate']*1e10
                    unit += ['p_kg_m_s', 'density_p_SI_inverse3', 'radial_measure_p_SI_inverse']
                    jac = hbar*1e10
                    columns += [p_si, data['density']/jac**3, data['measure']/jac]
                table(OUT/f'{name}_{label}.csv', unit, columns)
            s, marginal, line = projection
            table(OUT/f'{name}_g9_plus111.csv', ['s_A', 'marginal_A_inverse', 'on_axis_line_density_A_minus3', 'beyond_source_radius'], [s, marginal, line, abs(s) > .8])
            table(OUT/f'{name}_g9_minus111.csv', ['s_A', 'marginal_A_inverse', 'on_axis_line_density_A_minus3', 'beyond_source_radius'], [s, marginal[::-1], line[::-1], abs(s) > .8])
            dump(OUT/f'{name}_g8_joint.json', joint_states(spec, temperature))
            rows.append(summary)
            print(name, summary['position_density_1e6_A'], summary['momentum_density_1e6_A_inverse'], flush=True)
    dump(OUT/'single_site.json', rows)
    transitions()
    dump(OUT/'anchor_residuals.json', Well().anchors())


def tail_comparison(isotope, temperature, position, momentum):
    comparisons = {}
    for label, alternate in [('position', position), ('momentum', momentum)]:
        base = np.genfromtxt(OUT/f'{isotope}_{temperature}K_{label}.csv', delimiter=',', names=True)
        x = base[base.dtype.names[0]]
        for kind, column in [('density', base.dtype.names[1]), ('measure', base.dtype.names[2])]:
            y = base[column]; at = cutoff(x, y); threshold = float(y.max()*1e-6)
            delta = abs(alternate[kind]-y)
            at_delta = float(np.interp(at, x, delta))
            comparisons[f'{label}_{kind}'] = dict(cutoff_coordinate=at, threshold=threshold,
                empirical_absolute_delta_at_crossing=at_delta, delta_over_threshold=at_delta/threshold,
                integrated_tail_at_crossing=float(np.interp(at, x, base['tail_probability'])),
                integrated_tail_delta=float(np.interp(at, x, alternate['tail']-base['tail_probability'])),
                resolution='crossing empirically resolved for this comparison' if at_delta < threshold else 'crossing below comparison resolution')
    return comparisons


def controls():
    rows = json.loads((OUT/'single_site.json').read_text())
    controls = []
    for isotope in MASSES:
        coarse = solve(isotope, n=16, order=40)
        for temperature in (77, 300):
            base = next(row for row in rows if row['isotope'] == isotope and row['temperature_K'] == temperature)
            summary, pos, mom, _ = observables(coarse, temperature)
            deltas = {key: summary[key]-base[key] for key in ('mean_x2_A2', 'mean_kinetic_eV', 'position_density_1e6_A', 'momentum_density_1e6_A_inverse', 'probability_r_gt_080')}
            controls.append(dict(isotope=isotope, temperature_K=temperature, comparison='basis16 minus basis20, same quadrature40', deltas=deltas,
                                 tail_comparison=tail_comparison(isotope, temperature, pos, mom),
                                 low_cluster_excitation_delta_eV=(np.array([v['excitation_eV'] for v in levels(coarse)[:6]])-np.array([v['excitation_eV'] for v in base['levels'][:6]])).tolist()))
    # H 300K samples the widest positional support and is the deciding closure/tail control.
    base = next(row for row in rows if row['isotope'] == 'H' and row['temperature_K'] == 300)
    for label, spec, fine in [('quadrature48', solve('H', order=48), False),
                              ('fine_density_grid_and_angles', load_spectrum(OUT/'H_spectrum.npz'), True),
                              ('unmeasured_angular_and_exterior', solve('H', variant='closure'), False)]:
        summary, pos, mom, projection = observables(spec, 300, fine=fine)
        controls.append(dict(isotope='H', temperature_K=300, comparison=label, observables=summary,
                             deltas={key: summary[key]-base[key] for key in ('mean_x2_A2', 'mean_kinetic_eV', 'position_density_1e6_A', 'momentum_density_1e6_A_inverse', 'probability_r_gt_080')},
                             tail_comparison=tail_comparison('H', 300, pos, mom), low_clusters=levels(spec)[:6]))
        if label == 'fine_density_grid_and_angles':
            for name, data in [('position', pos), ('momentum', mom)]:
                table(OUT/f'H_300K_{name}_resolution.csv', ['coordinate', 'density', 'radial_measure', 'tail_probability'], list(data.values()))
        if label == 'unmeasured_angular_and_exterior':
            for name, data in [('position', pos), ('momentum', mom)]:
                table(OUT/f'H_300K_{name}_closure.csv', ['coordinate', 'density', 'radial_measure', 'tail_probability'], list(data.values()))
    dump(OUT/'numerical_controls.json', controls)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['produce', 'controls'])
    args = parser.parse_args()
    globals()[args.stage]()


if __name__ == '__main__':
    main()
