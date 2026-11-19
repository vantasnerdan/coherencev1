"""Read a mapped target-conditioned BO interaction; solve only with explicit preparation."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from .pairlaw import C, T, solve

SCHEMA = 'coherencev1-target-conditioned-pair-v1'


def vector(values, name):
    a = np.asarray(values, dtype=float)
    if a.ndim != 1 or not np.isfinite(a).all():
        raise ValueError(f'{name} must be a finite vector')
    return a


def require_text(fields, names):
    for name in names:
        if not isinstance(fields.get(name), str) or not fields[name].strip():
            raise ValueError(f'Provide explicit {name}')


def read_interaction(data):
    """On sampled support only; positive lowering is bare minus total."""
    if data.get('schema') != SCHEMA or data.get('observable') != 'target_conditioned_BO_total_interaction':
        raise ValueError('Require mapped total BO pair interaction, not U_phi, U_e, or drag')
    reference = data['reference']
    require_text(reference, ('source', 'configuration', 'energy_zero', 'periodic_to_isolated'))
    if reference.get('bare') != 'isolated_unit_charge_coulomb':
        raise ValueError('Bare reference must be the mapped isolated unit-charge Coulomb interaction')
    r = vector(data['R_A'], 'R_A')
    total = vector(data['V_pair_total_eV'], 'V_pair_total_eV')
    if len(r) < 2 or r[0] <= 0 or np.any(np.diff(r) <= 0) or total.shape != r.shape:
        raise ValueError('Require paired, strictly increasing positive radii and total energies')
    bare = C/r
    return r, total, dict(reference=reference, R_A=r.tolist(),
                          V_pair_total_eV=total.tolist(), V_bare_eV=bare.tolist(),
                          electronic_correction_eV=(total-bare).tolist(),
                          barrier_lowering_eV=(bare-total).tolist())


def prepared_ground_state(data, r, total, isotope, n):
    state = data.get('prepared_state')
    if state is None:
        raise ValueError('A collision curve alone does not specify a prepared pair; supply prepared_state')
    require_text(state, ('preparation', 'inner_continuation_source', 'outer_boundary'))
    confinement = vector(state['V_confinement_eV'], 'V_confinement_eV')
    if len(r) < 5 or confinement.shape != r.shape:
        raise ValueError('Ground state needs at least five radii with explicit confinement at each radius')
    origin = float(state['regular_potential_at_zero_eV'])
    if not np.isfinite(origin):
        raise ValueError('Supply a finite Coulomb-subtracted origin energy for the combined potential')
    if r[0] <= 1e-9:
        raise ValueError('First sample must be beyond the radial solver cusp radius (1e-9 Å)')
    rmax = float(r[-1])
    regular = CubicSpline(np.r_[0., r], np.r_[origin, total-C/r+confinement],
                          bc_type=((1, 0.), 'not-a-knot'), extrapolate=False)  # type: ignore[reportArgumentType]

    def potential(coordinate):
        # exp(log(rmax)) in the radial solver can round a few ulps past rmax.
        return C/coordinate + regular(np.minimum(coordinate, rmax))

    # The reused solver assumes a resolved single well and a matching region
    # on its inward forbidden side. Do not silently solve a clipped/multiwell input.
    probe = np.unique(np.r_[np.geomspace(1e-9, rmax, n), r, (r[:-1]+r[1:])/2])
    values = potential(probe)
    minima = np.flatnonzero((values[1:-1] < values[:-2]) & (values[1:-1] < values[2:]))
    if minima.size != 1 or not r[0] < probe[minima[0]+1] < rmax:
        raise ValueError('Expected one resolved well inside the sampled material curve')
    result = solve(potential, isotope, rmax, n=n)
    if (not np.isfinite(list(result.values())).all() or result['gap_eV'] <= 0
            or values[-1] <= result['E_eV']+4*result['gap_eV']
            or result['outer_tail_probability'] > 1e-8):
        raise ArithmeticError('Outer boundary does not resolve the confined ground state')
    return dict(isotope=isotope, state=state, solver=result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('curve', type=Path)
    parser.add_argument('--isotope', choices=tuple(T), help='Solve the explicitly prepared radial ground state')
    parser.add_argument('--grid', type=int, default=6500)
    args = parser.parse_args()
    if args.grid < 1000:
        parser.error('--grid must be at least 1000')
    raw = args.curve.read_bytes()
    try:
        data = json.loads(raw)
        r, total, out = read_interaction(data)
        out['source_sha256'] = hashlib.sha256(raw).hexdigest()
        if args.isotope is not None:
            out['prepared_ground_state'] = prepared_ground_state(data, r, total, args.isotope, args.grid)
    except (KeyError, TypeError, ValueError, ArithmeticError) as err:
        parser.exit(2, f'{parser.prog}: {err}\n')
    print(json.dumps(out, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
