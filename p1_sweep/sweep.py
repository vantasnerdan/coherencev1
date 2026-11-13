#!/usr/bin/env python3
"""Seeded pass-1 sweep of the SAME hypothetical pass0 equations, not LENR evidence.

Run beside the supplied pass0.py. Dependencies: Python stdlib, NumPy, SciPy.
OMP/BLAS threads=1. Outputs: sweep.csv and timeseries.npz. Seed=20261113.
Default 10000 rows follows the pilot CPU estimate; no equations are dropped.

rho_bare in [.1,10] is sampled as a DIMENSIONLESS input, z_eff uniformly in
(0,1], rho=z_eff*rho_bare. Both g and d use that SAME effective rho. This does
not infer screening microphysics. n_ref=1e28/m^3, Ec=.1eV are fixed dimensional
NORMALIZERS for power/temperature, not a claim that independently sampled rho
is the bare plasma frequency of that fixed material. During depletion,
rho(n)=rho_initial*sqrt(n/n_ref); r, f, eta and the bare gap stay fixed.

The underlying chain equations are unaltered; the supplied pass0 dependency
uses exact equilibrium-polarization coordinates with implicit Radau for the
sampled stiff rates, and the EXACT invariant zero-output normal branch.
'Long time' means t=20*max(tau_refill,tau_th), not an unobserved infinite limit.
Behaviour labels are finite-trajectory descriptions defined in behaviour().
No undefined locked frequency, domain size, or domain lifetime is output.

200 designed trajectories use moderate losses and small +/- distance from
threshold: a stated near-threshold slaved-regime subset for the SINDy consumer,
not a claim of a global one-variable closure. All five states and exact RHS
values are saved; no first-order approximation is used to generate them.
"""
import argparse
import csv
from dataclasses import dataclass
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
from pass0 import Model, Chain, chain, equilibrium, normal_matrix, rhs, e, k

SEED = 20261113
NREF = 1e28
EC_EV = 0.1
GROUP_NAMES = ["rho_bare", "z_eff", "rho", "r", "theta", "eta", "f", "kappa", "gamma"]


@dataclass(frozen=True)
class GroupModel(Model):
    """Only a group-to-unit adapter; inherited Hamiltonian and ODEs unchanged."""
    rho_initial: float = 1.0
    reference_density: float = NREF

    @property
    def rho(self):
        return self.rho_initial*np.sqrt(self.density_m3/self.reference_density)


def model(row):
    return GroupModel(NREF, Ec_eV=EC_EV, gap_eV=row["r"]*EC_EV,
                      T_K=row["theta"]*row["r"]*EC_EV*e/k,
                      kappa=row["kappa"], gamma=row["gamma"],
                      f=row["f"], eta=row["eta"], rho_initial=row["rho"])


def loguniform(rng, lo, hi):
    return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))


def draw(rng):
    rho_bare = loguniform(rng, .1, 10)
    z = float(1-rng.random())
    return dict(rho_bare=rho_bare, z_eff=z, rho=z*rho_bare,
                r=loguniform(rng, .3, 3), theta=loguniform(rng, .05, 2),
                eta=float(rng.random()), f=float(rng.uniform(.1, 1)),
                kappa=loguniform(rng, 1e-3, 1), gamma=loguniform(rng, 1e-3, 1))


def chain_inputs(row, rng):
    row.update(k_site_s=loguniform(rng, 1e-20, 1), Q_MeV=float(rng.uniform(1, 25)),
               heat_fraction=float(rng.random()), Cv_J_m3_K=loguniform(rng, 1e5, 1e7),
               tau_th_s=loguniform(rng, 1, 1e6), tau_refill_s=loguniform(rng, 1, 1e6))
    return row


def trajectory_design(rng, count):
    rows = []
    while len(rows) < count:
        i = len(rows)
        # Moderate-loss center-manifold regime, both sides; all are full ODEs.
        r, theta, f = loguniform(rng, .5, 1.5), loguniform(rng, .1, .6), float(rng.uniform(.4, 1))
        w = np.tanh(1/(2*theta))
        eta = float(rng.uniform(0, .5*f*w))
        mu = loguniform(rng, .002, .015) * (1 if i % 2 else -1)
        rho = np.sqrt((1+mu)/(f*w-eta))
        z = float(rng.uniform(.5, 1))
        row = dict(rho_bare=float(rho/z), z_eff=z, rho=float(rho), r=r,
                   theta=theta, eta=eta, f=f, kappa=loguniform(rng, .2, .8),
                   gamma=loguniform(rng, .1, .8), design="trajectory")
        poles = np.linalg.eigvals(normal_matrix(model(row)))
        j = int(np.argmax(poles.real))
        # Explicit reduction scope: real slow mode, >=4x faster stable poles.
        # This selects the time-series regime, not the broad sweep rows.
        if abs(poles[j].imag) > 1e-10 or min(-np.delete(poles, j).real) < 4*abs(poles[j].real):
            continue
        rows.append(chain_inputs(row, rng))
    return rows


def design(rows, trajectories, seed):
    rng = np.random.default_rng(seed)
    result = trajectory_design(rng, trajectories)
    remaining = rows-trajectories
    near = 2*(remaining//10)
    for _ in range(remaining-near):
        result.append(chain_inputs(draw(rng) | {"design": "bulk"}, rng))
    for i in range(near):
        while True:
            row = draw(rng)
            eta_star = row["f"]*np.tanh(1/(2*row["theta"]))-1/row["rho"]**2
            epsilon = loguniform(rng, 1e-4, .03)
            eta = eta_star + (epsilon if i % 2 else -epsilon)
            if 0 <= eta <= 1:
                break
        row.update(eta=float(eta), design="near_normal" if i % 2 else "near_ordered")
        result.append(chain_inputs(row, rng))
    return result


def behaviour(result):
    series = result["rows"]
    power = np.array([x["heat_W_m3"] for x in series])
    C = np.array([x["coherent_proxy"] for x in series])
    if np.all(C == 0):
        return "switched_off"
    tail = power[-8:]
    # Relative to observed power, with a stated 1e-12 W/m3 absolute floor.
    scale = max(float(np.max(abs(tail))), 1e-12)
    if np.ptp(tail) <= 1e-3*scale:
        return "steady"
    turns = np.diff(np.sign(np.diff(power[-24:])))
    if np.count_nonzero(turns) >= 4 and np.ptp(tail) > .01*scale:
        return "oscillating"
    if tail[-1] < tail[0]:
        return "decaying"
    # A rising transient is not one of the four authorized behavior labels.
    # Keep it explicit rather than silently calling it steady.
    return "rising_unsettled"


def evaluate(row, run_id):
    m = model(row)
    state, du, df, gain = equilibrium(m)
    w = np.tanh(1/(2*row["theta"]))
    denominator = row["f"]*w-row["eta"]
    critical = float(1/np.sqrt(denominator)) if denominator > 0 else None
    growth = float(np.linalg.eigvals(normal_matrix(m)).real.max())
    end = 20*max(row["tau_refill_s"], row["tau_th_s"])
    times = np.r_[0, np.geomspace(end*1e-8, end*.1, 32), np.linspace(end*.1, end, 33)[1:]]
    knobs = Chain(per_site_s=row["k_site_s"], Q_eV=row["Q_MeV"]*1e6,
                  heat_fraction=row["heat_fraction"], escape_fraction=.1,
                  refill_s=row["tau_refill_s"], thermal_s=row["tau_th_s"])
    output = chain(m, knobs, row["Cv_J_m3_K"], times)
    last = output["rows"][-1]
    return {"run_id": run_id, **row, "n_ref_m3": NREF, "Ec_eV": EC_EV,
            "T0_K": m.T_K, "threshold_margin": float(denominator*row["rho"]**2-1),
            "ordered": int(gain > 1), "rho_c": critical, "q": float(state[0]),
            "C": float(state[2]**2), "dU_Ec": du, "dF_Ec": df,
            "growth_per_tau": growth, "heat_long_W_m3": last["heat_W_m3"],
            "T_long_K": last["T_K"], "n_long_over_n0": last["availability_n_over_n0"],
            "C_long": last["coherent_proxy"], "long_time_s": end,
            "behaviour": behaviour(output), "energy_partition_relative": output["relative_energy_partition_residual"],
            "chain_rhs_evaluations": output["rhs_evaluations"]}


def full_trajectory(row, points):
    m = model(row)
    eq, _, _, gain = equilibrium(m)
    values, vectors = np.linalg.eig(normal_matrix(m))
    j = int(np.argmax(values.real))
    lam = values[j]
    if abs(lam.imag) > 1e-10:
        raise ValueError("Designed SINDy trajectory is not a real soft-mode branch")
    vector = vectors[:, j].real
    vector /= vector[0]
    seed = 1e-5
    initial = np.r_[seed*vector, 0.0]
    w = np.tanh(1/(2*row["theta"]))
    initial[4] = -np.sqrt(w*w-initial[2]**2-initial[3]**2)
    efolds = np.log(max(abs(eq[0])/seed, 1))+12 if gain > 1 else 16.0
    duration = max(16/min(m.kappa, m.gamma), efolds/abs(lam.real))
    times = np.linspace(0, duration, points)
    sol = solve_ivp(rhs, (0, duration), initial, args=(m,), method="BDF",
                    t_eval=times, rtol=1e-6, atol=1e-10)
    if not sol.success:
        raise RuntimeError(sol.message)
    derivative = np.array([rhs(t, y, m) for t, y in zip(times, sol.y.T)])
    return times, sol.y.T, derivative, float(eq[0]), float(lam.real), sol.nfev


COLUMN_DICTIONARY = {
    "run_id/design": "seeded row identity and bulk/near-threshold/trajectory sampling stratum",
    "rho_bare/z_eff/rho": "bare sampled coupling, effective-charge factor, rho=z_eff*rho_bare; d and g both use rho",
    "r/theta/eta/f/kappa/gamma": "Egap/Ec, kBT/Egap, diamagnetic fraction, oscillator strength, two losses/omega_c",
    "k_site_s/Q_MeV/heat_fraction": "per-active-site events/s, free MeV/event, event-energy fraction to heat",
    "Cv_J_m3_K/tau_th_s/tau_refill_s": "effective heat capacity, cooling time and availability-refill time",
    "n_ref_m3/Ec_eV/T0_K": "fixed power normalizers and group-derived initial temperature",
    "threshold_margin": "(f*tanh(1/(2theta))-eta)*rho^2-1; positive is ordered",
    "ordered/rho_c": "initial equilibrium order flag; effective critical rho, blank if no finite threshold",
    "q/C/dU_Ec/dF_Ec": "initial equilibrium field, squared transverse spin, internal/free coherent-minus-normal energies per site/Ec",
    "growth_per_tau": "largest real normal field-spin eigenvalue, dimensionless; longitudinal decoupled pole excluded",
    "heat_long_W_m3/T_long_K/n_long_over_n0/C_long/long_time_s": "actual finite-chain endpoint at 20*max(cooling,refill), not infinite-time certification",
    "behaviour": "finite trace: switched_off if C identically0; steady if last8 P range<=.001 max(P,1e-12 W/m3); oscillating if >=4 tail slope turns and >1% variation; otherwise decaying or explicit rising_unsettled",
    "energy_partition_relative/chain_rhs_evaluations": "Q partition residual and solver work; 0 evaluations denotes exact invariant normal branch"
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=10000)
    parser.add_argument("--trajectories", type=int, default=200)
    parser.add_argument("--points", type=int, default=513)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    if not 2 <= args.trajectories <= args.rows:
        raise ValueError("Require at least two trajectories and rows>=trajectories")
    rows = design(args.rows, args.trajectories, args.seed)
    start = time.perf_counter()
    table = [evaluate(row, i) for i, row in enumerate(rows)]
    sweep_s = time.perf_counter()-start
    with open("sweep.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    start = time.perf_counter()
    trajectories = [full_trajectory(row, args.points) for row in rows[:args.trajectories]]
    trajectory_s = time.perf_counter()-start
    np.savez_compressed("timeseries.npz", time=np.array([x[0] for x in trajectories]),
                        state=np.array([x[1] for x in trajectories]),
                        rhs=np.array([x[2] for x in trajectories]),
                        equilibrium_q=np.array([x[3] for x in trajectories]),
                        growth=np.array([x[4] for x in trajectories]),
                        groups=np.array([[row[n] for n in GROUP_NAMES] for row in rows[:args.trajectories]]),
                        group_names=np.array(GROUP_NAMES), state_names=np.array(["q", "p", "sx", "sy", "sz"]),
                        row_id=np.arange(args.trajectories))
    print(json.dumps({"rows": args.rows, "trajectories": args.trajectories, "points": args.points,
                      "seed": args.seed, "sweep_s": sweep_s, "trajectory_s": trajectory_s,
                      "ordered_rows": sum(row["ordered"] for row in table),
                      "behaviour_counts": {label: sum(row["behaviour"] == label for row in table)
                                            for label in sorted({row["behaviour"] for row in table})},
                      "max_energy_partition_relative": max(row["energy_partition_relative"] for row in table),
                      "max_trajectory_endpoint_q_error": float(max(abs(x[1][-1, 0]-x[3]) for x in trajectories)),
                      "trajectory_rhs_evaluations": sum(x[5] for x in trajectories),
                      "columns": COLUMN_DICTIONARY,
                      "five_sample_rows": [table[i] for i in np.linspace(0, args.rows-1, 5, dtype=int)]},
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

