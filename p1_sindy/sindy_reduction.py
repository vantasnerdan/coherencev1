#!/usr/bin/env python3
"""Actual PySINDy fit of a near-threshold real soft-mode equation.

Reads timeseries.npz from sweep.py: the FULL five-variable model generated
these data. Fits ONLY the declared slaved, real slow-mode subset after fast
relaxation, not a global first-order representation of oscillatory dynamics.
Time is tau=omega_c*t; q is the original field quadrature, not rescaled.

Physics-informed parameter features are supplied, not claimed discovered:
 w=tanh(1/(2theta)), g=rho*sqrt(f*r)/2, K=1+eta*rho^2,
 A=4*g^2*w/r-K,
 D=2*kappa+(4*g^2*w/r)*gamma/(r^2+gamma^2),
 B=32*g^4/r^3*[w-a*(1-w^2)], a=1/(2theta).
The slow-response expansion suggests q'=alpha*q-beta*q^3 with alpha=A/D,
beta=B/D; the coefficients of six competing library terms are actually fit
by PySINDy, not assigned to +/-1. Odd state terms reflect the exact Z2 parity.
Parameterized library: {1,alpha,beta} tensor {q,q^3}.

A fixed STLSQ threshold .05 is used once. No threshold tuning/holdout scheme.
Comparison runs are in-sample consumer demonstrations, not independent review.
States before 8/min(kappa,gamma) and |q|<=1e-8 are excluded from a FIRST-order
fit; the original untrimmed five-state data remain unchanged and available.
This declares the reduction/absolute-resolution scope, not a claim that fast
variables or tiny decaying amplitudes obey the reduced equation.
"""
import csv
import json
import time
import numpy as np
import pysindy as ps
from scipy.integrate import solve_ivp


def parameter_features(group, names):
    value = dict(zip(names, group))
    rho, r, theta = value["rho"], value["r"], value["theta"]
    w = np.tanh(1/(2*theta))
    g = rho*np.sqrt(value["f"]*r)/2
    K = 1+value["eta"]*rho*rho
    A = 4*g*g*w/r-K
    D = 2*value["kappa"]+(4*g*g*w/r)*value["gamma"]/(r*r+value["gamma"]**2)
    B = 32*g**4/r**3*(w-(1-w*w)/(2*theta))
    return float(A/D), float(B/D), value


def main():
    start = time.perf_counter()
    data = np.load("timeseries.npz")
    names = data["group_names"].tolist()
    xs, derivatives, controls, times = [], [], [], []
    selections, drivers = [], []
    for t, state, derivative, group in zip(data["time"], data["state"], data["rhs"], data["groups"]):
        alpha, beta, value = parameter_features(group, names)
        mask = (t >= 8/min(value["kappa"], value["gamma"])) & (abs(state[:, 0]) > 1e-8)
        if np.count_nonzero(mask) < 4:
            raise ValueError("Declared reduction regime has insufficient resolved samples")
        selections.append(mask)
        drivers.append((alpha, beta))
        xs.append(state[mask, :1])
        derivatives.append(derivative[mask, :1])
        controls.append(np.tile([alpha, beta], (np.count_nonzero(mask), 1)))
        times.append(t[mask])
    state_library = ps.CustomLibrary([lambda q: q, lambda q: q**3],
                                     function_names=[lambda q: q, lambda q: f"{q}^3"])
    library = ps.ParameterizedLibrary(parameter_library=ps.PolynomialLibrary(degree=1),
                                      feature_library=state_library, num_features=1, num_parameters=2)
    estimator = ps.SINDy(feature_library=library,
                         optimizer=ps.STLSQ(threshold=.05, alpha=1e-12,
                                            normalize_columns=False, unbias=True))
    estimator.fit(xs, t=times, x_dot=derivatives, u=controls,
                  feature_names=["q", "alpha", "beta"])
    features = estimator.get_feature_names()
    coefficients = estimator.coefficients()[0]
    prediction = estimator.predict(xs, u=controls)
    derivative_error = np.concatenate([a-b for a, b in zip(prediction, derivatives)])
    comparison = []
    run_metrics = []
    # Two sides, two parts of the sampled regime, all explicitly in-sample.
    for run in (0, 1, 100, 101):
        t = data["time"][run]
        state = data["state"][run]
        first = int(np.flatnonzero(selections[run])[0])
        alpha, beta = drivers[run]
        u = np.array([[alpha, beta]])

        def reduced_rhs(tau, q):
            return estimator.predict(np.array([[q[0]]]), u=u)[0]

        sol = solve_ivp(reduced_rhs, (t[first], t[-1]), state[first, :1],
                        t_eval=t[first:], method="DOP853", rtol=1e-7, atol=1e-10)
        if not sol.success:
            raise RuntimeError(sol.message)
        full = state[first:, 0]
        error = sol.y[0]-full
        run_metrics.append({"run_id": run, "alpha": alpha, "beta": beta,
                            "growth_full_per_tau": float(data["growth"][run]),
                            "start_tau": float(t[first]), "end_tau": float(t[-1]),
                            "q_full_end": float(full[-1]), "q_reduced_end": float(sol.y[0, -1]),
                            "max_abs_q_error": float(np.max(abs(error))),
                            "rms_q_error": float(np.sqrt(np.mean(error**2)))})
        for j in np.linspace(0, len(full)-1, 5, dtype=int):
            comparison.append({"run_id": run, "tau": float(t[first+j]),
                               "q_full": float(full[j]), "q_reduced": float(sol.y[0, j]),
                               "abs_error": float(abs(error[j]))})
    result = {"pysindy_version": ps.__version__, "trajectories": len(xs),
              "resolved_fitting_samples": int(sum(len(x) for x in xs)),
              "features": features, "coefficients": coefficients.tolist(),
              "equation": estimator.equations(precision=12)[0],
              "parameter_functions": {"alpha": "[f*rho^2*w-(1+eta*rho^2)]/[2*kappa+f*rho^2*w*gamma/(r^2+gamma^2)]",
                                      "beta": "[2*f^2*rho^4/r*(w-(1-w^2)/(2*theta))]/[2*kappa+f*rho^2*w*gamma/(r^2+gamma^2)]",
                                      "w": "tanh(1/(2*theta))"},
              "derivative_training_mse": float(np.mean(derivative_error**2)),
              "run_metrics": run_metrics, "comparison_table": comparison,
              "elapsed_s": time.perf_counter()-start,
              "limits": "physics-informed parameter features; near-threshold real slow branch after fast slaving only; not global oscillatory dynamics or physical LENR evidence"}
    with open("sindy-results.json", "w") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    with open("sindy-comparison.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
