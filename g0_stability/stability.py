"""Exact normal modes of the bright block, no rotating-wave approximation.

Source: Gamberale and Modanese, arXiv:2205.14015v3, Eqs. 2, 7, 14, 19, 20, 58-69.
Compares the reduction as printed with a completion that uses the same
vector potential in both the p.A and the A^2 terms.
Units: omega' = 1. Each bright block repeats in 3 Cartesian directions;
the 3(N-1) dark matter modes keep frequency^2 = 1. No fitted parameters.
"""
import json
import numpy as np

S = {"mode_sum": 8*np.pi/3, "fixed_axis_authors": 4.0}


def modes(eps, s):
    # Printed: H = .5(x*x+p*p+q*q+r*r) + sqrt(s)*eps*p*q, so nu^2 = 1 +/- sqrt(s)*eps.
    printed = 1 - np.sqrt(s)*eps
    # Same A: H = .5(P+gQ)^2 + .5X^2 + .5Pi^2 + .5w^2 Q^2,
    # g^2 = s*eps^2, w^2 = 1-eps^2. Roots of z^2 - Bz + C = 0.
    B = 2 + (s-1)*eps**2
    C = (1-eps)*(1+eps)
    disc = (s*eps**2)**2 + 2*s*eps**2*(2-eps**2) + eps**4
    consistent = 2*C/(B+np.sqrt(disc))
    return float(printed), float(consistent)


if __name__ == "__main__":
    grid = sorted(set([0, .1, .25, float(1/np.sqrt(S['mode_sum'])),
                       .4, .5, .52, .56, .69, .75, .9, .99, 1]))
    rows = []
    for e in grid:
        row = {"eps": e}
        for label, s in S.items():
            p, c = modes(e, s)
            row[label+"_printed_nu2_over_W2"] = p
            row[label+"_same_A_nu2_over_W2"] = c
        rows.append(row)
    print(json.dumps({"units": "nu^2 / omega_prime^2", "rows": rows,
        "printed_critical_eps": {k: float(1/np.sqrt(s)) for k, s in S.items()},
        "same_A": "strictly positive for 0<=eps<1; zero at bare-omega=0 endpoint",
        "physical_limit": "dipole-approximation quadratic model only, not a no-go theorem for a real hydride"},
        indent=2, allow_nan=False))
