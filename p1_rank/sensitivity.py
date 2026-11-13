"""p1_rank: OAT long-time heat of the existing pass-zero chain, no DFT.

Run from the repository: python3 p1_rank/sensitivity.py > p1_rank/ranking.json
Source model p0_chain/pass0.py, dbab821a353faba367aa9c936ce51ba38655f8f2.
Stationary reduction is approved by axiom; it changes no equilibrium equation.
Zeros are retained, never replaced with a logarithmic floor. Positive decimal
spans describe this declared sample, not a continuous/global sensitivity index.
"""
from dataclasses import dataclass, replace, asdict
from pathlib import Path
import json
import sys
import time

import numpy as np
from scipy.constants import e, k
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from p0_chain.pass0 import Model, Chain, preparation, equilibrium


@dataclass(frozen=True)
class ChargeModel(Model):
    z_eff: float = 1.0

    @property
    def rho(self):
        return self.z_eff * super().rho


def prepared(eta):
    base = preparation("Pd/D", 0.8, 300.0, 0.037, eta)
    return ChargeModel(**asdict(base))


def closure(m, knobs, Cv, log_C):
    trial_C = np.exp(log_C)
    x = 1 / (1 + knobs.refill_s * knobs.per_site_s * trial_C)
    rate = m.density_m3 * x * knobs.per_site_s * trial_C
    heat = knobs.heat_fraction * knobs.Q_eV * e * rate
    T = m.T_K + knobs.thermal_s * heat / Cv
    actual_C = equilibrium(replace(m, density_m3=m.density_m3*x, T_K=T))[0][2]**2
    return float(actual_C-trial_C), float(heat), float(T), float(x)


def stationary(m, knobs, Cv):
    initial_C = float(equilibrium(m)[0][2]**2)
    if initial_C == 0 or knobs.per_site_s == 0:
        return dict(heat_W_m3=0.0, T_K=m.T_K, availability=1.0)
    # C decreases as trial C depletes n and increases T: one passive branch.
    log_C = brentq(lambda v: closure(m, knobs, Cv, v)[0], -200.0, np.log(initial_C))
    _, heat, T, x = closure(m, knobs, Cv, log_C)
    return dict(heat_W_m3=heat, T_K=T, availability=x)


def grids():
    geom = lambda a,b,n=9: np.geomspace(a,b,n).tolist()
    lin = lambda a,b,n=9: np.linspace(a,b,n).tolist()
    return {
        "eta": (lin(0,1,17), "c", "README diagnostic interval"),
        "k_site": (geom(1e-20,1,11), "c", "README unmeasured rate interval"),
        "f": (geom(1e-3,1)+[0.9,0.95,0.99,0.999], "a*", "README bound; positive grid, f=0 excluded by Model"),
        "z_eff": (geom(1e-3,1)+[0.9,0.95,0.99,0.999], "c", "axiom (0,1]; positive grid includes a normal-state regime"),
        "Ec_eV": (geom(0.01,1,13), "c", "README free-mode guess"),
        "gap_eV": (lin(0.037,0.065), "a/b", "sourced PdD peak-to-shoulder band, not an error bar"),
        "loading": (lin(0.6,1), "b", "p0_scales fixed-cell loading scenarios"),
        "T_K": (lin(300,600), "b", "p0_scales bath-temperature scenarios"),
        "species": (["D","H"], "b", "isolated nuclear-mass switch at fixed n and bare gap"),
        "Q_MeV": (geom(1,25), "c", "README no reaction selected"),
        "heat_fraction": (lin(0,1), "c", "README definition interval"),
        "escape_fraction": (lin(0,1,3), "c", "README definition interval"),
        "Cv": (geom(1e5,1e7), "a/b", "README loaded-hydride heat-capacity guess"),
        "refill_s": (geom(1,1e6), "c", "README unspecified transport/geometry"),
        "thermal_s": (geom(1,1e6), "c", "README unspecified thermal geometry"),
        "kappa": (geom(1e-3,1), "c", "README phenomenological field-loss guess"),
        "gamma": (geom(1e-3,1), "c", "README phenomenological matter-loss guess"),
    }


def varied(name, value, m, knobs, Cv):
    if name == "loading":
        return replace(m, density_m3=4*value/(4.084e-10)**3), knobs, Cv
    if name == "k_site":
        return m, replace(knobs, per_site_s=value), Cv
    if name == "Q_MeV":
        return m, replace(knobs, Q_eV=value*1e6), Cv
    if name == "Cv":
        return m, knobs, value
    if name in ("heat_fraction","escape_fraction","refill_s","thermal_s"):
        return m, replace(knobs, **{name:value}), Cv
    return replace(m, **{name:value}), knobs, Cv


def summarize(samples):
    powers = [x["heat_W_m3"] for x in samples]
    positive = [x for x in powers if x > 0]
    return dict(positive_decades=(float(np.log10(max(positive)/min(positive)))
                                  if positive else None),
                zero_nonzero_switch=bool(positive and min(powers) == 0),
                heat_min_W_m3=min(powers), heat_max_W_m3=max(powers),
                positive_min_W_m3=min(positive) if positive else None,
                sampled_values=samples)


def rank_case(eta, knobs, Cv, near_eta):
    m = prepared(eta)
    base = stationary(m, knobs, Cv)
    rows = []
    for name, (values, kind, basis) in grids().items():
        if name == "eta":
            values = values + [0.0, near_eta]
        samples = []
        for value in sorted(set(values)):
            result = stationary(*varied(name,value,m,knobs,Cv))
            samples.append(dict(value=value, **result))
        rows.append(dict(parameter=name,kind=kind,basis=basis,**summarize(samples)))
    rows.sort(key=lambda r: (-(r["positive_decades"] or 0),r["parameter"]))
    return dict(eta=eta,initial_model=asdict(m),central=base,ranking=rows,
                restricted=[r for r in rows if r["kind"].startswith(("a","b"))])


def main():
    start = time.perf_counter()
    m = prepared(0)
    critical_eta = m.f*np.tanh(m.gap*e/(2*k*m.T_K))-1/m.rho**2
    near_eta = 0.99*critical_eta
    knobs = Chain(per_site_s=9.1e-11,Q_eV=25e6)
    Cv = 3e6
    if "--pilot" in sys.argv:
        output = dict(critical_eta=critical_eta,near_eta=near_eta,
                      central_off=stationary(m,knobs,Cv),
                      central_near=stationary(prepared(near_eta),knobs,Cv),
                      high_rate=stationary(m,replace(knobs,per_site_s=1.0),Cv))
    else:
        output = dict(critical_eta=float(critical_eta),near_eta=float(near_eta),
                      chain_knobs=asdict(knobs),Cv=Cv,
                      cases=[rank_case(eta,knobs,Cv,near_eta) for eta in (0.0,near_eta)])
    output["elapsed_s"] = time.perf_counter()-start
    print(json.dumps(output,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
