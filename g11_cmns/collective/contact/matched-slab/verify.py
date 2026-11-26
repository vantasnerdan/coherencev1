"""Independent finite native-path and long-join numerical cross-check, not an enclosure.

The explicit Hermite basis, separately assembled Hamiltonian and u=R exp(z)
Coulomb propagation do not call the producer or G11 solve/contact_log.
"""
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh_tridiagonal
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
NATIVE = ROOT / "g10_quantum_site/pdd-covered-slab-pair/native"


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(ok, message):
    if not ok:
        raise ValueError(message)


def hermite(q, x, y, slopes, derivative=0):
    """Explicit interval basis, independently of CubicHermiteSpline."""
    x, y, slopes, q = map(np.asarray, (x, y, slopes, q))
    j = np.clip(np.searchsorted(x, q, side="right")-1, 0, len(x)-2)
    width = x[j+1]-x[j]
    t = (q-x[j])/width
    if derivative:
        return ((6*t*t-6*t)*y[j]+(-6*t*t+6*t)*y[j+1])/width + (3*t*t-4*t+1)*slopes[j] + (3*t*t-2*t)*slopes[j+1]
    return (2*t**3-3*t*t+1)*y[j]+(-2*t**3+3*t*t)*y[j+1]+width*((t**3-2*t*t+t)*slopes[j]+(t**3-t*t)*slopes[j+1])


def ray(q, data, derivative=0):
    x = [0., *data["q_A"]]
    y = [0., *data["E_eV"]]
    s = [data["dE_dq_at_zero_eV_A"], *data["dE_dq_eV_A"]]
    check(data["origin_cut_A"] == 0., "Independent ray requires zero cut")
    end = x[-1]
    t = np.maximum(np.asarray(q)-end, 0.)
    outside = s[-1]+data["K_eV_A2"]*t if derivative else y[-1]+s[-1]*t+.5*data["K_eV_A2"]*t*t
    return np.where(np.asarray(q)>end, outside, hermite(q, x, y, s, derivative))


def outer(r, desc, derivative=0):
    u = (desc["R0_A"]-np.asarray(r))/2
    q = np.abs(u)
    w = desc["connected"]
    end = w["u_A"][-1]
    wc = hermite(np.clip(u, 0., end), [0., *w["u_A"]], [0., *w["W_eV"]], [0., *w["dW_du_eV_A"]], derivative)
    wc = np.where(u>end, w["dW_du_eV_A"][-1] if derivative else w["W_eV"][-1]+w["dW_du_eV_A"][-1]*(u-end), wc)
    wc = np.where(u<=0., 0., wc)
    total = np.zeros_like(u)
    for site in desc["sites"]:
        model = desc["wells"][site["well"]]
        a, b = ray(q, model["inward"], derivative), ray(q, model["outward"], derivative)
        total += np.where(u<0., -b if derivative else b, a)
    return -.5*(total+wc) if derivative else total+wc


def main():
    result = read(HERE / "result.json")
    for name, digest in result["source_sha256"].items():
        check(sha(ROOT/name) == digest, "Changed bound source " + name)
    spec = importlib.util.spec_from_file_location("g11_for_outer_probe", HERE.parent / "contact.py")
    g11 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g11)
    native_errors, consumer_errors, hermite_errors, zero_rows = [], [], [], []
    for mesh in (3, 2):
        desc = read(HERE / f"pair-m{mesh}.json")
        pair = g11.SitePair(HERE / f"pair-m{mesh}.json")
        p = next(p for p in pair.potentials() if p.include_flagged)
        zero = read(NATIVE / f"m{mesh}-zero-0.json")
        n = np.array(zero["point"]["axis_A_to_B"])
        base = np.array(zero["positions_A"])
        for u in (.06, .12, .24, .36, .48):
            a, b, ab = [read(NATIVE / f"m{mesh}-{role}-{round(100*u)}.json") for role in ("a", "b", "ab")]
            va, vb, total = [v["energy_free_eV"]-zero["energy_free_eV"] for v in (a,b,ab)]
            da = -np.dot(a["forces_eV_A"][26], n)
            db = np.dot(b["forces_eV_A"][28], n)
            dt = -np.dot(ab["forces_eV_A"][26], n)+np.dot(ab["forces_eV_A"][28], n)
            expected = base.copy(); expected[26] += u*n; expected[28] -= u*n
            check(np.max(np.abs(expected-ab["positions_A"])) < 3e-12, "AB path differs")
            r = pair.r0-2*u
            native_errors.extend([abs(float(outer(r, desc))-total), abs(-2*float(outer(r, desc, 1))-dt)])
            consumer_errors.extend([abs(float(p.lattice(r))-total), abs(-2*float(p.lattice(r,1))-dt)])
            for well, energy, slope in (("A",va,da),("B",vb,db)):
                native_errors.extend([abs(float(ray(u,desc["wells"][well]["inward"]))-energy),
                                      abs(float(ray(u,desc["wells"][well]["inward"],1))-slope)])
        for r in np.linspace(.9, 4.1, 151):
            hermite_errors.extend([abs(float(outer(r,desc))-float(p.lattice(r))), abs(float(outer(r,desc,1))-float(p.lattice(r,1)))])
        for well, idx, sign in (("A",26,1),("B",28,-1)):
            expected = -sign*float(np.dot(zero["forces_eV_A"][idx],n))
            inward = desc["wells"][well]["inward"]["dE_dq_at_zero_eV_A"]
            outward = desc["wells"][well]["outward"]["dE_dq_at_zero_eV_A"]
            check(inward == expected and outward == -expected, "Lost signed native zero force/C1")
            zero_rows.append(dict(mesh=mesh, site=well, inward_slope_eV_A=inward, outward_slope_eV_A=outward))
    check(max(native_errors+consumer_errors+hermite_errors) < 5e-12, "Finite native/consumer proof failed")

    # Validate both preserved outward samples against their own mesh2 zero, not inward mirrors.
    surface = read(ROOT/"g10_quantum_site/pdd-surface-vacancy/site-wells.json")["surface"]
    outward_errors = []
    for mesh in (2,3):
        desc = read(HERE/f"pair-m{mesh}.json")
        nz = desc["axis"][2]
        for well,label,sign in (("A","subsurface",1),("B","adsorbed",-1)):
            rows = {p["point"]["displacement_A"]:p for p in surface["directions"][label] if p["point"]["axis"] == "normal"}
            for z in (.12,.36):
                native = rows[-sign*z]
                q = z/nz
                energy = native["energy_free_eV"]-surface["zero"]["energy_free_eV"]
                slope = sign*nz*native["force_parallel_eV_A"]
                outward_errors.extend([abs(float(ray(q,desc["wells"][well]["outward"]))-energy),
                    abs(float(ray(q,desc["wells"][well]["outward"],1))-slope)])
    check(max(outward_errors) < 2e-12, "Outward native normal proxy not reproduced")
    refined_control = result["single_finer_mesh_control"]["mesh_4_triplet"]
    zero4,a4,b4,ab4 = [read(NATIVE/f"m4-{role}-{amp}.json") for role,amp in (("zero",0),("a",36),("b",36),("ab",36))]
    w4 = math.fsum((ab4["energy_free_eV"],-a4["energy_free_eV"],-b4["energy_free_eV"],zero4["energy_free_eV"]))
    check(abs(w4-refined_control["connected_eV"]) < 2e-12, "Finer control changed")

    # Separate numerical route for the named primary long join only.
    desc = read(HERE/"pair-m3.json")
    join = 1.6
    coeff, coulomb, kb = [result["constants"][key] for key in ("kinetic_eV_A2","C_eV_A","KB_eV_K")]
    cr = np.array(desc["inner"]["R_A"])
    el = PchipInterpolator(cr, np.array(desc["inner"]["V_total_eV"])-coulomb/cr)
    def inner(r):
        r = np.asarray(r)
        a = cr[0]
        smooth = float(el(a))+float(el(a,1))*(r*r-a*a)/(2*a)
        return coulomb/r+np.where(r<a, smooth, el(np.clip(r,a,cr[-1])))
    shift = float(outer(join,desc)-inner(join))
    def potential(r):
        r = np.asarray(r)
        return np.where(r>=join,outer(r,desc),inner(r)+shift)
    h, rmax, ecut = .00075, 9.5, 2.
    ngrid = round(rmax/h)-1; h = rmax/(ngrid+1)
    r = np.arange(1,ngrid+1)*h
    energy, vec = eigh_tridiagonal(2*coeff/h**2+potential(r), np.full(ngrid-1,-coeff/h**2),
        select="v", select_range=(-1.,ecut), tol=1e-11, lapack_driver="stebz")
    contacts = []
    for j,en in enumerate(energy):
        # u=R exp(z), v=dz/d(log R): v'=R^2(U-E)/t-v-v^2.
        amplitude = vec[:,j]/math.sqrt(h)
        match = np.flatnonzero((r<join+1.) & (potential(r)>en+.22))[-1]
        check(abs(amplitude[match])>1e-8*np.max(np.abs(amplitude)), "Unresolved independent match")
        start = 1e-9; a = coulomb*start/(2*coeff)
        def ode(x,zv):
            rr = math.exp(x); z,v = zv
            return [v,rr*rr*(float(potential(rr))-en)/coeff-v-v*v]
        sol = solve_ivp(ode,(math.log(start),math.log(r[match])),[math.log1p(a),a/(1+a)],method="DOP853",rtol=5e-10,atol=1e-10)
        check(sol.success,"Independent regular-origin propagation failed")
        log_slope = math.log(abs(amplitude[match]))-math.log(r[match])-sol.y[0,-1]
        contacts.append(2*log_slope-math.log(4*math.pi))
    named = next(v for v in result["variants"] if v["case"]=="m3" and v["model"]["join_A"]==join)
    numerical = []
    for prep,temp in (("cold_ground",None),("thermal_77K",77.),("thermal_300K",300.)):
        logc = contacts[0]/math.log(10) if temp is None else (logsumexp(-(energy-energy[0])/(kb*temp)+contacts)-logsumexp(-(energy-energy[0])/(kb*temp)))/math.log(10)
        reported = next(row["log10_contact_A_minus3"] for row in named["observables"] if row["preparation"]==prep)
        numerical.append(dict(preparation=prep, independently_computed_log10_rho_A_minus3=float(logc), reported=reported, change_dex=float(logc-reported)))
    check(max(abs(row["change_dex"]) for row in numerical) < .001,"Named long-join contact not reproduced independently")
    proof = dict(schema="covered-slab-independent-finite-proof-v1", result_sha256=sha(HERE/"result.json"), verifier_sha256=sha(Path(__file__)),
        descriptor_sha256={str((HERE/name).relative_to(ROOT)):sha(HERE/name) for name in ("pair-m3.json","pair-m2.json","pair-short-terminal.json")},
        explicit_native_path_max_error_eV_or_eV_A=max(native_errors), actual_G11_native_max_error_eV_or_eV_A=max(consumer_errors),
        explicit_basis_vs_G11_max_error_eV_or_eV_A=max(hermite_errors), outward_proxy_max_error_eV_or_eV_A=max(outward_errors), signed_zero_slopes=zero_rows,
        mesh4_W_eV=w4, independent_numerics=dict(join_A=join,h_A=h,rmax_A=rmax,ecut_eV=ecut,levels=len(energy),observables=numerical),
        scope="Explicit finite native/consumer agreement and separate finite numerical long-join check, not rigorous interval bounds, continuum proof, electronic error bound or physical contact measurement.")
    (HERE/"proof.json").write_text(json.dumps(proof,indent=2,allow_nan=False)+"\n")
    print(json.dumps(proof,indent=2))


if __name__ == "__main__":
    main()
