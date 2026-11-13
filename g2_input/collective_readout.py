"""Fit the measured matched-sampling collective grids; retain raw energies/forces."""
import hashlib
import json
from pathlib import Path

import numpy as np
from mesh_readout import load_case


def read_grid(metal, mesh, state, direction):
    path = Path("collective_"+metal)/f"k{mesh}_{state}_{direction}"/"result.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    if (data["metal"], data["kmesh"], data["strain_case"], data["direction"], data["cutoff_eV"], data["smearing_eV"]) != (metal, mesh, state, direction, 500, 0.1):
        raise ValueError("Wrong collective grid profile")
    expected = [0., .06, .12, .18, .24, .30] if direction == "100" else [.06, .12, .18, .24, .30]
    if [p["u_A"] for p in data["points"]] != expected:
        raise ValueError("Incomplete collective displacement grid")
    return data, dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def fit_state(metal, mesh, state):
    base, base_source = load_case(metal, mesh, case="hydro_plus" if state == "hydro_plus" else None)
    K = base["K_Gamma_eV_A2"]
    grids, sources = {}, [base_source]
    for direction in ("100", "111"):
        grid, source = read_grid(metal, mesh, state, direction)
        if grid["deformation"] != base["deformation"] or grid["K2_force_eV_A2"] != K:
            raise ValueError("Collective grid and force curvature do not share one profile")
        grids[direction] = grid
        sources.append(source)
    E0 = grids["100"]["points"][0]["energy_free_eV"]
    fits = {}
    for direction, grid in grids.items():
        points = [p for p in grid["points"] if p["u_A"] > 0]
        u = np.array([p["u_A"] for p in points])
        energy = np.array([p["energy_free_eV"]-E0 for p in points])
        force = np.array([p["force_parallel_eV_A"] for p in points])
        feature = u**4/4
        K4 = float(feature @ (energy-K*u*u/2)/(feature @ feature))
        energy_error = energy-(K*u*u/2+K4*u**4/4)
        force_error = force+K*u+K4*u**3
        fits[direction] = dict(K2_force_eV_A2=K, K4_eV_A4=K4,
                               energy_RMS_meV=float(1000*np.sqrt(np.mean(energy_error**2))),
                               energy_max_error_meV=float(1000*np.max(np.abs(energy_error))),
                               force_prediction_RMS_eV_A=float(np.sqrt(np.mean(force_error**2))),
                               finite_step_quartic_curvature_estimate_eV_A2=K4*0.01**2,
                               points=[dict(u_A=p["u_A"], energy_free_eV=p["energy_free_eV"],
                                            energy_relative_eV=p["energy_free_eV"]-E0,
                                            force_parallel_eV_A=p["force_parallel_eV_A"]) for p in grid["points"]])
    return dict(metal=metal, kmesh=mesh, strain_case=state, deformation=base["deformation"],
                K2_force_eV_A2=K, zero_energy_free_eV=E0, directions=fits), sources


def main():
    states, sources = [], []
    for metal, mesh, state in (("Pd",12,"zero"),("Pd",12,"hydro_plus"),("Ni",8,"zero")):
        result, ss = fit_state(metal, mesh, state)
        states.append(result)
        sources.extend(ss)
        print(metal,state,'K2',result['K2_force_eV_A2'])
        for direction, fit in result['directions'].items():
            print(direction,'K4',fit['K4_eV_A4'],'energy RMS meV',fit['energy_RMS_meV'],'held-out force RMS',fit['force_prediction_RMS_eV_A'])
    output = dict(definition="Uniform collective displacement energy per primitive formula unit, fixed metal; not an isolated-ion potential.",
                  fit="K2 fixed by +/-0.01A central force response; one-parameter least-squares quartic fit to energies through0.30A. Forces are an independent prediction check. Ideal reference-site inversion supplies V(-u)=V(u).",
                  normalization="For N uniformly displaced formula units, total energy and collective-coordinate mass both scale with N; any single-ion/Einstein quantization is a separate model closure.",
                  states=states, source_files=sources)
    Path('collective-wells.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
