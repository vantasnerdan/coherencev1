"""Four driven lossy spins: QuTiP, local quimb tensor gates, and dense reference."""
from functools import reduce
import json

import numpy as np
import quimb
import quimb.tensor as qtn
from scipy.linalg import expm
import qutip as qt


def kron(*matrices):
    return reduce(np.kron, matrices)


def liouvillian(h, *jumps):
    """Row-major rho[a,b]; kets and bras stay separately indexed."""
    eye = np.eye(len(h))
    generator = -1j * (kron(h, eye) - kron(eye, h.T))
    for jump in jumps:
        n = jump.conj().T @ jump
        generator += kron(jump, jump.conj()) - 0.5 * (kron(n, eye) + kron(eye, n.T))
    return generator


def main():
    sites, dt, steps = 4, 0.025, 80
    eye = np.eye(2)
    down = np.array([[0., 1.], [0., 0.]])
    x = down + down.T
    h = 0.35 * x + np.diag([0., 0.20])
    jump = np.sqrt(0.27) * down
    pair_h = 0.7 * (kron(down, down.T) + kron(down.T, down))

    def place(op, location, width=1):
        return kron(*([eye] * location), op, *([eye] * (sites-location-width)))

    full_h = sum(place(h, i) for i in range(sites))
    full_h += sum(place(pair_h, i, 2) for i in range(sites-1))
    jumps = [place(jump, i) for i in range(sites)]
    psi = np.eye(2**sites)[8]
    rho0 = np.outer(psi, psi)
    end = dt * steps
    reference = (expm(end * liouvillian(full_h, *jumps)) @ rho0.ravel()).reshape(rho0.shape)
    solution = qt.mesolve(qt.Qobj(full_h), qt.Qobj(rho0), [0, end],
                          [qt.Qobj(j) for j in jumps],
                          options={"rtol": 1e-11, "atol": 1e-13})
    qutip_rho = solution.states[-1].full()

    half_site = expm(dt / 2 * liouvillian(h, jump))

    def pair_gate(t):
        u = expm(-1j * t * pair_h)
        return kron(u, u.conj()).reshape([2] * 8).transpose(0,2,1,3,4,6,5,7).reshape(16,16)

    half_even, odd = pair_gate(dt / 2), pair_gate(dt)
    state = qtn.MPS_product_state(np.eye(4, dtype=complex)[[3,0,0,0]])
    for _ in range(steps):
        for site in range(sites):
            state.gate_(half_site, site, contract=True)
        for i, gate in ((0,half_even),(2,half_even),(1,odd),(2,half_even),(0,half_even)):
            state.gate_split_(gate, (i,i+1), cutoff=1e-13, cutoff_mode="abs")
        for site in range(sites):
            state.gate_(half_site, site, contract=True)
    order = [*range(0,2*sites,2), *range(1,2*sites,2)]
    mpdo = state.to_dense().reshape([2]*(2*sites)).transpose(order).reshape(rho0.shape)
    qutip_error = float(np.max(np.abs(qutip_rho-reference)))
    mpdo_error = float(np.max(np.abs(mpdo-reference)))
    excited = np.diag([0.,1.])
    populations = [float(np.trace(mpdo @ place(excited,i)).real) for i in range(sites)]
    assert qutip_error < 1e-8, qutip_error
    trace_error = float(abs(np.trace(mpdo)-1))
    assert mpdo_error < 1e-4 and trace_error < 1e-8, (mpdo_error, trace_error)
    print(json.dumps({"sites": sites, "qutip": qt.__version__, "quimb": quimb.__version__,
                      "qutip_dense_max_error": qutip_error, "mpdo_dense_max_error": mpdo_error,
                      "mpdo_max_bond": state.max_bond(), "mpdo_trace_error": trace_error,
                      "excited_populations": populations}, indent=2))


if __name__ == "__main__":
    main()
