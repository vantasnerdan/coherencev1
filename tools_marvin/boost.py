"""Lowest-order positive-energy Dirac term, DOI 10.1088/1361-6455/acf3be Eq. 27."""
import json
import numpy as np


def projected_boost(energies_MeV, spin_cross_internal_pc_MeV, P_over_Mc, mass_MeV):
    """Array H[f,i]; components of Σ_j σ_j×(π_j c) in the SAME eigenbasis.

    P/(M c) can be a lattice matrix element, not a free-particle transition rate.
    mass_MeV is the equal nucleon rest energy, not the total nuclear rest energy.
    """
    energies = np.asarray(energies_MeV)
    projected = np.tensordot(P_over_Mc, spin_cross_internal_pc_MeV, axes=1)
    return -1j * (energies[None,:]-energies[:,None]) * projected / (2*mass_MeV)


def main():
    # Two-state HO witness: |0s,T0spin> -> |1px,S0spin>, fixed T=0 isospin.
    # The transverse orbital projections of p are zero in this subspace.
    m, quantum = 938.92, 20.0
    spin_x = np.array([[0,1],[1,0]])
    spin_y = np.array([[0,-1j],[1j,0]])
    spin_z = np.diag([1,-1])
    unit = np.eye(2)
    triplet = np.array([0,1,1,0])/np.sqrt(2)
    singlet = np.array([0,1,-1,0])/np.sqrt(2)
    p0 = np.sqrt((m/2)*quantum/2)
    p_x = p0*np.array([[0,-1j],[1j,0]])
    states = np.stack([np.kron(unit[:,0],triplet), np.kron(unit[:,1],singlet)],axis=1)
    spin_diff = [np.kron(s,unit)-np.kron(unit,s) for s in (spin_x,spin_y,spin_z)]
    tensors = (np.zeros((8,8)), np.kron(p_x,spin_diff[2]), -np.kron(p_x,spin_diff[1]))
    o = np.array([states.conj().T @ t @ states for t in tensors])
    velocity = 1e-4  # declared dimensionless matrix element, not a material input
    h = projected_boost([0,quantum],o,[0,velocity,0],m)
    target = -quantum*p0*velocity/m
    assert np.allclose(h, [[0,target],[target,0]], rtol=1e-14, atol=1e-16), h
    print(json.dumps({"model": "analytic two-nucleon HO test; not a physical matrix element",
                      "off_diagonal_MeV": float(h[1,0].real), "reference_MeV": target,
                      "P_over_Mc": velocity}, indent=2))


if __name__ == "__main__":
    main()
