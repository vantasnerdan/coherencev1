"""Effective single-isotope cubic well, not a many-body Gamma-mode quantization.
Run: python anharmonic.py computed_inputs.json > anharmonic.json
K2 and quartic coordinate/source meanings are supplied explicitly in the input.
This reduced single-isotope identification is not a bulk phonon calculation.
NumPy only; parity blocks, two small bases, exact projected coordinate powers.
"""
import json
import sys
import numpy as np
from domain import HC, HBAR


def kron3(a, b, c):
    return np.kron(np.kron(a, b), c)


def groups(energies, strengths, count):
    out = []
    for energy, strength in zip(energies, strengths):
        if out and abs(energy-out[-1]['energy_eV']) < 1e-8:
            out[-1]['raw_strength_A2'] += float(strength)
            out[-1]['block_degeneracy'] += 1
        else:
            if len(out) == count:
                break
            out.append(dict(energy_eV=float(energy), raw_strength_A2=float(strength),
                            block_degeneracy=1))
    return out


def spectrum(k, a, b, mass, npar, Jsum):
    c = HC**2/mass  # hbar^2/m, eV A^2
    roots = np.roots([1., 0., -c*k, -(3*a+2*b)*c*c/2])
    Eb = max(z.real for z in roots if abs(z.imag)<1e-10 and z.real>0)
    # Extra states make P Q^4 P exact; (P Q P)^4 is not used.
    nfull = 2*npar+4
    ladder = np.diag(np.sqrt(np.arange(1, nfull)), 1)
    Q = np.sqrt(c/(2*Eb))*(ladder+ladder.T)
    Q2 = Q@Q
    Q4 = Q2@Q2
    h1 = np.diag(Eb*(np.arange(nfull)+.5)) + .5*(k-Eb*Eb/c)*Q2 + a*Q4/4
    I = np.eye(npar)
    ix = [np.arange(p, 2*npar, 2) for p in (0, 1)]

    def block(parities):
        hs = [h1[np.ix_(ix[p], ix[p])] for p in parities]
        qs = [Q2[np.ix_(ix[p], ix[p])] for p in parities]
        H = kron3(hs[0],I,I)+kron3(I,hs[1],I)+kron3(I,I,hs[2])
        H += b/2*(kron3(qs[0],qs[1],I)+kron3(qs[0],I,qs[2])+kron3(I,qs[1],qs[2]))
        return np.linalg.eigh(H)

    e0, v0 = block((0,0,0))
    e1, v1 = block((1,0,0))
    exy, vxy = block((1,1,0))
    e111, v111 = block((1,1,1))
    qeo = Q[np.ix_(ix[0],ix[1])]
    x1 = kron3(qeo,I,I)@v1[:,0]
    y1 = kron3(I,qeo.T,I)@v1[:,0]
    wparallel = abs(v0.T@x1)**2
    wperp = abs(vxy.T@y1)**2
    higher = e0 > e1[0]
    plines = groups(e0[higher], wparallel[higher], 2)
    tlines = groups(exy, wperp, 1)
    E01 = float(e1[0]-e0[0])
    x01sq = float(wparallel[0])
    for line in plines+tlines:
        line['one_to_two_meV'] = 1e3*(line['energy_eV']-e1[0])
        line['spacing_change_meV'] = line['one_to_two_meV']-1e3*E01
    # Hellmann-Feynman response to diagonal spring tensor, with quartics held.
    def q2_mean(vector, parities, axis):
        factors = [I, I, I]
        p = parities[axis]
        factors[axis] = Q2[np.ix_(ix[p], ix[p])]
        return float(vector @ (kron3(*factors) @ vector))
    dK = [.5*(q2_mean(v1[:,0], (1,0,0), axis) -
               q2_mean(v0[:,0], (0,0,0), axis)) for axis in range(3)]
    scalar_spring_sensitivity = k*sum(dK)/E01
    u20=sum(q2_mean(v0[:,0],(0,0,0),axis) for axis in range(3))/3
    thermal_Z=0.
    thermal_r2=0.
    for en, vec, par, degeneracy in (
            (e0,v0,(0,0,0),1),(e1,v1,(1,0,0),3),
            (exy,vxy,(1,1,0),3),(e111,v111,(1,1,1),1)):
        qs=[Q2[np.ix_(ix[p],ix[p])] for p in par]
        R2=kron3(qs[0],I,I)+kron3(I,qs[1],I)+kron3(I,I,qs[2])
        means=np.sum(vec*(R2@vec),axis=0)
        weights=np.exp(-(en-e0[0])/(8.617333262e-5*300.))
        thermal_Z+=degeneracy*weights.sum()
        thermal_r2+=degeneracy*np.dot(weights,means)
    u2300=float(thermal_r2/(3*thermal_Z))
    # Lattice RPA: chi(q,E)=1/[chi_local(E)^(-1)-J(q)].
    # Keep ALL finite-basis dipole poles, not just a normalized first line.
    ground_dipole = kron3(qeo.T,I,I) @ v0[:,0]
    dipole_weights = abs(v1.T @ ground_dipole)**2
    gaps = e1-e0[0]
    numerator = 2*gaps*dipole_weights
    def chi(z):
        return np.sum(numerator/(gaps*gaps-z))
    chi0 = float(chi(0.))
    denominator = 1-Jsum*chi0
    collective = None
    residue = None
    if denominator>0 and Jsum>0:
        lo, hi = 0., E01*E01
        for _ in range(64):
            mid=(lo+hi)/2
            if Jsum*chi(mid)>1:
                hi=mid
            else:
                lo=mid
        collective=np.sqrt((lo+hi)/2)
        prime=np.sum(numerator/(gaps*gaps-collective**2)**2)
        residue=float(chi(collective**2)**2/prime)
    return dict(parity_basis_per_axis=npar, largest_block=npar**3,
        basis_quantum_eV=float(Eb), ground_energy_eV=float(e0[0]),
        zero_to_one_meV=1e3*E01, zero_to_one_THz=E01/(2*np.pi*HBAR*1e12),
        parametric_stiffness_clock_THz=E01/(np.pi*HBAR*1e12),
        raw_zero_to_one_x_strength_A2=x01sq,
        zero_to_one_oscillator_strength=2*E01*x01sq/c,
        single_pole_effective_spring_eV_A2=E01/(2*x01sq),
        dE01_dKx_A2=dK[0], dE01_dKy_A2=dK[1], dE01_dKz_A2=dK[2],
        dlnE01_dlnK_fixed_quartic=float(scalar_spring_sensitivity),
        single_pole_mass_over_bare=c/(2*E01*x01sq),
        parallel_one_to_two_lines=plines, perpendicular_one_to_two_lines=tlines,
        total_parallel_upward_strength_A2=float(wparallel[higher].sum()),
        reported_parallel_strength_A2=float(sum(x['raw_strength_A2'] for x in plines)),
        local_static_susceptibility_A2_eV=chi0,
        RPA_static_denominator=float(denominator),
        RPA_critical_J_eV_A2=1/chi0,
        RPA_Gamma_quantum_eV=float(collective) if collective is not None else None,
        RPA_Gamma_residue_eV_A2=residue,
        RPA_Gamma_x_strength_A2=float(residue/(2*collective)) if collective is not None else None,
        RPA_Gamma_oscillator_strength=float(residue/c) if residue is not None else None,
        local_mean_square_axis_A2_T0=float(u20),
        local_mean_square_axis_A2_T300=u2300,
        cross_manifold_shift_meV=1e3*float(exy[0]-2*e1[0]+e0[0]))


def main():
    data=json.load(open(sys.argv[1]))
    rows=[]
    for source in data['materials']:
        name=source['name']
        mass=source['mass_eVc2']
        k=source['K2_eV_A2']
        a=source['K4_100_eV_A4']
        along111=source['K4_111_eV_A4']
        b=(3*along111-a)/2
        if min(a,along111)<=0:
            raise ValueError('Computed cubic quartic is not bounded; no forced positive saturation')
        small=spectrum(k,a,b,mass,8,source['J_sum_eV_A2'])
        actual=spectrum(k,a,b,mass,10,source['J_sum_eV_A2'])
        rows.append(dict(material=name, K2_eV_A2=k, K4_100_eV_A4=a,
                         K4_111_eV_A4=along111, cross_quartic_eV_A4=b,
                         harmonic_quantum_meV=1e3*np.sqrt(HC**2*k/mass),
                         basis_change_zero_to_one_meV=actual['zero_to_one_meV']-small['zero_to_one_meV'],
                         basis_change_parallel_12_meV=[y['one_to_two_meV']-x['one_to_two_meV']
                             for x,y in zip(small['parallel_one_to_two_lines'],actual['parallel_one_to_two_lines'])],
                         basis_change_cross_shift_meV=actual['cross_manifold_shift_meV']-small['cross_manifold_shift_meV'],
                         spectrum=actual))
    print(json.dumps(dict(
        assumption=data['effective_well_assumption'],
        potential='K2 sum(u_i^2)/2 + A sum(u_i^4)/4 + B sum_i<j(u_i^2 u_j^2)/2; B=(3 K4_111-A)/2',
        coordinate_meaning='physical single-isotope displacement in A; mass in eV/c^2',
        rows=rows, Pd_H_over_D=rows[1]['spectrum']['zero_to_one_meV']/rows[2]['spectrum']['zero_to_one_meV']),
        indent=2,allow_nan=False))


if __name__=='__main__':
    main()
