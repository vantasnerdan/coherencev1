"""Minimal binary-collision Monte Carlo for D implanted into Pd, Ti and their hydrides/oxides.

Keeps: ZBL universal potential with exact (Gauss-Mehler) scattering integral, Lindhard-Scharff electronic
stopping (local, continuous, Bragg sum, one scale factor per element), amorphous target with free path N^-1/3
(TRIM convention), NRT/Kinchin-Pease vacancy count per primary recoil (Lindhard partition for damage energy).
Drops: recoil cascades are not followed (NRT counts them), no crystal channelling, no surface binding for the
projectile, no sputtering (taken from literature fits in transport.py).
Units: eV, nm.
"""
import numpy as np

E2 = 1.439964      # e^2/(4 pi eps0), eV nm
A0 = 0.0529177     # Bohr radius, nm
M_D, Z_D = 2.0141, 1

# electronic stopping scale factors relative to Lindhard-Scharff for hydrogen projectiles (estimate, see
# transport.py docstring; tuned only by the range check against the paper's SRIM depths)
# Pd 1.54, Ti 1.52: ICRU Report 49 proton stopping (Geant4 G4BraggModel coefficients) over LS at 5-10 keV/u.
# H, O, Al, Pt, Ag: estimates.
SE_FACTOR = {1: 0.75, 8: 1.0, 22: 1.52, 46: 1.54, 78: 1.4, 13: 1.0, 47: 1.4}


def zbl_phi(x):
    return (0.18175 * np.exp(-3.19980 * x) + 0.50986 * np.exp(-0.94229 * x)
            + 0.28022 * np.exp(-0.40290 * x) + 0.028171 * np.exp(-0.20162 * x))


def zbl_dphi(x):
    return -(0.18175 * 3.19980 * np.exp(-3.19980 * x) + 0.50986 * 0.94229 * np.exp(-0.94229 * x)
             + 0.28022 * 0.40290 * np.exp(-0.40290 * x) + 0.028171 * 0.20162 * np.exp(-0.20162 * x))


def se_ls(E, Z2, Z1=Z_D, M1=M_D):
    """Lindhard-Scharff electronic stopping cross section, eV nm^2, E in eV (lab)."""
    Z = (Z1 ** (2 / 3) + Z2 ** (2 / 3)) ** 1.5
    v_v0 = np.sqrt(np.maximum(E, 0) / (25e3 * M1 / 1.00728))
    return 8 * np.pi * E2 * A0 * Z1 ** (7 / 6) * Z2 / Z * v_v0 * SE_FACTOR.get(Z2, 1.0)


def sn_zbl(E, Z2, M2, Z1=Z_D, M1=M_D):
    """ZBL universal nuclear stopping cross section, eV nm^2."""
    a = 0.8854 * A0 / (Z1 ** 0.23 + Z2 ** 0.23)
    eps = a * E * M2 / (M1 + M2) / (Z1 * Z2 * E2)
    sn = np.where(eps <= 30, np.log(1 + 1.1383 * eps) / (2 * (eps + 0.01321 * eps ** 0.21226 + 0.19593 * eps ** 0.5)),
                  np.log(np.maximum(eps, 1.0)) / (2 * np.maximum(eps, 1e-30)))
    return 4 * np.pi * a * Z1 * Z2 * E2 * M1 / (M1 + M2) * sn


_NQ = 24
_XQ = np.cos((2 * np.arange(1, _NQ // 2 + 1) - 1) * np.pi / (2 * _NQ))


def cm_angle(eps, B):
    """CM scattering angle for ZBL potential, reduced energy eps and reduced impact parameter B (arrays)."""
    R = 1 / (2 * eps) + np.sqrt(1 / (4 * eps ** 2) + B ** 2)        # unscreened turning point (upper bound)
    for _ in range(40):
        x = R
        ph = zbl_phi(x)
        f = 1 - ph / (eps * x) - (B / x) ** 2
        fp = -(zbl_dphi(x) * x - ph) / (eps * x ** 2) + 2 * B ** 2 / x ** 3
        R = np.maximum(R - f / fp, 0.5 * R)
        if np.max(np.abs(f)) < 1e-9:
            break
    u = _XQ[None, :]
    Rr = R[:, None]
    r = Rr / u
    F = 1 - zbl_phi(r) / (eps[:, None] * r) - (B[:, None] * u / Rr) ** 2
    g = np.sqrt(np.clip(1 - u ** 2, 0, None) / np.clip(F, 1e-30, None))
    th = np.pi * (1 - 2 * B / (_NQ * R) * g.sum(1))
    return np.clip(th, 0, np.pi)


def nrt(T, Z2, M2, Ed):
    """NRT displacements for a recoil of energy T (eV)."""
    EL = 30.724 * Z2 * Z2 * np.sqrt(2 * Z2 ** (2 / 3)) * (2 * M2) / M2
    kL = 0.1337 * Z2 ** (1 / 6) * np.sqrt(Z2 / M2)
    e = T / EL
    Td = T / (1 + kL * (3.4008 * e ** (1 / 6) + 0.40244 * e ** 0.75 + e))
    return np.where(T < Ed, 0.0, np.where(Td < 2 * Ed / 0.8, 1.0, 0.8 * Td / (2 * Ed)))


def run(E0, comps, N, nions=4000, Ecut=5.0, zmax=400.0, dz=0.25, seed=1):
    """comps: list of (Z2, M2, atomic fraction, Ed). N: atoms per nm^3. Returns dict of histograms per ion."""
    rng = np.random.default_rng(seed)
    Zs = np.array([c[0] for c in comps], float)
    Ms = np.array([c[1] for c in comps], float)
    fr = np.array([c[2] for c in comps], float)
    fr = fr / fr.sum()
    Eds = np.array([c[3] for c in comps], float)
    L = N ** (-1 / 3)
    pmax = 1 / np.sqrt(np.pi * N * L)
    a = 0.8854 * A0 / (Z_D ** 0.23 + Zs ** 0.23)
    nb = int(zmax / dz)
    stop_h = np.zeros(nb)
    vac_h = np.zeros(nb)
    E = np.full(nions, float(E0))
    z = np.zeros(nions)
    cz = np.ones(nions)
    alive = np.ones(nions, bool)
    nrefl = 0
    Erefl = 0.0
    while alive.any():
        idx = np.nonzero(alive)[0]
        Ei, zi, czi = E[idx], z[idx], cz[idx]
        # electronic loss over the free path (Bragg sum), then move
        Se = sum(f * se_ls(Ei, Zz) for f, Zz in zip(fr, Zs))
        Ei = Ei - N * L * Se
        zi = zi + L * czi
        # collision partner
        k = rng.choice(len(fr), size=idx.size, p=fr)
        M2, Z2, ak = Ms[k], Zs[k], a[k]
        Ec = Ei * M2 / (M_D + M2)
        eps = np.maximum(ak * Ec / (Z_D * Z2 * E2), 1e-6)
        p = pmax * np.sqrt(rng.random(idx.size))
        th = cm_angle(eps, p / ak)
        gam = 4 * M_D * M2 / (M_D + M2) ** 2
        T = gam * Ei * np.sin(th / 2) ** 2
        # vacancies from host recoils
        nv = np.zeros(idx.size)
        for j in range(len(fr)):
            sel = k == j
            if Eds[j] > 0 and sel.any():
                nv[sel] = nrt(T[sel], Zs[j], Ms[j], Eds[j])
        ok = (zi > 0) & (nv > 0)
        np.add.at(vac_h, np.clip((zi[ok] / dz).astype(int), 0, nb - 1), nv[ok])
        Ei = Ei - T
        psi = np.arctan2(np.sin(th), np.cos(th) + M_D / M2)
        phi = 2 * np.pi * rng.random(idx.size)
        czi = czi * np.cos(psi) + np.sqrt(np.clip(1 - czi ** 2, 0, 1)) * np.sin(psi) * np.cos(phi)
        refl = zi < 0
        nrefl += refl.sum()
        Erefl += Ei[refl].sum()
        stopped = (Ei < Ecut) & ~refl
        np.add.at(stop_h, np.clip((zi[stopped] / dz).astype(int), 0, nb - 1), 1)
        E[idx], z[idx], cz[idx] = Ei, zi, czi
        alive[idx[refl | stopped]] = False
    zc = (np.arange(nb) + 0.5) * dz
    imp = stop_h / nions
    Rp = (zc * imp).sum() / imp.sum()
    dRp = np.sqrt(((zc - Rp) ** 2 * imp).sum() / imp.sum())
    return dict(z=zc, imp=imp, vac=vac_h / nions, RN=nrefl / nions, RE=Erefl / nions / E0,
                Rp=Rp, dRp=dRp, nv_per_ion=vac_h.sum() / nions)


if __name__ == "__main__":
    import time
    from scipy.integrate import quad
    # check the quadrature against direct integration of the scattering integral
    for eps, B in [(0.01, 0.5), (0.3, 1.0), (3.0, 0.2), (30.0, 2.0)]:
        th = cm_angle(np.array([eps]), np.array([B]))[0]
        R0 = 1 / (2 * eps) + np.sqrt(1 / (4 * eps ** 2) + B ** 2)
        from scipy.optimize import brentq
        R0 = brentq(lambda r: 1 - zbl_phi(r) / (eps * r) - (B / r) ** 2, 1e-6, R0 + 1e-9)
        I = quad(lambda r: 1 / (r ** 2 * np.sqrt(1 - zbl_phi(r) / (eps * r) - (B / r) ** 2)), R0, np.inf, limit=400)[0]
        print(f"eps={eps} B={B}: GM {th:.5f}  quad {np.pi - 2 * B * I:.5f}")
    t = time.time()
    for name, comps, N in [("Pd", [(46, 106.42, 1, 40)], 68.0), ("Ti", [(22, 47.87, 1, 30)], 56.7)]:
        for E in [500, 2000, 5000, 13000]:
            r = run(E, comps, N, nions=2000)
            print(f"{name} E_lab={E}: Rp={r['Rp']:.1f} dRp={r['dRp']:.1f} nm RN={r['RN']:.2f} vac/ion={r['nv_per_ion']:.2f}")
    print("time", time.time() - t)
