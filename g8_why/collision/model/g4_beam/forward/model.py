"""Forward model of the relative thick-target D-D yield R(E) reported by Karahadian et al.,
Nat. Commun. 17:8845 (2026), Fig. 2, with an explicit measurement layer.

Layers
  reaction   : Bosch & Hale (1992) S-factors for D(d,p)T and D(d,n)3He; screened cross section of
               the paper's Eq. 7; Czerski (2026 v3) thermal-track rate as a separate source.
  stopping   : ICRU-49 (Andersen-Ziegler type) velocity-proportional electronic stopping + ZBL universal
               nuclear stopping, Bragg additivity for MD_x. Path range (CSDA) and mean projected range
               (first-moment transport, sigma_tr = 2 S_n / (gamma E) for M2 >> M1).
  target     : thick-target integral along the slowing-down path, with constant Ue, a step Ue(z) surface
               layer (z = mean projected depth), or a per-beam-energy Ue (paper's sigmoid, SI Table 2).
  beam       : 85 % D+ / 15 % D2+ by current (D2+ -> two D at E/2), 5 % FWHM Gaussian energy spread.
  measurement: Faraday-cup charge normalisation (I or I^2; collected-fraction exponents), proton/neutron
               stitch factor below 4.5 keV, energetic-neutral component at fixed energy.

All energies: E_lab of the deuteron in keV unless named E_cm (E_cm = E_lab / 2 for d on d at rest).
Stopping cross sections: eV / (1e15 host atoms / cm^2) per host atom (compound counted per host atom).
NumPy only. Importable; the fitting driver is fits.py.
"""
import numpy as np

# ----------------------------------------------------------------------------------------------
# Reaction layer
# ----------------------------------------------------------------------------------------------
BG = 31.3970                       # sqrt(keV), Bosch & Hale (1992) Gamow constant for D+D
EG = BG**2                         # keV, Gamow energy as used in Eq. 7 / Eq. 12
BH = {'p': (5.5576e4, 2.1054e2, -3.2638e-2, 1.4987e-6, 1.8181e-10),    # D(d,p)T, keV mb
      'n': (5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9)}    # D(d,n)3He, keV mb


def S_BH(E_cm, ch='p'):
    """Bosch-Hale astrophysical S-factor (keV mb), E_cm in keV."""
    A1, A2, A3, A4, A5 = BH[ch]
    E = np.asarray(E_cm, float)
    return A1 + E*(A2 + E*(A3 + E*(A4 + E*A5)))


def sigma(E_cm, Ue=0.0, ch='p'):
    """Screened cross section, paper Eq. 7 (mb): S(E)/sqrt(E(E+Ue)) exp(-BG/sqrt(E+Ue))."""
    E = np.asarray(E_cm, float)
    return S_BH(E, ch)/np.sqrt(E*(E + Ue))*np.exp(-BG/np.sqrt(E + Ue))


def reaction_S(E_cm, ch='p', reaction='proton'):
    """Common fixed nuclear background, source v3 literal Eqs14/17, keV mb.

    Gamma_p=40meV, phase115deg, a=7fm, E_R->0. Neutron parameters
    are unknown: baseline 'proton' leaves them Bosch-Hale. 'mirror' is
    the explicitly authorized fixed-neutron sensitivity; 'eq7' is effective.
    """
    if reaction not in ('proton', 'mirror', 'eq7'):
        raise ValueError(reaction)
    sb = S_BH(E_cm, ch)
    if reaction == 'eq7' or (ch == 'n' and reaction == 'proton'):
        return sb
    e = np.asarray(E_cm, float)
    sr = np.sqrt(2*986.)*197326.9804**3*(40.e-6)/(7.*937806.5**1.5*e**2)*10.
    return sb + sr + 2*np.sqrt(sb*sr/3)*np.cos(np.deg2rad(115.))


def reaction_sigma(E_cm, Ue=0., ch='p', reaction='proton'):
    """Eq7 penetration with the fixed common reaction S-factor, in mb."""
    e = np.asarray(E_cm, float)
    return reaction_S(e, ch, reaction)/np.sqrt(e*(e+Ue))*np.exp(-BG/np.sqrt(e+Ue))




# ----------------------------------------------------------------------------------------------
# Stopping layer
# ----------------------------------------------------------------------------------------------
M_D, Z_D = 2.01410, 1
ELEM = {'Pd': dict(Z=46, M=106.42, A1=5.238, n0=6.80e22, dV=0.159),   # dV: fractional volume per D/M
        'Ti': dict(Z=22, M=47.867, A1=4.858, n0=5.67e22, dV=0.142),
        'H':  dict(Z=1,  M=2.0141, A1=1.254)}
# A1: ICRU Report 49 (1993) proton coefficients, S_e = A1 sqrt(T[keV/amu]) for T < 10 keV/amu
# (values as tabulated in Geant4 G4BraggModel::ElectronicStoppingPower). dV: lattice expansion,
# beta-PdD0.65 a=4.02 A vs 3.89 A; delta-TiD1.7 a=4.44 A (4 Ti/cell) vs alpha-Ti 5.67e22 cm^-3.

MATERIALS = {   # host, near-surface D/M used in the stopping (assumed; see report)
    'Pd':  dict(host='Pd', x=0.10),
    'PdD': dict(host='Pd', x=0.70),
    'Ti':  dict(host='Ti', x=1.00),
    'TiD': dict(host='Ti', x=1.70),
}


def Se_atom(E_lab, el):
    """Electronic stopping of D (lab energy keV) on one atom of element el, eV/(1e15 at/cm^2)."""
    T = np.asarray(E_lab, float)*1.007276/M_D        # keV/amu (proton-equivalent energy)
    return ELEM[el]['A1']*np.sqrt(T)


def Sn_atom(E_lab, Z2, M2, Z1=Z_D, M1=M_D):
    """ZBL universal nuclear stopping (Ziegler, Biersack, Littmark 1985), eV/(1e15 at/cm^2)."""
    E = np.asarray(E_lab, float)
    a = Z1**0.23 + Z2**0.23
    eps = 32.53*M2*E/(Z1*Z2*(M1 + M2)*a)
    with np.errstate(divide='ignore', invalid='ignore'):
        sn = np.where(eps <= 30,
                      np.log1p(1.1383*eps)/(2*(eps + 0.01321*eps**0.21226 + 0.19593*np.sqrt(eps))),
                      np.log(eps)/(2*eps))
    return 8.462*Z1*Z2*M1*sn/((M1 + M2)*a)


def stopping(E_lab, mat):
    """(S_e, S_n) per host atom of M D_x by Bragg additivity, eV/(1e15 host atoms/cm^2)."""
    m = MATERIALS[mat]; h = ELEM[m['host']]; x = m['x']
    se = Se_atom(E_lab, m['host']) + x*Se_atom(E_lab, 'H')
    sn = Sn_atom(E_lab, h['Z'], h['M']) + x*Sn_atom(E_lab, 1, M_D)
    return se, sn


def host_density(mat):
    m = MATERIALS[mat]; h = ELEM[m['host']]
    return h['n0']/(1 + h['dV']*m['x'])              # host atoms / cm^3


# ----------------------------------------------------------------------------------------------
# Target: grids and thick-target integrals
# ----------------------------------------------------------------------------------------------
class Target:
    """Precomputes stopping, path and projected-depth functions on a log grid of lab energies.
    E_cut (lab keV) is the lower end of the slowing-down integral (convention, see report)."""

    def __init__(self, mat, E_cut=0.1, E_max=32.0, N=1600, reaction='proton'):
        self.mat = mat
        self.reaction = reaction
        self.x = MATERIALS[mat]['x']
        self.n = host_density(mat)
        self.E = np.geomspace(E_cut, E_max, N)
        se, sn = stopping(self.E, mat)
        self.S = se + sn                                     # eV/(1e15/cm2) per host atom
        self.Sn_host = Sn_atom(self.E, ELEM[MATERIALS[mat]['host']]['Z'], ELEM[MATERIALS[mat]['host']]['M'])
        self.dEdx = self.S*1e-15*self.n*1e-7                 # eV/nm
        M2 = ELEM[MATERIALS[mat]['host']]['M']
        gam = 4*M_D*M2/(M_D + M2)**2
        sig_tr = 2*self.Sn_host*1e-15/(gam*self.E*1e3)        # cm^2 (M2 >> M1)
        inv_ltr = self.n*sig_tr*1e-7                         # 1/nm
        dE = np.diff(self.E)*1e3                             # eV
        w = 1.0/self.dEdx                                    # nm/eV
        self.path = np.concatenate([[0], np.cumsum(0.5*(w[1:] + w[:-1])*dE)])     # nm from E_cut
        g = inv_ltr*w                                        # 1/eV
        self.Lam = np.concatenate([[0], np.cumsum(0.5*(g[1:] + g[:-1])*dE)])
        q = np.exp(self.Lam)*w
        self.Q = np.concatenate([[0], np.cumsum(0.5*(q[1:] + q[:-1])*dE)])
        self._lnE = np.log(self.E)

    # --- ranges ---
    def path_range(self, E0):
        return np.interp(np.log(E0), self._lnE, self.path)

    def proj_range(self, E0):
        L = np.interp(np.log(E0), self._lnE, self.Lam); Q = np.interp(np.log(E0), self._lnE, self.Q)
        return np.exp(-L)*Q

    def energy_at_depth(self, E0, d):
        """Energy (lab keV) at which an ion started at E0 reaches mean projected depth d (nm)."""
        E0 = np.asarray(E0, float)
        L0 = np.interp(np.log(E0), self._lnE, self.Lam); Q0 = np.interp(np.log(E0), self._lnE, self.Q)
        Qd = Q0 - d*np.exp(L0)
        return np.exp(np.interp(Qd, self.Q, self._lnE, left=np.log(self.E[0]))), Qd <= 0

    # --- cumulative yield integrals ---
    def cumulative(self, Ue, ch='p'):
        """C(E) = int x sigma(E'/2)/S(E') dE': reactions per incident D."""
        f = self.x*reaction_sigma(self.E/2, Ue, ch, self.reaction)/self.S*1e-9
        return np.concatenate([[0], np.cumsum(0.5*(f[1:] + f[:-1])*np.diff(self.E))])

    def Y_const(self, E0, Ue, ch='p', C=None):
        """Thick-target probability per incident D at lab energy E0 (keV)."""
        if C is None:
            C = self.cumulative(Ue, ch)
        return np.interp(np.log(E0), self._lnE, C)

    def Y_rows(self, E0, Ue, ch='p'):
        """Yield per incident D for each starting energy E0[i] with its own constant Ue[i]
        (used for the per-beam-energy Ue_eff of the paper's sigmoid)."""
        E0 = np.asarray(E0, float); shp = E0.shape
        e0 = E0.ravel(); u = np.broadcast_to(np.asarray(Ue, float), shp).ravel()
        Ecm = self.E[None, :]/2
        f = self.x*reaction_sigma(Ecm, u[:, None], ch, self.reaction)/self.S[None, :]*1e-9
        dE = np.diff(self.E)[None, :]
        C = np.concatenate([np.zeros((len(e0), 1)), np.cumsum(0.5*(f[:, 1:] + f[:, :-1])*dE, 1)], 1)
        le = np.log(e0)
        j = np.clip(np.searchsorted(self._lnE, le) - 1, 0, len(self.E) - 2)
        t = (le - self._lnE[j])/(self._lnE[j + 1] - self._lnE[j])
        r = np.arange(len(e0))
        result = (C[r, j]*(1 - t) + C[r, j + 1]*t).reshape(shp)
        return np.where(E0 <= self.E[0], 0.0, result)

    def Y_layer(self, E0, Ue_s, Ue_b, d, ch='p'):
        """Step Ue(z): Ue_s for mean projected depth z < d (nm), Ue_b beyond."""
        Cs = self.cumulative(Ue_s, ch); Cb = self.cumulative(Ue_b, ch)
        Ed, inside = self.energy_at_depth(E0, d)
        lnE0 = np.log(E0); lnEd = np.log(Ed)
        y = (np.interp(lnE0, self._lnE, Cs) - np.interp(lnEd, self._lnE, Cs) + np.interp(lnEd, self._lnE, Cb))
        return np.where(inside, np.interp(lnE0, self._lnE, Cs), y)


# ----------------------------------------------------------------------------------------------
# Beam layer
# ----------------------------------------------------------------------------------------------
F_D1, F_D2 = 0.85, 0.15            # current fractions D+, D2+
SPREAD_FWHM = 0.05
_gh_x, _gh_w = np.polynomial.hermite_e.hermegauss(5)
_gh_w = _gh_w/_gh_w.sum()


def beam_average(yfun, E0):
    """Per unit charge: 0.85 Y(E0) + 2*0.15 Y(E0/2), averaged over a Gaussian spread (5 % FWHM).
    yfun maps an array of lab energies to yield per incident deuteron."""
    E0 = np.atleast_1d(np.asarray(E0, float))
    sg = SPREAD_FWHM/2.3548
    Es = E0[:, None]*(1 + sg*_gh_x[None, :])
    y = F_D1*yfun(Es) + 2*F_D2*yfun(Es/2)
    return (y*_gh_w[None, :]).sum(1)


# ----------------------------------------------------------------------------------------------
# Measurement layer
# ----------------------------------------------------------------------------------------------
# SI Fig. 3 Faraday-cup IV (V_lab kV : I uA), constant extraction voltage
IV_V = np.arange(0.5, 15.01, 0.5)
IV_I = np.array([38, 68, 94, 110, 118, 123, 130, 134, 136, 138, 140, 145, 152, 156, 160, 165, 171, 176, 180,
                 182, 192, 207, 217, 222, 228, 238, 251, 263, 275, 283], float)


def I_cup(E_cm):
    return np.interp(2*np.asarray(E_cm, float), IV_V, IV_I)


I_INT = I_cup(2.5)     # 138 uA, top of the internal-bias regime (5 kV)
I_EXT0 = I_cup(2.75)   # 140 uA, onset of external bias (5.5 kV)


def collected_fraction(E_cm, b_hi=0.0, b_lo=0.0):
    """c = I_cup / I_true.  b_hi = 1: true current constant above 5.5 kV (SI Fig 3 caption, focusing);
    b_lo = 1: true current constant below 5 kV (falling IV is a cup-collection effect)."""
    E = np.asarray(E_cm, float); I = I_cup(E)
    return np.where(E >= 2.75, (I/I_EXT0)**b_hi, (I/I_INT)**b_lo)


def plotted(E_cm, Yion, Yneu, meas):
    """Per-energy quantity that the figure is proportional to (before normalisation).
    rate = I_true(E) Yion(E) + Phi_n Yneu,  I_true = I_cup / c(E),
    Phi_n = f_n I_true (delta=0, neutrals made from the delivered beam) or f_n I_INT (delta=1,
    neutrals made from the constant extracted current), plotted = rate / I_cup^(1+k).
    Yion: per-charge beam yield (beam_average); Yneu: per-D yield at the neutral energy E_n."""
    E = np.asarray(E_cm, float)
    k = meas.get('k', 0); bh = meas.get('b_hi', 0.0); bl = meas.get('b_lo', 0.0)
    fn = meas.get('f_n', 0.0); dl = meas.get('delta', 0)
    Ic = I_cup(E); It = Ic/collected_fraction(E, bh, bl)
    phi = fn*(It if dl == 0 else I_INT)
    return (It*Yion + phi*Yneu)/Ic**(1 + k)


def yield_function(T, physics, p):
    """Per-incident-D physical probability; only the confirmed candidates.

    'const' is the no-plateau control and the ion part of fast neutrals.
    Thermal sources are supplied per incident charge to relative_yield.
    """
    if physics == 'const':
        C = {ch: T.cumulative(p['Ue'], ch) for ch in ('p', 'n')}
        return lambda E, ch: np.interp(np.log(E), T._lnE, C[ch], left=0.0)
    if physics == 'layer':
        Cs = {ch: T.cumulative(p['Ue_s'], ch) for ch in ('p', 'n')}
        Cb = {ch: T.cumulative(p['Ue'], ch) for ch in ('p', 'n')}

        def yf(E, ch):
            Ed, inside = T.energy_at_depth(E, p['d'])
            l0, ld = np.log(E), np.log(Ed)
            ys = np.interp(l0, T._lnE, Cs[ch])
            return np.where(inside, ys, ys - np.interp(ld, T._lnE, Cs[ch]) + np.interp(ld, T._lnE, Cb[ch]))
        return yf
    if physics == 'sigmoid_tt':
        return lambda E, ch: T.Y_rows(E, Ue_sigmoid(np.asarray(E)/2, p['Ue_s'], p['Ue'], p['E0']), ch)
    raise ValueError(physics)


def relative_yield(T, E_cm, physics, p, meas):
    """Model of the plotted relative yield R(E) = Y(E)/Y(6.5 keV) for one material.
    Proton channel for E <= 4.0 keV, stitched to the neutron channel at 4.5 keV and multiplied by
    the stitch factor s; neutron channel for E >= 4.5 keV.  meas may hold k, b_hi, b_lo, f_n,
    E_n (lab keV), delta and s. No arbitrary additive relative-yield floor."""
    E = np.asarray(E_cm, float)
    yf = yield_function(T, 'const' if physics == 'thermal' else physics, p)
    def Yb(EE, ch):
        direct = beam_average(lambda x: yf(x, ch), 2*np.atleast_1d(EE))
        if physics == 'thermal' and ch == 'p':
            direct = direct + p['thermal_per_charge'](np.atleast_1d(EE))
        return direct
    En = meas.get('E_n', 5.0)
    yn = {ch: float(np.atleast_1d(yf(np.array([En]), ch))[0]) for ch in ('p', 'n')}
    pl = lambda EE, ch: plotted(EE, Yb(EE, ch), yn[ch], meas)
    ref = pl(np.array([6.5]), 'n')[0]
    stitch = pl(np.array([4.5]), 'n')[0]/pl(np.array([4.5]), 'p')[0]
    R = np.where(E <= 4.0, pl(E, 'p')*stitch*meas.get('s', 1.0), pl(E, 'n'))/ref
    return R


# ----------------------------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------------------------
MATS = ['Pd', 'PdD', 'Ti', 'TiD']


def load_data(path):
    """Returns {mat: (E_cm, R, sigma_lnR)} for E < 6.5 from errbars.csv."""
    out = {m: ([], [], []) for m in MATS}
    for l in open(path):
        if l.startswith('#') or l.startswith('E_cm'):
            continue
        E, m, R, h, sl = l.split(',')[:5]
        if float(E) >= 6.5:
            continue
        out[m][0].append(float(E)); out[m][1].append(float(R)); out[m][2].append(float(sl))
    return {m: tuple(np.array(v) for v in out[m]) for m in MATS}


# SI Table 1: S:B with background 0.007 cph; I in uA
SI1_E = np.array([3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.25])
SI1_I = np.array([145, 138, 134, 123, 110, 68, 38.])
SI1_SB = {'Pd': [443, 121, 83, 53, 44, 15, 21], 'PdD': [465, 177, 138, 92, 43, 32, 30],
          'Ti': [1525, 225, 82, 106, 69, 47, 15], 'TiD': [5427, 1082, 368, 261, 112, 62, 56]}
SI2 = {'Pd': (3.28, 0.8, 0.6), 'PdD': (3.68, 1.7, 0.6), 'Ti': (3.21, 0.3, 0.6), 'TiD': (2.75, 0.9, 0.6)}


def counts_per_C_table(mat):
    """Counts per coulomb as defined in SI Table 1 (rate = yield/C x I x 3600)."""
    return np.array(SI1_SB[mat], float)*0.007/(SI1_I*1e-6*3600)


def Ue_sigmoid(E_cm, Ue_s, Ue_b, E0, Delta=0.5):
    """Paper Eq. 16."""
    return Ue_b + (Ue_s - Ue_b)/(1 + np.exp((np.asarray(E_cm) - E0)/Delta))
