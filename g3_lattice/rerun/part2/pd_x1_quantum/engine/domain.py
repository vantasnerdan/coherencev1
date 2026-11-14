"""Constructive displacement mean field: 5 switches, 3 measured host references.
Units: E0=INS reference energy, q=u/sqrt(hbar^2/(m E0)), t=E0*tSI/hbar.
Hamiltonian retains (p+g A)^2/2; screened dipoles, positive quartic,
passive electron response with static counterterm, drive/loss, gradient.
Decision: which sets bound and order; static versus drive-locked is explicit.
Run: python domain.py > domain.json. NumPy only; expected <2s/one thread.
"""
import json
import numpy as np

HBAR = 6.582119569e-16  # eV s
HC = 1973.269804       # eV A
E2 = 4*np.pi/137.035999084
MH = 938.27208943e6
ME = .51099895069e6
A0 = .529177210544
KB = 8.617333262e-5
NAMES = ['dipole', 'quartic', 'electron', 'drive_loss', 'gradient']
MATERIALS = [
    dict(name='NiH', a=3.732, xref=1., xmax=1., E=.089, mass=MH,
         exclusion_A=3.732/np.sqrt(2), source='p0_scales entry101; optic peak, low-T NiH'),
    dict(name='TiH1.84', a=4.454, xref=1.84, xmax=2., E=.1494, mass=MH,
         exclusion_A=4.454/2, source='g0_scales entry313; fitted P=0 INS intercept; lattice from separate near-stoichiometric sample'),
    dict(name='PdD', a=4.084, xref=1., xmax=1., E=.037, mass=1875.612945e6,
         exclusion_A=4.084/np.sqrt(2), source='p0_scales entry101; optic peak, low-T PdD'),
]
DEFAULT = dict(beta=1., electron_strength=1., gamma=.01, force=.001,
               gradient=.1, photon_ratio=10., electron_count=1., T_K=300.)

def reference(m, p):
    nmetal = 4/m['a']**3  # A^-3; loading is isotope/host atom
    ne = p['electron_count']*nmetal
    kF = (3*np.pi**2*ne)**(1/3)
    lamTF = np.sqrt(np.pi*A0/(4*kF))
    epi1 = np.sqrt(E2*nmetal*HC**3/m['mass'])
    epe = np.sqrt(E2*ne*HC**3/ME)
    return nmetal, lamTF, epi1, epe/m['E']

def screening(d, lam):
    y = d/lam
    return (1+y)*np.exp(-y)

def parameters(m, p):
    nm, tf, ep1, we = reference(m, p)
    x = p.get('loading', m['xref'])
    lam = p.get('screen_A', tf)
    wc = p['photon_ratio']
    chi = p['electron_strength']
    # Match FULL disordered proton pole to measured E0 at xref and TF screening.
    # Thus measured interaction/electronic curvature is not softened twice.
    gref2 = (ep1/m['E'])**2*m['xref']
    Jref = gref2*screening(m['exclusion_A'], tf)/3
    rref = 1+gref2/(wc*wc-1)+chi/(we*we-1)
    onsite = rref+Jref
    Jperx = (ep1/m['E'])**2*screening(m['exclusion_A'], lam)/3
    return dict(nm=nm, tf=tf, ep1=ep1, we=we, x=x, lam=lam,
                onsite=onsite, Jperx=Jperx, g=np.sqrt(x)*ep1/m['E'])

def canonical(k, g, wc, chi, we, gamma):
    # z=(q,A,y,p,PiA,py); electron potential .5*we^2*(y-sqrt(chi)*q/we)^2.
    K = np.diag([k+chi, wc*wc+g*g, we*we, 1., 1., 1.])
    K[1,3] = K[3,1] = g
    K[0,2] = K[2,0] = -we*np.sqrt(chi)
    J = np.block([[np.zeros((3,3)), np.eye(3)],
                  [-np.eye(3), np.zeros((3,3))]])
    L = J@K
    # Physical viscous forces: -gamma*qdot, -gamma*Adot, -gamma*ydot.
    L[3,3] -= gamma
    L[3,1] -= gamma*g
    L[4,4] -= gamma
    L[5,5] -= gamma
    return K, L

def impedance(nu, k, g, wc, chi, we, gamma):
    da = wc*wc-nu*nu-1j*gamma*nu
    de = we*we-nu*nu-1j*gamma*nu
    return k-nu*nu-g*g*nu*nu/da+chi-chi*we*we/de-1j*gamma*nu

def evaluate(m, mask=31, settings=None):
    p = DEFAULT | (settings or {})
    on = dict(zip(NAMES, [bool(mask & (1<<i)) for i in range(5)]))
    z = parameters(m, p)
    r = z['onsite']-(z['Jperx']*z['x'] if on['dipole'] else 0.)
    beta = p['beta'] if on['quartic'] else 0.
    chi = p['electron_strength'] if on['electron'] else 0.
    gamma = p['gamma'] if on['drive_loss'] else 0.
    force = p['force'] if on['drive_loss'] else 0.
    c = p['gradient']*m['exclusion_A']**2 if on['gradient'] else 0.
    xc = z['onsite']/z['Jperx'] if on['dipole'] else np.inf
    ans = dict(mask=mask, r=float(r), loading=float(z['x']),
        threshold_loading=float(xc) if np.isfinite(xc) else None,
        threshold_density_m3=float(xc*z['nm']*1e30) if np.isfinite(xc) else None)
    if r <= 0 and beta == 0:
        return ans | dict(state='unbounded_quadratic', bounded=False)
    q0 = np.sqrt(-r/beta) if r < 0 else 0.
    k = -2*r if r < 0 else r
    wc, we, g = p['photon_ratio'], z['we'], z['g']
    K, L = canonical(k, g, wc, chi, we, gamma)
    ev = np.linalg.eigvals(L)
    poles = sorted([v for v in ev if v.imag>1e-8], key=lambda v:v.imag)
    pole = poles[0]
    nu = float(pole.imag)
    A = 0.
    en = -r*r/(4*beta) if r < 0 else 0.
    if force:
        D = impedance(nu, k, g, wc, chi, we, gamma)
        # First-harmonic Duffing: s[(ReD+3 beta*s/4)^2+ImD^2]=f^2.
        co = [.75**2*beta**2, 1.5*beta*D.real, abs(D)**2, -force**2]
        roots = np.roots(np.trim_zeros(co, 'f'))
        stable = []
        for root in roots:
            if abs(root.imag)<1e-9 and root.real>0:
                s = float(root.real)
                slope = (D.real+.75*beta*s)**2+D.imag**2+1.5*beta*s*(D.real+.75*beta*s)
                if slope>0:
                    stable.append(s)
        s = min(stable)
        A = np.sqrt(s)
        aph = 1j*nu*g*A/(wc*wc-nu*nu-1j*gamma*nu)
        yph = we*np.sqrt(chi)*A/(we*we-nu*nu-1j*gamma*nu)
        en += (nu*nu+k)*s/4+3*beta*s*s/32
        en += (nu*nu+wc*wc)*abs(aph)**2/4
        en += nu*nu*abs(yph)**2/4+we*we*abs(yph-np.sqrt(chi)*A/we)**2/4
        diameter = 2*np.sqrt(2*c/(-D.imag)) if c else None
        state = 'static_order_plus_forced_response' if q0 else 'drive_locked_response'
    else:
        diameter = None  # Positive gradient does not select a finite bulk domain.
        state = 'static_order' if q0 else 'disordered'
    uunit = HC/np.sqrt(m['mass']*m['E'])
    life = -HBAR/(m['E']*pole.real) if pole.real < -1e-14 else None
    thermal_q2 = 1/(2*nu*np.tanh(m['E']*nu/(2*KB*p['T_K'])))
    nu_THz = float(nu*m['E']/(2*np.pi*HBAR)/1e12)
    ans |= dict(state=state, bounded=True, amplitude_mode_THz=nu_THz,
        locked_frequency_THz=nu_THz if force else (0. if q0 else None),
        static_displacement_A=float(q0*uunit), driven_amplitude_A=float(A*uunit),
        energy_per_isotope_eV=float(en*m['E']), response_diameter_nm=float(diameter/10) if diameter else None,
        static_bulk_domain='not selected by positive gradient',
        correlation_length_nm=float(np.sqrt(c/k)/10) if c else None,
        ringdown_ps=float(life/1e-12) if life else None,
        ordered_state_lifetime='infinite deterministic persistence' if q0 else ('maintained while driven' if force else 'no ordered state'),
        thermal_rms_A=float(np.sqrt(thermal_q2)*uunit),
        screening_A=float(z['lam']), mode_energy_eV=float(nu*m['E']))
    return ans

if __name__ == '__main__':
    cases = []
    for m in MATERIALS:
        for label, settings in [('TF_reference',{}), ('weak_screen_5A',{'screen_A':5.})]:
            rows = [evaluate(m, mask, settings) for mask in range(32)]
            cases.append(dict(material=m['name'], screening=label, rows=rows))
    output = dict(switch_bits=NAMES, defaults=DEFAULT, materials=MATERIALS,
        selected_mechanism='screened longitudinal dipole plus positive quartic for autonomous static order; drive retained as a probe',
        cases=cases,
        definitions={'threshold':'static mean-field instability, not drive onset; forced response has no finite onset',
          'frequency':'static order locks at zero; amplitude pole and driven pump frequency are separate',
          'size':'localized resonant-response e-fold diameter; static domain size unselected',
          'energy':'stored excess energy; negative only for static broken-symmetry minimum',
          'lifetime':'ringdown amplitude e-fold; ordered minimum persists deterministically'},
        model_assumptions='Lorentz/Yukawa exclusion sphere; TF one-electron host; measured disordered pole calibrated once; beta=1, gamma=.01, force=.001, stiffness=.1, passive electron strength1 and photon/E0=10 are exploratory.')
    print(json.dumps(output, indent=2, allow_nan=False))
