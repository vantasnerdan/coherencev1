# g1_rank/rank.py (replacement with same-mode metric)
"""Selected screened-dipole + quartic mechanism, using domain.py unchanged.
Local: python rank.py > rank.json. Published node directories:
PYTHONPATH=../g0_model python rank.py > rank.json
Decision: which measurable/computable input changes the five domain outcomes?
Expected <5s, one BLAS/OpenMP thread, <50MB. No fit or regression.
"""
import json
import numpy as np
from domain import MATERIALS, DEFAULT, MH, HBAR, evaluate, parameters, canonical, KB

# Completed g0_scales entry330, preserving its Ti mixed-sample descriptor.
MATS = [dict(m) for m in MATERIALS]
MATS[1].update(name='delta-TiH2_reference', xref=2., E=.149,
              source='entry330: TiH2 a4.454A; rounded149meV INS descriptor from TiH1.84')
FIELDS = ['threshold_loading', 'threshold_density_m3', 'amplitude_mode_THz',
          'static_displacement_A', 'energy_per_isotope_eV',
          'correlation_length_nm', 'ringdown_ps']
BASE = dict(screen_A=5., force=0.)

def ranges(m):
    ni, ti = m['name']=='NiH', m['name']=='delta-TiH2_reference'
    return {
      'screen_A': (.2, 5., 'electronic-structure-computable', 'exploratory finite-q screening proxy'),
      'beta': (.1, 10., 'electronic-structure-computable', 'exploratory quartic'),
      'gamma': (.001, .1, 'lookup', 'exploratory; actual mode linewidth/ringdown needed'),
      'gradient': (.01, 1., 'electronic-structure-computable', 'exploratory dispersion stiffness'),
      'electron_strength': (0., 4., 'electronic-structure-computable', 'exploratory passive response strength'),
      'photon_ratio': (5., 20., 'neither', 'exploratory EM boundary-mode choice'),
      'electron_count': (.5, 2., 'electronic-structure-computable', 'exploratory one-band effective count'),
      'loading': ((.01 if ni else 1.74 if ti else .6), m['xmax'], 'lookup',
                  'Ni phase endpoints, Ti entry313 range, Pd entry101 scenario; fixed-cell interpolation'),
      'a': ((3.726 if ni else 4.453 if ti else 4.080),
            (3.738 if ni else 4.455 if ti else 4.088), 'lookup',
            'Ni/Ti reported lattice precision; Pd reference-cell spread scenario'),
      'E': ((.089 if ni else .1493 if ti else .037),
            (.130 if ni else .1495 if ti else .065), 'lookup',
            'Ni/Pd peak-to-shoulder MODE CHOICE, not uncertainty; Ti fitted-intercept interval'),
      'T_K': (5., 600., 'lookup', 'temperature scenario; fixed force/loss coefficients'),
    }

def collect(m, key, value):
    mm, pp = dict(m), BASE.copy()
    if key in ('a', 'E'):
        mm[key] = value
        if key=='a':
            mm['exclusion_A'] *= value/m['a']
    else:
        pp[key] = value
    return evaluate(mm, 31, pp)

def spread(rows, key):
    vals = [abs(r[key]) for r in rows if r.get(key) is not None]
    if not vals:
        return {'decades':None, 'reason':'not defined'}
    if any(v==0 for v in vals):
        return {'decades':0. if all(v==0 for v in vals) else None,
                'reason':'identically zero' if all(v==0 for v in vals) else 'zero/phase boundary; no log floor'}
    return {'decades':float(np.log10(max(vals)/min(vals))),
            'min_abs':float(min(vals)), 'max_abs':float(max(vals))}

def critical_screen(m):
    def margin(lam):
        z = parameters(m, DEFAULT | {'screen_A':lam})
        return z['onsite']-z['Jperx']*m['xmax']
    lo, hi = .01, 1e6
    if margin(hi)>=0:
        return {'material':m['name'], 'reachable_loading':m['xmax'],
                'lambda_critical_A':None, 'reason':'even unscreened coupling insufficient'}
    for _ in range(80):
        mid = (lo+hi)/2
        if margin(mid)>0:
            lo = mid
        else:
            hi = mid
    return {'material':m['name'], 'reachable_loading':m['xmax'],
            'lambda_critical_A':float((lo+hi)/2)}

def gain_branch(m, changes=None):
    # dn/dt=(G-kappa)n-s*n^2; phase rotates with the optic mode.
    # Pump conversion cross sections are explicit UNKNOWN coefficients.
    p = dict(gamma=.01, gradient=.1, saturation_q=.1, pump_ratio=2.,
             sigma_e=1e-20, sigma_H=1e-20, sigma_Q=1e-20,
             thermal_W_mK=10., diffusion_m2_s=1e-7, T_K=300.,
             qd=np.pi/2) | (changes or {})
    base = DEFAULT | {'force':0., 'gamma':p['gamma']}
    z = parameters(m, base)
    d = m['exclusion_A']*1e-10
    q = p['qd']/d
    r = z['onsite']-z['Jperx']*z['x']
    kq = r+4*p['gradient']*np.sin(p['qd']/2)**2
    _, L = canonical(kq,z['g'],base['photon_ratio'],
                     base['electron_strength'],z['we'],p['gamma'])
    pole = min([v for v in np.linalg.eigvals(L) if v.imag>1e-8],
               key=lambda v:v.imag)
    Omega0 = m['E']/HBAR
    Omega, kappa = pole.imag*Omega0, -2*pole.real*Omega0
    # Differentiate the SAME hybrid-mode dispersion D(nu,q)=0.
    # M_dyn=-dD/d(nu^2); it also normalizes the mode action.
    nu2 = pole.imag**2
    Mdyn = (1+z['g']**2*base['photon_ratio']**2/
            (base['photon_ratio']**2-nu2)**2+
            base['electron_strength']*z['we']**2/(z['we']**2-nu2)**2)
    vg = Omega0*p['gradient']*d*np.sin(p['qd'])/(pole.imag*Mdyn)
    vphase = Omega/q
    ne = z['nm']*1e30  # one-electron host scenario, NOT measured n_e
    charge = 1.602176634e-19
    Eph_J = HBAR*Omega*charge
    # G_e=sigma_e*(j/e-ne*vphase), G_H=sigma_H*Phi_H,
    # G_Q=sigma_Q*kthermal*|gradT|/Eph.
    jkin = charge*ne*vphase
    jth = jkin+charge*kappa/p['sigma_e']
    Hth = kappa/p['sigma_H']
    gradT = kappa*Eph_J/(p['sigma_Q']*p['thermal_W_mK'])
    gradx = Hth/(p['diffusion_m2_s']*z['nm']*1e30)
    mu = p['pump_ratio']-1
    ans = dict(material=m['name'], q_m_inv=float(q),
        mode_frequency_THz=float(Omega/(2*np.pi*1e12)),
        phase_velocity_m_s=float(vphase), group_velocity_m_s=float(vg),
        passive_energy_decay_s_inv=float(kappa),
        current_kinematic_A_m2=float(jkin), current_threshold_A_m2=float(jth),
        H_flux_threshold_m2_s=float(Hth), temperature_gradient_threshold_K_m=float(gradT),
        loading_gradient_threshold_m_inv=float(gradx),
        after_pump_amplitude_decay_ps=float(2/kappa/1e-12),
        pump_ratio=p['pump_ratio'], coefficients=p)
    if mu<=0:
        return ans | dict(state='below_or_at_gain_threshold', gain_length_nm=None,
                          phonon_number=0., linewidth_FWHM_GHz=None, locked_frequency_THz=None)
    length = 2*vg/(kappa*mu)  # amplitude spatial e-fold length
    Nc = z['nm']*z['x']*1e30*np.pi*length**3/6
    nstar = .5*Nc*pole.imag*Mdyn*p['saturation_q']**2*mu
    nth = 1/np.expm1(HBAR*Omega/(KB*p['T_K']))
    G = p['pump_ratio']*kappa
    noise = kappa*(nth+.5)/2+G/4
    # Duffing amplitude-to-phase conversion, with nonlinear damping saturation.
    henry = 3*base['beta']*p['saturation_q']**2*Omega0/(4*kappa*pole.imag*Mdyn)
    fwhm = noise*(1+henry**2)/(2*np.pi*nstar)
    # Analytic logistic growth from nstar/100, sampled over 10 gain times.
    times = np.array([0.,1.,5.,10.])/(kappa*mu)
    growth = nstar/(1+99*np.exp(-kappa*mu*times))
    shift = 3*base['beta']*Omega0*p['saturation_q']**2*mu/(8*pole.imag*Mdyn)
    stop_t = np.array([0.,1.,2.,5.])*2/kappa
    stop_n = nstar*np.exp(-kappa*stop_t)/(1+mu*(1-np.exp(-kappa*stop_t)))
    return ans | dict(state='self_sustained_limit_cycle',
        gain_length_nm=float(length/1e-9), coherent_isotopes=float(Nc),
        phonon_number=float(nstar), hybrid_action_metric=float(Mdyn),
        linewidth_FWHM_GHz=float(fwhm/1e9) if nstar>=10 else None,
        linewidth_regime='large-occupation phase diffusion' if nstar>=10 else 'occupation<10; no semiclassical linewidth assigned',
        locked_frequency_THz=float((Omega+shift)/(2*np.pi*1e12)),
        nonlinear_frequency_shift_GHz=float(shift/(2*np.pi*1e9)),
        radial_relaxation_ps=float(1/(kappa*mu)/1e-12),
        logistic_time_ps=(times/1e-12).tolist(), logistic_phonons=growth.tolist(),
        after_stop_time_ps=(stop_t/1e-12).tolist(), after_stop_phonons=stop_n.tolist())

def gain_sweep():
    spans = {
      'sigma_e':(1e-24,1e-18,'electronic-structure-computable'),
      'sigma_H':(1e-24,1e-18,'neither'),
      'sigma_Q':(1e-24,1e-18,'neither'),
      'diffusion_m2_s':(1e-12,1e-7,'lookup'),
      'thermal_W_mK':(1.,100.,'lookup'),
      'gamma':(.001,.1,'lookup'),
      'gradient':(.01,1.,'electronic-structure-computable'),
      'saturation_q':(.03,.3,'neither'),
      'pump_ratio':(1.01,4.,'neither'),
      'T_K':(5.,600.,'lookup'),
      'qd':(.2,2.8,'electronic-structure-computable')}
    fields = ['current_threshold_A_m2','H_flux_threshold_m2_s',
              'temperature_gradient_threshold_K_m','loading_gradient_threshold_m_inv',
              'locked_frequency_THz','gain_length_nm','linewidth_FWHM_GHz',
              'after_pump_amplitude_decay_ps']
    out = []
    for key,(lo,hi,kind) in spans.items():
        per_material = []
        for m in MATS:
            rows = [gain_branch(m,{key:float(v)}) for v in np.geomspace(lo,hi,21)]
            per_material.append(dict(material=m['name'],
                metrics={k:spread(rows,k) for k in fields},
                small_occupation_points=sum(r.get('phonon_number',0)<10 for r in rows)))
        score = max(v.get('decades') or 0 for row in per_material
                    for v in row['metrics'].values())
        out.append(dict(input=key, range=[lo,hi], classification=kind,
            range_origin='exploratory, not measured intervals',
            per_material=per_material, max_finite_decades=score))
    return sorted(out,key=lambda r:r['max_finite_decades'],reverse=True)

if __name__ == '__main__':
    scans = []
    baselines = []
    for m in MATS:
        baselines.append({'material':m['name'], 'g0_scales330':evaluate(m,31,BASE),
                          'TF_reference':evaluate(m,31,{'force':0.})})
        for key, (lo,hi,kind,origin) in ranges(m).items():
            grid = np.geomspace(lo,hi,41) if lo>0 else np.linspace(lo,hi,41)
            rows = [collect(m,key,float(v)) for v in grid]
            states = sorted(set(r['state'] for r in rows))
            metrics = {k:spread(rows,k) for k in FIELDS}
            finite = [v['decades'] for v in metrics.values() if v['decades'] is not None]
            scans.append(dict(material=m['name'], input=key, range=[lo,hi],
                classification=kind, range_origin=origin, states=states,
                phase_crossing=any(r['r']<0 for r in rows) and any(r['r']>=0 for r in rows),
                metrics=metrics, max_finite_decades=max(finite) if finite else None))
    ranking = []
    for key in ranges(MATS[0]):
        matches = [r for r in scans if r['input']==key]
        ranking.append(dict(input=key, classification=matches[0]['classification'],
            max_finite_decades=max(r['max_finite_decades'] or 0 for r in matches),
            phase_crossing=any(r['phase_crossing'] for r in matches)))
    ranking.sort(key=lambda r:r['max_finite_decades'], reverse=True)
    # Matched-isotope comparison: physical harmonic/quartic coefficients fixed,
    # photon/electron energies fixed, physical force zero, fractional damping fixed.
    isotope = []
    for m in MATS:
        ratio = (1875.612945e6/MH) if m['mass']==MH else MH/1875.612945e6
        a = evaluate(m,31,BASE)
        z = parameters(m, DEFAULT | BASE)
        # Hold the calibrated physical spring, rather than recalibrating E_INS
        # after changing the mass. Reuse the original canonical matrix builder.
        k = -2*a['r'] if a['r']<0 else a['r']
        scale = np.sqrt(ratio)
        _, L = canonical(k,z['g'],DEFAULT['photon_ratio']*scale,
            DEFAULT['electron_strength'],z['we']*scale,DEFAULT['gamma'])
        pole = min([v for v in np.linalg.eigvals(L) if v.imag>1e-8],key=lambda v:v.imag)
        energy = m['E']/scale
        b = dict(threshold_loading=a['threshold_loading'],
            static_displacement_A=a['static_displacement_A'],
            energy_per_isotope_eV=a['energy_per_isotope_eV'],
            correlation_length_nm=a['correlation_length_nm'],
            amplitude_mode_THz=float(pole.imag*energy/(2*np.pi*HBAR)/1e12),
            ringdown_ps=float(-HBAR/(energy*pole.real)/1e-12))
        isotope.append(dict(host=m['name'], mass_ratio=float(ratio), original=a, substituted=b))
    print(json.dumps(dict(selected='screened dipole plus quartic; force off for autonomous branch',
        five_number_meanings={'threshold':'xcrit,ncrit', 'frequency':'static lock0; amplitude mode separate',
          'domain_size':'not selected; correlation length ranked separately',
          'energy':'static minimum relative to disordered',
          'lifetime':'infinite deterministic persistence; ringdown ranked separately'},
        inactive_controls={'force':0., 'switch_mask':31},
        ranking=ranking, scans=scans, baselines=baselines, matched_isotopes=isotope,
        critical_screening=[critical_screen(m) for m in MATS],
        gain_reference=[gain_branch(m) for m in MATS],
        gain_threshold_cases=[gain_branch(m,{'pump_ratio':p}) for m in MATS for p in [.5,1.,2.,4.]],
        gain_ranking=gain_sweep(),
        domain_radius_result='unselected at every input; gradient alone does not fix a bulk radius',
        top_tool_input='finite-q longitudinal screening/force-constant kernel in the actual loaded phase'),
        indent=2, allow_nan=False))

