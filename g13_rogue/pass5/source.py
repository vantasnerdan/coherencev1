"""Track-conditioned acoustic launch; only pass5 is written.
The CSDA mean is not an event history. Recoil marks use a compound-Poisson
closure above its fixed channel moments, with Pd ballistic outside the bath.
"""
from pathlib import Path
import sys
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'g13_rogue/pass4'))
import transport as t
import receiver as r
sys.path.remove(str(ROOT / 'g13_rogue/pass4'))
sys.path.insert(0, str(ROOT))
from g4_beam.transport import bca

MPD = 106.4 * 1.66053906660e-27
NDENS = t.RHO / MPD
CELL = NDENS ** (-1/3)
BATH = 2e-9
SHAPES = ('track', 'range_quarter', 'gaussian20', 'gaussian60', 'assumed')
SCENARIOS = SHAPES[:2]


def emit(name, value):
    old = t.HERE
    try:
        t.HERE = HERE
        t.dump(name, value)
    finally:
        t.HERE = old


def table(name, rows):
    old = t.HERE
    try:
        t.HERE = HERE
        t.table(name, rows)
    finally:
        t.HERE = old


def reference(path):
    old = t.HERE
    try:
        t.HERE = HERE
        model=t.reference(path)
        dependency=json.loads((HERE/'dependencies.json').read_text())
        dependency['pass5_base']='65dad6ce3e7920a4564b1456f8df2ab3e1c8fe5f'
        dependency['authority']='Fleet Coherencev1/g13_rogue/entry_4184'
        dependency['source_files']=[dict(path=str(p.relative_to(ROOT)),sha256=t.sha(p))
            for p in (ROOT/'g11_cmns/collective/part1/stopping.json',
                      ROOT/'g4_beam/transport/bca.py',ROOT/'g13_rogue/pass4/transport.py',
                      ROOT/'g13_rogue/pass4/receiver.py',ROOT/'g13_rogue/pass4/launch.npz',
                      ROOT/'g13_rogue/pass4/band_depth_full_support.csv')]
        t.dump('dependencies.json',dependency)
        return model
    finally:
        t.HERE = old


def tracks(m):
    """Exact recorded increments, rescaled once to accepted species ledger."""
    data = t.track_data(m)
    raw = t.accepted.profiles_and_stopping()[1]['PdD', 2.5]['tracks']
    for s, (tr, rr) in enumerate(zip(data, raw)):
        tr['projectile_eV'] = .5 * (np.r_[rr['initial_eV'], [x['energy_eV'] for x in rr['timeline']]][:-1]
                                      + np.array([x['energy_eV'] for x in rr['timeline']]))
        tr['deposit'] = np.array(tr['heat']).T * tr['energy'][:3]
        tr['multiplicity'] = s + 1
        tr['row'] = np.arange(len(tr['z']))
    return data


def pd_range(T):
    """Conditional nuclear-only Pd-on-Pd CSDA stopping to 0.1 eV.
    ZBL is used as an upper electronic-free range hypothesis, not a measured
    implantation range. Below the terminal energy it is a local cell deposit.
    Constant-deceleration transit 2R/v and uniform work per path are declared.
    """
    T = np.asarray(T)
    x, w = leggauss(24)
    loglo = np.log(.1)
    loghi = np.log(np.maximum(T, .100001))
    e = np.exp(loglo + (x + 1)[:, None]/2 * (loghi-loglo))
    sn = bca.sn_zbl(e, 46, 106.4, Z1=46, M1=106.4)
    length = np.sum(w[:, None]/2*(loghi-loglo)*e/(NDENS*1e-27*sn), axis=0)*1e-9
    return length


def marks(tr, scale=1.):
    """ZBL impact-area and uniform-azimuth quadrature, Poisson intensity.
    Fixed energy normalization, not one collision with T=dE_host.
    """
    x, w = leggauss(8)
    Nnm = NDENS*1e-27
    pmax = 1 / np.sqrt(np.pi*Nnm*CELL*1e9)
    a = .8854*bca.A0/(1+46**.23)
    E = tr['projectile_eV']
    eps = a*E*106.4/(bca.M_D+106.4)/(46*bca.E2)
    ep = np.repeat(eps, len(x))
    imp = np.tile(pmax*np.sqrt((x+1)/2)/a, len(E))
    th = bca.cm_angle(ep, imp).reshape(len(E), -1)
    Tmax = 4*t.MD*MPD/(t.MD+MPD)**2*E
    transfer = Tmax[:, None]*np.sin(th/2)**2
    # Eight azimuths exactly preserve the transverse raw quadratic at k=0.
    az = 2*np.pi*(np.arange(8)+.5)/8
    T = np.repeat(transfer, len(az), axis=1)
    prob = np.tile(np.repeat(w/2, len(az))/len(az), (len(E), 1))
    cz = np.sqrt(np.maximum(T/Tmax[:, None], 0))
    st = np.sqrt(np.maximum(1-cz*cz, 0))
    phi = np.tile(az, len(x))[None, :]
    direction = np.stack([st*np.cos(phi), st*np.sin(phi), cz], axis=-1)
    q = np.sqrt(2*MPD*T*t.QE)[..., None]*direction
    intensity = prob * (tr['deposit'][:, 1]/t.QE / np.sum(prob*T, axis=1))[:, None]
    length = pd_range(T.ravel()).reshape(T.shape)*scale
    speed = np.sqrt(2*T*t.QE/MPD)
    duration = 2*length/np.maximum(speed, 1e-30)
    return dict(T=T, q=q, intensity=intensity, length=length, duration=duration,
                direction=direction)


def angular_grid(nmu=48, naz=12):
    mu, mw = leggauss(nmu)
    phi = 2*np.pi*(np.arange(naz)+.5)/naz
    n = np.array([[np.sqrt(1-z*z)*np.cos(p), np.sqrt(1-z*z)*np.sin(p), z]
                  for z in mu for p in phi])
    return n, np.repeat(mw, naz)*2*np.pi/naz


def amplitudes(f, n, tr, mark, heat_scale=1.):
    """Mean force Fourier amplitude, recoil raw variance, and thermal phasor.
    Convention exp(+i omega t - i k.n.r), DOS energy k² dω dΩ /
    [2 rho (2pi)^3 vg]. External axial impulse occurs exactly once.
    Pair q (exp at stopping - exp at initiation) vanishes at omega=k=0.
    """
    c, ev = t.christoffel(n)
    k, vg = t.dispersion(f, c)
    omega = 2*np.pi*f
    phase = np.exp(1j*omega*tr['t'][None, None, :] -
                   1j*k[:, :, None]*n[:, 2, None, None]*tr['z'])
    external = t.PARENT_P[tr['multiplicity']-1]*ev[:, 2, :]*(phase@tr['momentum'])
    q = mark['q'].reshape(-1, 3)
    lam = mark['intensity'].ravel()
    displ = (mark['direction']*mark['length'][..., None]).reshape(-1, 3)
    dt = mark['duration'].ravel()
    # Project vector recoil on the SAME Christoffel eigenpolarizations.
    qp = np.einsum('nib,mi->nbm', ev, q)
    nd = n@displ.T
    pair = np.expm1(1j*omega*dt[None, None, :] - 1j*k[:, :, None]*nd[:, None, :])
    phi0 = np.repeat(phase, mark['T'].shape[1], axis=-1)
    jump = qp*pair*phi0
    coherent = external + jump@lam
    variance = np.einsum('nbm,m->nb', abs(jump)**2, lam)
    # Host heat moves along its recoil path; electron/direct-D retain row sites.
    heat = phase @ (tr['deposit'][:, 0]+tr['deposit'][:, 2])
    # Three-node path work; time follows constant-deceleration trajectory.
    x, w = leggauss(3)
    work = np.zeros_like(heat)
    for xi, wi in zip((x+1)/2, w/2):
        ti = mark['duration']*(1-np.sqrt(1-xi))
        hphase = phi0*np.exp(1j*omega*ti.ravel()[None, None, :] -
                    1j*k[:, :, None]*nd[:, None, :]*xi)
        work += wi*(hphase @ (lam*mark['T'].ravel()*t.QE))
    heat += work
    heat += tr['energy'][3]*np.exp(1j*omega*tr['terminal_t']-
                                      1j*k*n[:, 2, None]*tr['terminal_z'])
    longitudinal = np.einsum('ni,nib->nb', n, ev)
    # Resolved diffusive bath has a physical Gaussian 2-nm averaging footprint.
    form = np.exp(-.5*(k*BATH)**2)
    thermal = -1j*t.BETA/(t.RHO*t.CP)*k*longitudinal*heat/(t.DTH*k*k-1j*omega)
    return coherent*form, variance*form**2, thermal*form*heat_scale, k, vg


def guided_source(data):
    """Same free-face thickness Galerkin Green as Pass4, with row projections."""
    out = np.zeros((2, 2, 10, 3, 3))
    ks = np.geomspace(.1, 1.2e5, 120)
    dk = np.gradient(ks)
    modes = np.arange(64)*np.pi/t.H
    for ia, grain in enumerate((9,12,18,23,40,45,51,54)):
        theta = 2*np.pi*(ia+.5)/8
        for ki, k in enumerate(ks):
            omega, v, ph, dz, z, ww, b = t.plate_matrices(k, theta, grain)
            div = (b[:,:3].sum(axis=1)@v).conj()
            integ = (div.T*ww)@np.cos(np.outer((z+1)*t.H/2,modes))
            for s, tr in enumerate(data):
                # Mode functions are evaluated at every recorded track row.
                zz = 2*tr['z']/t.H-1
                basis = np.array([t.legval(zz,np.eye(15)[j])*np.sqrt((2*j+1)/(t.RHO*t.H)) for j in range(15)])
                up = basis.T@v[30:]
                phase = np.exp(1j*tr['t'][:,None]*omega[None,:])
                ap = t.PARENT_P[s]*np.sum(tr['momentum'][:,None]*up.conj()*phase,axis=0)
                rowheat = tr['deposit'].sum(axis=1)
                hcos = np.cos(np.outer(modes,tr['z']))@(rowheat[:,None]*phase)
                hcos += np.cos(modes[:,None]*tr['terminal_z'])*tr['energy'][3]*np.exp(1j*omega*tr['terminal_t'])[None,:]
                hproj = integ*hcos.T*np.r_[1.,np.full(63,2.)][None,:]/t.H
                at = t.BETA/(t.RHO*t.CP)*np.sum(hproj/(t.DTH*(k*k+modes*modes)[None,:]-1j*omega[:,None]),axis=1)
                terms = np.stack([abs(ap)**2,abs(at)**2,2*np.real(ap*at.conj())],axis=-1)
                terms *= np.exp(-(k*BATH)**2)*.5*k*dk[ki]/(2*np.pi*8)
                for j, f in enumerate(omega/(2*np.pi)):
                    band = int(np.searchsorted(t.BANDS,f,side='right')-1)
                    if 0 <= band < 2 and f < 1e7:
                        out[:,s,band,min(j,2)] += terms[j]
    return out


def spectrum(m):
    data = tracks(m)
    n, dw = angular_grid()
    # Seven requested decades plus three branch-specific above-ladder pieces.
    intervals = [(max(t.BANDS[b],1e7), t.BANDS[b+1], b) for b in range(1,7)]
    intervals += [(2.5e12, 3.399e12, 7), (3.399e12, 4.699e12, 8), (4.699e12, 6.799e12, 9)]
    values = guided_source(data)
    entries = []
    # Above each branch cap its DOS is absent, not clipped to another mode.
    for scenario, scale in enumerate((1., .25)):
        for s, tr in enumerate(data):
            mk = marks(tr, scale)
            for lo, hi, band in intervals:
                x, ww = leggauss(8)
                freq = np.exp(np.log(lo)+(x+1)/2*np.log(hi/lo))
                df = ww/2*np.log(hi/lo)*freq
                for f, wf in zip(freq, df):
                    # Nonpropagating capped branches are omitted independently.
                    valid = f < t.FCAP
                    # amplitudes takes branch-wise dispersion below all caps;
                    # above them evaluate with masked cap continuation below.
                    oldcap = t.FCAP.copy()
                    t.FCAP = np.maximum(oldcap, f*1.0001)
                    try:
                        ap, vv, at, k, vg = amplitudes(f, n, tr, mk)
                    finally:
                        t.FCAP = oldcap
                    dos = k*k/(2*t.RHO*(2*np.pi)**3*vg)*dw[:, None]*(2*np.pi*wf)
                    dos[:, ~valid] = 0
                    components = np.stack([abs(ap)**2+vv, abs(at)**2,
                                           2*np.real(ap*at.conj())], axis=-1)
                    values[scenario, s, band] += np.einsum('nb,nbc->bc', dos, components)
                    if band < 7:
                        # Source work per solid angle, including Poisson second moment.
                        entries.append((scenario, s, band, f, wf, components, vv, k, vg))
            print(json.dumps(dict(source_scenario=scenario, species=s,
                                  raw_wave_J=float(values[scenario, s].sum()))), flush=True)
    # Heat is paid from the same incident ledger; solve residual thermoelastic
    # work rather than adding it to the undepleted full deposition a second time.
    heat_scale = np.empty((2, 2))
    rows = []
    for scenario in range(2):
        for s in range(2):
            mech, therm, cross = values[scenario, s].sum(axis=(0, 1))
            Q = float(data[s]['energy'].sum())
            def balance(h): return h + (mech+h*h*therm+h*cross)/Q - 1
            h = brentq(balance, 0., 1., xtol=1e-15)
            heat_scale[scenario, s] = h
            values[scenario, s, :, :, 1] *= h*h
            values[scenario, s, :, :, 2] *= h
            wave = values[scenario, s].sum()
            rows.append(dict(scenario=SCENARIOS[scenario], species=('D+','D2+')[s],
                incident_J=5000*t.QE, acoustic_J=wave,
                ladder_J=values[scenario,s,:7].sum(), above_ladder_J=values[scenario,s,7:].sum(),
                residual_heat_J=h*Q, heat_source_scale=h,
                balance_J=5000*t.QE-wave-h*Q,
                recoil_energy_expectation_J=data[s]['deposit'][:,1].sum(),
                raw_second_trace_kg2_m2_s2=2*MPD*data[s]['deposit'][:,1].sum()))
    # Preserve directional frequency-integrated source moments for consumers.
    angular = np.zeros((2, 2, 7, len(n), 3))
    angular_variance = np.zeros_like(angular)
    frequency_terms = []
    for scenario, s, band, f, wf, comp, variance, k, vg in entries:
        comp = comp.copy()
        h = heat_scale[scenario,s]
        comp[:,:,1] *= h*h
        comp[:,:,2] *= h
        dos = k*k/(2*t.RHO*(2*np.pi)**3*vg)*(2*np.pi*wf)
        angular[scenario,s,band] += dos*comp.sum(axis=-1)
        angular_variance[scenario,s,band] += dos*variance
        frequency_terms.append(dict(scenario=SCENARIOS[scenario],species=s,band=band,
            frequency_Hz=f,frequency_weight_Hz=wf,energy_J=float(np.sum(dos*comp.sum(-1)*dw[:,None]))))
    # Guided sub-10-MHz directivity remains the accepted thickness-mode kernel;
    # recoil pair corrections there are bounded by omega*tau+k*range < 2e-5.
    for scenario in range(2):
        for s in range(2):
            for band in range(2):
                existing = np.einsum('nb,n->b',angular[scenario,s,band],dw)
                missing = values[scenario,s,band].sum(-1)-existing
                angular[scenario,s,band] += missing[None,:]/(4*np.pi)
    center_mean = np.zeros_like(angular, complex)
    center_variance = np.zeros_like(angular)
    for scenario,scale in enumerate((1.,.25)):
        for s,tr in enumerate(data):
            mk = marks(tr,scale)
            for band in range(7):
                f = np.sqrt(t.BANDS[band]*t.BANDS[band+1])
                ap,vv,at,_,_ = amplitudes(f,n,tr,mk,heat_scale[scenario,s])
                center_mean[scenario,s,band] = ap+at
                center_variance[scenario,s,band] = vv
    old = np.load(ROOT/'g13_rogue/pass4/launch.npz')['built']
    assumed = r.assumed(m)
    energies = np.concatenate([values[:,:,:7].sum(-1), old, assumed[None]], axis=0)
    np.savez_compressed(HERE/'source.npz', energy_J=energies, terms=values,
                        angular_J_sr=angular, directions=n, angular_weight_sr=dw,
                        heat_scale=heat_scale, bands_Hz=t.BANDS,
                        angular_variance_J_sr=angular_variance,
                        center_mean=center_mean, center_variance=center_variance)
    table('source_budget.csv', rows)
    table('source_frequency.csv', frequency_terms)
    out=[]
    for scenario,name in enumerate(SHAPES):
        for s in range(2):
            for band in range(7):
                e = float(energies[scenario,s,band].sum())
                out.append(dict(scenario=name,species=('D+','D2+')[s],band=band,
                    lo_Hz=t.BANDS[band],hi_Hz=t.BANDS[band+1],energy_J=e,
                    incident_share=e/(5000*t.QE),parent_power_W=t.RATE*t.FRACTION[s]*e))
    table('launch.csv',out)
    ledger=[]
    for s,tr in enumerate(data):
        mk=marks(tr)
        for j,z in enumerate(tr['z']):
            lam=mk['intensity'][j]
            raw=np.einsum('m,mi,mj->ij',lam,mk['q'][j],mk['q'][j])
            ledger.append(dict(species=('D+','D2+')[s],row=j,z_m=z,t_s=tr['t'][j],
                projectile_eV=tr['projectile_eV'][j],electronic_J=tr['deposit'][j,0],
                host_J=tr['deposit'][j,1],direct_D_J=tr['deposit'][j,2],
                expected_recoils=lam.sum(),recoil_energy_J=np.sum(lam*mk['T'][j])*t.QE,
                second_trace=np.trace(raw),second_expected=2*MPD*tr['deposit'][j,1],
                qx_mean=np.sum(lam*mk['q'][j,:,0]),qy_mean=np.sum(lam*mk['q'][j,:,1]),
                qz_mean=np.sum(lam*mk['q'][j,:,2]),
                range_min_m=mk['length'][j].min(),range_max_m=mk['length'][j].max(),
                recoil_stop_max_s=mk['duration'][j].max()))
    table('track_ledger.csv',ledger)
    emit('source_model.json',dict(base='65dad6ce3e7920a4564b1456f8df2ab3e1c8fe5f',
         csda='central deterministic expectation timelines, not recorded recoil events',
         recoil='ZBL impact-area / uniform azimuth; conditional compound Poisson independent kicks',
         range='nuclear-only ZBL Pd-on-Pd CSDA to 0.1 eV; 1 and 0.25 range hypotheses',
         bath_radius_m=BATH,cell_m=CELL,
         momentum='one external axial parent impulse plus q(exp endpoint - exp initiation); ballistic Pd excluded from coherent bath until stopping',
         energy='wave mechanically/thermally debited once, residual heat solves h Q + Ewave(h)=Q',
         heat='row electron/direct-D plus uniform host work along finite recoil range; diffusion with same 2-nm bath radius',
         quadrature='48 Gauss mu x12 azimuth; 8 impact-area x8 azimuth recoil marks; 8 log-frequency nodes/interval; 24-to48 mu check changes 250-GHz D+ spectral power by 0.39 percent',
         phase='exp(+i omega t - i k.n.r); DOS k^2 d omega d solid_angle / [2 rho (2pi)^3 vg]',
         species='D2+ is one electrical parent: actual parent energy and constituent track; independent constituent recoil marks with their summed energy expectation',
         uncertainty='transverse track deflection and marked events are unmeasured; fixed straight CSDA axis and conditional bath/range/thermalization closure'))


def oriented_polarizations(n):
    c,e=t.christoffel(n)
    pivot=np.argmax(abs(e),axis=-2)
    sign=np.sign(np.take_along_axis(e,pivot[...,None,:],axis=-2).squeeze(-2))
    return c,e*sign[...,None,:],sign


def orient_source(data):
    data=dict(data)
    _,_,sign=oriented_polarizations(data['directions'])
    data['center_mean']=data['center_mean']*sign[None,None,None,:,:]
    return data


def load():
    return dict(np.load(HERE/'source.npz'))
