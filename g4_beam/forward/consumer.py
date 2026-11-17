"""Physical yield/C and R(E) from energy-indexed transport profiles.

Profiles: callable(E_cm_keV) -> (layers, nads_m2, state), or exact-key dict.
Layer = (z_top_m,z_bottom_m,n_host_m3,n_D_m3,Z,M_amu,n_metal_m3).
Straight-ahead CSDA, source profiles frozen for each nominal beam setting.
No transport evolution is performed here. No detector factor enters physical yield/C.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from functools import lru_cache
from scipy.interpolate import PchipInterpolator
from . import model as M
from . import thermal

QE = 1.602176634e-19
COULOMB_EV_M = 1.4399645478e-9
WKB_NUCLEAR_RADIUS_M = 4.e-15       # fixed model inner boundary, NOT the 7 fm resonance radius


def _potential_cross_section(potential, max_E_cm, reaction):
    """s-wave WKB action enhancement of bare Bosch-Hale, same nuclear radius.

    potential=(r_m,V_eV): TOTAL DD repulsive potential, not induced potential.
    Tables must cover the fixed 4 fm WKB nuclear radius and the outer turning
    point at 50 eV CM. Samples below 4 fm do not move that physical boundary.
    No inward extrapolation is permitted. Pre-exponential stays bare; unlike
    Eq7 this is a declared WKB potential model, not an inferred constant Ue.
    """
    r,v = map(lambda a:np.asarray(a,float),potential)
    if r.ndim != 1 or r.shape != v.shape or len(r)<4 or np.any(~np.isfinite(r)) or np.any(np.diff(r)<=0) or r[0]<=0 or np.any(~np.isfinite(v)):
        raise ValueError('Potential requires increasing finite positive r_m and finite V_eV')
    if r[0]>WKB_NUCLEAR_RADIUS_M:
        raise ValueError('Potential grid must support the fixed 4 fm WKB nuclear radius; short-range extrapolation is forbidden')
    if v[-1]>50 or COULOMB_EV_M/r[-1]>50:
        raise ValueError('Potential grid must include the outer turning point at 50 eV')
    interp=PchipInterpolator(np.log(r),v,extrapolate=False)
    interior=r[(r>WKB_NUCLEAR_RADIUS_M)&(r<r[-1])]
    points=np.log(interior[::max(1,len(interior)//30)])
    energy=np.geomspace(50,max(100,max_E_cm*1000),220)
    hbar=1.054571817e-34
    prefactor=2*np.sqrt(2*thermal.MU_KG*QE)/hbar
    # Integrating in log radius keeps the femtometre-to-angstrom scales resolved.
    def action(e, screened):
        f=lambda t:np.sqrt(max((float(interp(t)) if screened else COULOMB_EV_M/np.exp(t))-e,0))*np.exp(t)
        return prefactor*quad(f,np.log(WKB_NUCLEAR_RADIUS_M),np.log(r[-1]),epsabs=1e-19,epsrel=2e-5,limit=250,points=points)[0]
    delta=np.array([action(e,False)-action(e,True) for e in energy])
    def xs(e,ch):
        return M.reaction_sigma(e,0,ch,reaction)*np.exp(np.interp(np.log(np.maximum(e, .05)*1000),np.log(energy),delta))
    return xs


def stopping_eV_per_m(E_lab_keV, layer):
    """ICRU Pd/Ti/H electronic plus ZBL nuclear, densities in SI.

    For other hosts the caller must supply a stopping function; no guessed
    coefficient or silent generic element fallback is used.
    """
    _,_,nh,nd,z,m,_=layer
    a1={1:1.254,22:4.858,46:5.238}.get(z)
    if a1 is None:
        raise ValueError('Supply stopping_fn for host Z='+str(z))
    e=max(float(E_lab_keV),1e-8)
    se=a1*np.sqrt(e*1.007276/M.M_D)
    host=(se+float(M.Sn_atom(e,z,m)))*1e-19
    deut=(float(M.Se_atom(e,'H'))+float(M.Sn_atom(e,1,M.M_D)))*1e-19
    return nh*host+nd*deut

_GX, _GW = np.polynomial.legendre.leggauss(12)


@lru_cache(maxsize=2048)
def _csda_segments(e0, layers, stopping_fn, cutoff):
    """Reusable stopping solution, independent of nuclear/screening parameters."""
    energy=float(e0); segments=[]
    for layer in layers:
        if energy<=cutoff:break
        zt,zb,_,nd,_,_,_=layer
        def rhs(z,state):
            return [-stopping_fn(max(float(state[0]),cutoff),layer)/1000]
        def event(z,state):return state[0]-cutoff
        event.terminal=True;event.direction=-1
        sol=solve_ivp(rhs,(zt,zb),[energy],rtol=2e-7,atol=1e-8,
                      events=event,max_step=10e-9,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        end=float(sol.t[-1])
        segments.append((zt,end,nd,sol.sol))
        energy=float(sol.y[0,-1])
        if sol.t_events[0].size:break
    return segments


@lru_cache(maxsize=2048)
def _csda_quadrature(e0,frozen,stopping_fn,cutoff):
    parts=_csda_segments(e0,frozen,stopping_fn,cutoff)
    energies=[];depths=[];weights=[];intervals=[]
    for za,zb,nd,solution in parts:
        cuts=np.linspace(za,zb,max(1,int(np.ceil((zb-za)/10e-9)))+1)
        mid=(cuts[:-1]+cuts[1:])/2;half=np.diff(cuts)/2
        z=(mid[:,None]+half[:,None]*_GX).ravel()
        energies.append(solution(z)[0]);depths.append(z)
        weights.append((half[:,None]*_GW*nd*1e-31).ravel())
        intervals.extend((a,b,nd,solution) for a,b in zip(cuts[:-1],cuts[1:]))
    if not parts:return np.array([]),np.array([]),np.array([]),0.,[]
    return (np.concatenate(energies),np.concatenate(depths),
            np.concatenate(weights),parts[-1][1],intervals)


def csda_nodes(e0,layers,stopping_fn=stopping_eV_per_m,cutoff=.1,break_m=None):
    """Reusable authoritative stopping ODE and 12-point Gauss reaction nodes.

    Return lab keV, depth m, n_D*dz*1e-31 weights, and terminal depth.
    Split a screening step exactly, replacing only its intersected interval.
    """
    frozen=tuple(tuple(float(v) for v in row) for row in layers)
    en,z,w,length,intervals=_csda_quadrature(float(e0),frozen,stopping_fn,float(cutoff))
    if break_m is not None:
        for za,zb,nd,solution in intervals:
            if za<break_m<zb:
                cuts=np.array([za,break_m,zb])
                mid=(cuts[:-1]+cuts[1:])/2;half=np.diff(cuts)/2
                zz=(mid[:,None]+half[:,None]*_GX).ravel()
                keep=(z<za)|(z>zb)
                return (np.r_[en[keep],solution(zz)[0]],np.r_[z[keep],zz],
                        np.r_[w[keep],(half[:,None]*_GW*nd*1e-31).ravel()],length)
    return en,z,w,length


def forward_yields(E_cm_keV, profiles, *, Ue_keV=None, potential=None,
                   candidate='none', parameters=None, measurement=None,
                   reference_E_cm_keV=6.5, stopping_fn=stopping_eV_per_m,
                   E_cut_lab_keV=.1, reaction='proton'):
    """Return arrays R, physical_yield_per_C, counts_per_C, reported_signal.

    Energies are nominal per-deuteron CM keV. Exactly one of Ue_keV (scalar
    or callable(E_cm,z_m)) and potential=(r_m,total_V_eV) is required. Potential
    tables must cover the fixed 4 fm WKB inner boundary and outer turning point;
    a reporting-only radial window is rejected, never extrapolated inward.
    Candidate: none, layer, sigmoid_tt, fast_neutrals, thermal.
    Params: Ue_s/d_nm or E0; neutral f_n/E_n_lab_keV; thermal requires
    explicit Ue_hot, A (effective V*tau*branch m3 s), E_A_eV, T_track_K,
    E_R_eV. A is fitted at a declared ER/T convention; do not also rescale A
    when varying explicit ER/T. With V(r),
    layer/sigmoid use an explicitly supplied delta_Ue_surface_keV instead.
    Measurement: efficiency (scalar/callable), collected_fraction (ditto),
    current_A (ditto), extra_current_power, channel ('p' or 'n'), D1_fraction,
    D2_fraction, spread_FWHM. Fractions are electrical-current fractions;
    D2+ contributes TWO incident deuterons. Neutral delta=1 uses explicit
    neutral_reference_current_A; otherwise f_n is per incident charge.
    Counts/C is per collected cup coulomb; physical yield/C per TRUE charged
    coulomb, including uncounted neutrals. Additional current power affects
    only reported_signal and R, never either physical or counts/C output.
    """
    if (Ue_keV is None)==(potential is None):
        raise ValueError('Provide exactly one of Ue_keV and potential')
    if candidate not in ('none','layer','sigmoid_tt','fast_neutrals','thermal'):
        raise ValueError('Unconfirmed candidate '+str(candidate))
    p=dict(parameters or {});meas=dict(measurement or {})
    if candidate=='thermal':
        required={'Ue_hot','A','E_A_eV','T_track_K','E_R_eV'}
        if not required.issubset(p):raise ValueError('Explicit thermal parameters required: '+str(sorted(required)))
    energies=np.atleast_1d(np.asarray(E_cm_keV,float))
    if np.any(~np.isfinite(energies)) or np.any(energies<=0) or reference_E_cm_keV<=0:
        raise ValueError('Positive finite energies required')
    ch=meas.get('channel','p')
    if ch not in ('p','n'):raise ValueError('Explicit p or n channel required')
    eval_setting=lambda val,e:float(val(e) if callable(val) else val)
    f1,f2=meas.get('D1_fraction',.85),meas.get('D2_fraction',.15)
    if min(f1,f2)<0 or abs(f1+f2-1)>1e-10:raise ValueError('Current fractions must sum to one')
    gx,gw=np.polynomial.hermite_e.hermegauss(5);gw/=gw.sum()
    spread=meas.get('spread_FWHM',.05)/2.3548
    maxe=max(np.max(energies),reference_E_cm_keV,p.get('E_n_lab_keV',5)/2)*(1+spread*np.max(gx))+max(0.,float(p.get('delta_Ue_surface_keV',0.)))
    potential_xs=_potential_cross_section(potential,maxe,reaction) if potential is not None else None
    if potential is not None and candidate in ('layer','sigmoid_tt'):
        if 'delta_Ue_surface_keV' not in p:raise ValueError('Specify potential surface shift explicitly')
    def at_energy(e):
        layers,nads,_=profiles(e) if callable(profiles) else profiles[float(e)]
        layers=np.asarray(layers,float)
        if layers.ndim!=2 or layers.shape[1]!=7 or not len(layers) or np.any(~np.isfinite(layers)) or nads<0 or not np.isfinite(nads):
            raise ValueError('Seven-field SI layers and finite nonnegative adsorbate required')
        if np.any(layers[:,1]<=layers[:,0]) or np.any(layers[:,2:4]<0) or abs(layers[0,0])>1e-15 or np.any(np.abs(layers[1:,0]-layers[:-1,1])>1e-15):
            raise ValueError('Contiguous nonnegative-density layers must start at zero')
        def ion(e0):
            def xs(energy,z):
                ec=np.maximum(np.asarray(energy)/2,.05)
                if potential_xs is not None:base=0.
                elif callable(Ue_keV):
                    base=np.array([float(Ue_keV(e,zz)) for zz in np.atleast_1d(z)])
                else:base=float(Ue_keV)
                weight=0.
                if candidate=='layer':weight=np.asarray(z)<p['d_nm']*1e-9
                elif candidate=='sigmoid_tt':weight=1/(1+np.exp((e0/2-p['E0'])/.5))
                u=base+weight*((p.get('Ue_s',base)-base) if potential_xs is None else p.get('delta_Ue_surface_keV',0.))
                extra=p.get('additional_Ue_keV',0.)
                u+=np.array([float(extra(e,zz)) for zz in np.atleast_1d(z)]) if callable(extra) else float(extra)
                if np.any(u<0):raise ValueError('Screening energy must be nonnegative')
                if potential_xs is None:return M.reaction_sigma(ec,u,ch,reaction)
                return potential_xs(ec+u,ch)*M.reaction_S(ec,ch,reaction)/M.reaction_S(ec+u,ch,reaction)*(ec+u)/ec
            y=float(nads)*float(np.asarray(xs(e0,0)).item())*1e-31
            en,z,w,length=csda_nodes(e0,layers,stopping_fn,E_cut_lab_keV,
                                    p['d_nm']*1e-9 if candidate=='layer' else None)
            if len(en):y+=float(np.sum(w*xs(en,z)))
            if candidate=='thermal' and ch=='p':
                from .profile_thermal import density_squared_average
                # Volume-average along the frozen target's traversed track.
                nd2=density_squared_average(layers,length)
                y+=float(thermal.probability(np.sqrt(nd2),Ue_hot_keV=p['Ue_hot'],E_A_eV=p['E_A_eV'],exposure_m3_s=p['A'],T_K=p['T_track_K'],E_R_eV=p['E_R_eV']))
            return y
        yi=sum(w*(f1*ion(2*e*(1+spread*g))+2*f2*ion(e*(1+spread*g))) for g,w in zip(gx,gw))
        current=eval_setting(meas.get('current_A',1.),e)
        frac=eval_setting(meas.get('collected_fraction',1.),e)
        efficiency=eval_setting(meas.get('efficiency',1.),e)
        if current<=0 or frac<=0 or efficiency<0:raise ValueError('Invalid measurement current/fraction/efficiency')
        yn=0.
        if candidate=='fast_neutrals':
            ratio=p['f_n']
            if p.get('delta',0)==1:ratio*=meas['neutral_reference_current_A']/current
            yn=ratio*ion(p['E_n_lab_keV'])
        physical=(yi+yn)/QE
        counts=physical*efficiency/frac
        signal=counts/(current*frac)**meas.get('extra_current_power',0)
        return physical,counts,signal,yi/QE,yn/QE
    values=np.array([at_energy(float(e)) for e in np.r_[energies,reference_E_cm_keV]])
    if values[-1,2]<=0:raise ValueError('Reference signal is zero; R is undefined')
    return dict(E_cm_keV=energies,R=values[:-1,2]/values[-1,2],
                physical_yield_per_C=values[:-1,0],counts_per_C=values[:-1,1],
                reported_signal=values[:-1,2],ion_yield_per_C=values[:-1,3],
                neutral_yield_per_C=values[:-1,4],reference_signal=values[-1,2],
                reference_E_cm_keV=reference_E_cm_keV,channel=ch)
