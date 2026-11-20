"""One-off g7_gap estimates. NumPy/SciPy only; not a fitted experimental model.
Run from checkout root: python g7_gap/motion_once.py
Writes selective results.json beside this script and exercises model guards once.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaln, logsumexp, roots_genlaguerre, roots_legendre

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from g4_beam.forward import thermal
from g4_beam.forward import model as DD
from g4_beam.merit.api import reference
from g7_micro.model import FrozenResponse

OUT = Path(__file__).resolve().parent
QE = 1.602176634e-19
KB_EV = 8.617333262145e-5
MD, MT = 2.01410, 3.0160493
B_DD = 31.3970  # sqrt(keV), retained Bosch-Hale convention
MU_DD = MD / 2
LABS = (.5, 1., 2., 5., 13.)


def gamow_b(m1, m2, zprod=1):
    return B_DD * zprod * math.sqrt((m1*m2/(m1+m2))/MU_DD)


def log_sigma(e, b, u=0.):
    """Constant S=1; Eq7 with frozen S. Energies keV, cross-section arbitrary units."""
    e = np.maximum(np.asarray(e), 1e-300)
    return -.5*np.log(e*(e+u)) - b/np.sqrt(e+u)


def log_shift_ratio(e, u, b):
    return .5*math.log(e/(e+u)) + b*(1/math.sqrt(e)-1/math.sqrt(e+u))


def apparent_u(logf, e, b):
    """Smallest nonnegative Eq7 shift, not a fit; null if above its monotone branch."""
    if logf <= 0:
        return 0. if abs(logf) < 1e-12 else None
    upper = b*b-e
    if upper <= 0 or logf > log_shift_ratio(e, upper, b):
        return None
    return 1000*brentq(lambda u: log_shift_ratio(e, u, b)-logf,
                       0., upper, xtol=1e-13)


class Motion:
    def __init__(self, n=48):
        x, w = roots_genlaguerre(n, .5)
        c, a = roots_legendre(n)
        self.x = x[:, None]
        self.c = c[None, :]
        self.logw = np.log(w[:, None]) - gammaln(1.5) + np.log(a[None, :]/2)

    def average(self, lab, theta_eV, m1=MD, m2=MD, zprod=1, u=0., flux=True):
        alpha, beta = m2/(m1+m2), m1/(m1+m2)
        e0 = alpha*lab
        if theta_eV == 0:
            return float(log_sigma(e0, gamow_b(m1,m2,zprod),u))
        k = self.x * theta_eV/1000  # Maxwell target kinetic energy
        er = e0 + beta*k - 2*np.sqrt(e0*beta*k)*self.c
        logh = log_sigma(er, gamow_b(m1,m2,zprod),u)
        if flux:
            logh = logh + .5*np.log(er/e0)  # g/v_projectile
        return float(logsumexp(logh+self.logw))

    def fixed_k(self, lab, k_eV, flux=True):
        e0 = lab/2
        er = e0 + k_eV/2000 - 2*np.sqrt(e0*k_eV/2000)*self.c.ravel()
        _, w = roots_legendre(self.c.size)
        h = log_sigma(er, B_DD)
        if flux:
            h += .5*np.log(er/e0)
        return float(logsumexp(h+np.log(w/2)))

    def source_average(self, labs, theta_eV, u=.040, channel='p'):
        """Actual common proton/BH neutron kernel, not frozen S; chunked in energy."""
        labs=np.atleast_1d(labs)
        if theta_eV==0:
            return DD.reaction_sigma(labs/2,u,channel,'proton')
        result=np.empty(len(labs))
        k=self.x*theta_eV/1000
        for start in range(0,len(labs),64):
            e=labs[start:start+64,None,None]/2
            er=e+.5*k-2*np.sqrt(e*.5*k)*self.c
            h=DD.reaction_sigma(er,u,channel,'proton')*np.sqrt(er/e)
            result[start:start+64]=np.sum(h*np.exp(self.logw),axis=(1,2))
        return result


def source_yield(response, rule, theta, channel='p'):
    value=0.
    for lab,share,en,weights in response.tracks:
        value+=share*(float(weights@rule.source_average(en,theta,channel=channel))+
                      response.nads*1e-31*rule.source_average([lab],theta,channel=channel)[0])
    return value/QE


def motion_row(rule, label, lab, theta, m2=MD):
    b = gamow_b(MD,m2)
    e0 = m2/(MD+m2)*lab
    base = float(log_sigma(e0,b))
    ls = rule.average(lab,theta,m2=m2,flux=False)-base
    lr = rule.average(lab,theta,m2=m2,flux=True)-base
    return dict(state=label, target_mass_amu=m2, lab_keV=lab, cm_stationary_keV=e0,
                theta_component_eV=theta, logF_sigma=ls, logF_rate=lr,
                Uapp_sigma_eV=apparent_u(ls,e0,b), Uapp_rate_eV=apparent_u(lr,e0,b))


def thick(rule, lab, theta, m2=MD, n=48, cutoff=.1):
    """Homogeneous sqrt(E_lab) stopping, no sheet, pure D+, finite .1-keV lab cutoff."""
    x,w = roots_legendre(n)
    energies = cutoff+(lab-cutoff)*(x+1)/2
    logs = np.array([rule.average(e,theta,m2=m2) for e in energies])
    weights = w*(lab-cutoff)/2/np.sqrt(energies)
    return float(logsumexp(logs+np.log(weights)))


def thick_shift(lab, u, n=64, cutoff=.1):
    x,w = roots_legendre(n)
    energies = cutoff+(lab-cutoff)*(x+1)/2
    return float(logsumexp(log_sigma(energies/2,B_DD,u)+
                          np.log(w*(lab-cutoff)/2/np.sqrt(energies))))


def thick_u(lab, logy):
    if logy <= thick_shift(lab,0):
        return 0.
    return 1000*brentq(lambda u: thick_shift(lab,u)-logy,0.,100.,xtol=1e-12)


def zp_theta(hw_eV, temp_K):
    return hw_eV/(2*math.tanh(hw_eV/(2*KB_EV*temp_K)))


def source_hashes():
    paths = ['g7_micro/README.md','g6_plateau/discriminators.csv',
             'g6_collision/README.md','g4_beam/priordata/experiments.csv',
             'g4_beam/forward/thermal.py','g4_beam/forward/bulk_screening.csv',
             'g4_beam/czerski/normalization.csv','g7_gap/motion_once.py',
             'g3_lattice/lattice_inputs.json',
             'g7_gap/candidate_ledger.csv','g7_gap/experiment_rows_addendum.csv',
             'g4_beam/merit/profiles.json','g4_beam/forward/model.py',
             'g4_beam/forward/consumer.py']
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    packet=ROOT/'g7_gap/source_observables.md'
    if packet.is_file():
        hashes['g7_gap/source_observables.md']=hashlib.sha256(packet.read_bytes()).hexdigest()
    return hashes


def main():
    start = time.monotonic()
    rule = Motion()
    pd = json.loads((ROOT/'g3_lattice/lattice_inputs.json').read_text())['materials']['PdD']
    hw_pd = 1973.269804*math.sqrt(pd['Kii_eV_A2']/pd['mass_eVc2'])
    rows = []
    for temp in (263.15,300.,323.15,373.15,423.15,473.15,1080.):
        for lab in LABS:
            rows.append(motion_row(rule,f'classical_{temp:g}K',lab,KB_EV*temp))
    for hw in (.030,hw_pd,.060,.100):
        for mass,name in ((MD,'D'),(MT,'T')):
            omega = hw*math.sqrt(MD/mass)
            for lab in LABS:
                rows.append(motion_row(rule,f'{name}_harmonic_hwD{hw:g}eV_300K',
                                       lab,zp_theta(omega,300),m2=mass))
    for theta in (.1,1.,10.,100.):
        for lab in LABS:
            rows.append(motion_row(rule,f'conditional_hot_theta{theta:g}eV',lab,theta))
    for lab in LABS:
        rows.append(motion_row(rule,'T_classical_300K',lab,KB_EV*300,m2=MT))

    thick_rows=[]
    for theta in (KB_EV*300,zp_theta(hw_pd,300),1.,10.,100.):
        ref = thick(rule,13.,theta)-thick(rule,13.,0.)
        for lab in LABS:
            logy = thick(rule,lab,theta)
            logf = logy-thick(rule,lab,0.)
            thick_rows.append(dict(lab_keV=lab,theta_component_eV=theta,
                                  logF_rate_thick=logf,Uapp_thick_eV=thick_u(lab,logy),
                                  logF_normalized_to_13lab=logf-ref,
                                  stopping='sqrt(E_lab); homogeneous; pure D+; cutoff .1keV lab'))

    # Source-specific normalization: actual PdD frozen profiles, sheet, 85/15
    # electrical fractions, 5% spread and original .1-keV lab stopping cutoff.
    responses={lab:FrozenResponse(lab/2,reference('PdD').profiles(lab/2)) for lab in LABS}
    source_rows=[]
    source_identity=0.
    for theta in (KB_EV*300,zp_theta(hw_pd,300)):
        cold_ref=responses[13.].yield_per_C('two_site',0.,.040,channel='n',exact=True)
        hot_ref=source_yield(responses[13.],rule,theta,channel='n')
        for lab,response in responses.items():
            channel='n' if lab>=9. else 'p'
            cold=response.yield_per_C('two_site',0.,.040,channel=channel,exact=True)
            hot=source_yield(response,rule,theta,channel=channel)
            source_identity=max(source_identity,abs(source_yield(response,rule,0.,channel)/cold-1))
            lognorm=math.log(hot/cold)-math.log(hot_ref/cold_ref)
            def normalized(u):
                y=response.yield_per_C('two_site',1.,u,channel=channel,exact=True)
                ref=responses[13.].yield_per_C('two_site',1.,u,channel='n',exact=True)
                return math.log(y/cold)-math.log(ref/cold_ref)
            delta=None if lab==13. else 1000*(brentq(lambda u:normalized(u)-lognorm,
                                                      .040,.060,xtol=1e-13)-.040)
            source_rows.append(dict(material='PdD',lab_keV=lab,channel=channel,
                                    theta_eV=theta,baseline_U_eV=40.,
                                    stationary_reactions_per_true_C=cold,
                                    moving_reactions_per_true_C=hot,
                                    logF_absolute=math.log(hot/cold),
                                    logF_own_reference_shape=lognorm,
                                    delta_Uapp_normalized_eV=delta,
                                    meaning='conditional motion increment in actual frozen source consumer; no refit; fixed stitch cancels enhancement'))

    mixtures=[]
    e0=.25
    for theta in (.1,1.,10.,100.):
        logf=rule.average(.5,theta)-float(log_sigma(e0,B_DD))
        for target_u in (300.,1000.,3000.):
            logtarget=log_shift_ratio(e0,target_u/1000,B_DD)
            fraction=math.expm1(logtarget)/math.expm1(logf)
            mixtures.append(dict(lab_keV=.5,theta_hot_eV=theta,target_Uapp_eV=target_u,
                                 logF_all_hot=logf,required_hot_number_fraction=fraction,
                                 feasible_fraction=fraction<=1.,
                                 meaning='cold theta=0 plus independent Maxwell hot target; not measured'))
    fixed=[]
    for lab in (.5,2.):
        for k in (1.,10.,100.,500.):
            base=float(log_sigma(lab/2,B_DD))
            lf=rule.fixed_k(lab,k)-base
            fixed.append(dict(lab_keV=lab,target_K_eV=k,logF_rate_isotropic=lf,
                              Uapp_eV=apparent_u(lf,lab/2,B_DD)))

    charges=[]
    for name,m1,m2,zp in [('DD',MD,MD,1),('D_Li6',MD,6.01512,3),
                          ('Li7_p',7.016004,1.007825,3),
                          ('F19_p',18.998403,1.007825,9),
                          ('F19_D',18.998403,MD,9)]:
        b=gamow_b(m1,m2,zp)
        for fixed_in in ('Ecm5keV','vrel_DD_at5keV'):
            mu=m1*m2/(m1+m2)
            ec=5. if fixed_in=='Ecm5keV' else 5.*mu/MU_DD
            lab=ec*(m1+m2)/m2
            lf=rule.average(lab,KB_EV*300,m1,m2,zp)-float(log_sigma(ec,b))
            charges.append(dict(reaction=name,comparison=fixed_in,zprod=zp,mu_amu=mu,
                                cm_keV=ec,logF_rate_300K=lf,Uapp_eV=apparent_u(lf,ec,b)))
    source_charge=[]
    for name,m1,m2,zp,window in [
        ('DD',MD,MD,1,(2.5,5.,13.)),
        ('D_Li6',MD,6.01512,3,(30.,50.,75.)),
        ('Li7_p',7.016004,1.007825,3,(437.,1000.,2071.)),
        ('F19_p',18.998403,1.007825,9,(6410.,6466.,6819.)),
        ('F19_D',18.998403,MD,9,(3089.,5250.,9200.))]:
        b=gamow_b(m1,m2,zp)
        for lab in window:
            ec=m2/(m1+m2)*lab
            lf=rule.average(lab,KB_EV*300,m1,m2,zp)-float(log_sigma(ec,b))
            source_charge.append(dict(projectile_target=name,lab_keV=lab,cm_keV=ec,
                                      logF_rate_300K=lf,Uapp_eV=apparent_u(lf,ec,b),
                                      meaning='smooth constant-S thin proxy at source kinematics, not a resonance-strength or thick-target refit'))
    empirical=[]
    for name,u,err,z in [('DD_Pd_Kasagi2002',310.,20.,1),
                         ('D_Li_Pd_Kasagi2004',1500.,310.,3),
                         ('Li7_p_hardPd_Cvetinovic2023',2860.,190.,3),
                         ('F19_p_hardPd_Cvetinovic2023',18700.,1500.,9),
                         ('F19_D_hardPd_Cvetinovic2023',18200.,3300.,9),
                         ('F19_D_softPd_Cvetinovic2023',3200.,1900.,9)]:
        empirical.append(dict(source_key=name,U_reported_eV=u,error_eV=err,zprod=z,
                              U_over_Z_eV=u/z,U_over_Z2_eV=u/z**2,
                              error_meaning=('statistical only; additional systematic 50eV'
                                             if name=='DD_Pd_Kasagi2002' else 'source reported'),
                              warning='not matched sites/kinematics/reference; ratios do not establish a charge law'))

    overlap=[]
    for label,flux in [('author_flux',1.1e19),('140uA_assumed_3mm_diameter',140e-6/QE/(math.pi*(.0015)**2))]:
        for r_nm in (1.,10.,100.):
            for tau_s in (1e-14,1e-12,1e-10,1e-9):
                lam=flux*math.pi*(r_nm*1e-9)**2*tau_s
                overlap.append(dict(flux_assumption=label,flux_m2_s=flux,radius_nm=r_nm,
                                    lifetime_s=tau_s,mean_extra_arrivals=lam,
                                    probability_at_least_one_extra=-math.expm1(-lam),
                                    lifetime_for_mean_one_s=1/(flux*math.pi*(r_nm*1e-9)**2)))

    track=[]
    with (ROOT/'g4_beam/czerski/normalization.csv').open() as f:
        for row in csv.DictReader(f):
            mat=row['material'];nd=float(row['nD_m3'])
            for temp in (800.,1080.,1500.):
                for er in (.1,1.,2.):
                    u={'Pd':.400,'Ti':.370,'Zr':.340}[mat]
                    p=float(thermal.probability(nd,Ue_hot_keV=u,E_A_eV=thermal.BARRIERS[mat],
                            exposure_m3_s=1.87e-34,T_K=temp,E_R_eV=er))
                    track.append(dict(material=mat,nD_m3=nd,T_K=temp,ER_eV=er,Uhot_eV=u*1000,
                                      A_m3_s=1.87e-34,branch_assumed=1.,reactions_per_incident_D=p))
    nh=6e28;volume=1870e-27;deltaT=1080.-300.
    heat=3*nh*volume*KB_EV*deltaT
    energy=dict(host_number_density_m3=nh,volume_nm3=1870.,initial_K=300.,final_K=1080.,
                classical_host_heat_eV=heat,
                volume_nm3_heatable_by_2keV=2000/(3*nh*KB_EV*deltaT)*1e27,
                meaning='uniform classical host only; illustrative, no electrons or fuel heat; not a measured track')

    recoil=[]
    for ec in (.25,.5,1.):
        max_ec=2*ec
        ls=float(log_sigma(max_ec,B_DD,.04)-log_sigma(ec,B_DD,.04))
        deficit=float(log_sigma(max_ec,B_DD,.04)-log_sigma(ec,B_DD,3.))
        recoil.append(dict(nominal_DD_cm_keV=ec,cold_single_ion_energy_keV=2*ec,
                           maximum_two_D_relative_energy_keV=max_ec,
                           local_sigma_ratio_max_vs_40eV=math.exp(ls),
                           local_sigma_max_over_3keV_shift=math.exp(deficit),
                           required_additional_encounter_prefactor=math.exp(-deficit),
                           limitation='local S frozen, no flux/path multiplicity bound; same binary pair stays at original Erel'))

    # One settled numerical/model guard: protects semantics and resolves adequate quadrature.
    finer=Motion(72)
    checkcases=[(.5,KB_EV*300),(.5,zp_theta(.100,300)),(.5,100.),(13.,100.)]
    diff=max(abs(rule.average(e,t)-finer.average(e,t)) for e,t in checkcases)
    thickdiff=abs(thick(rule,.5,100.,n=48)-thick(finer,.5,100.,n=72))
    identities=[]
    for lab in LABS:
        identities.append(abs(rule.average(lab,0.)-float(log_sigma(lab/2,B_DD))))
    # equal-mass collision: v'=(cos^2 a, sin a cos a), u=(sin^2 a,-sin a cos a).
    binary=max(abs(np.dot(np.array([math.cos(a)**2,math.sin(a)*math.cos(a)])-
                                 np.array([math.sin(a)**2,-math.sin(a)*math.cos(a)]),
                                 np.array([math.cos(a)**2,math.sin(a)*math.cos(a)])-
                                 np.array([math.sin(a)**2,-math.sin(a)*math.cos(a)]))-1)
               for a in np.linspace(0,math.pi/2,17))
    # Compare exact low-T average to analytic expansion for rate weighting h=E^-1/2 exp(-b/sqrtE).
    ec=1.;theta=1e-5;beta=.5
    expected=beta*(theta/1000)*(-B_DD/(2*ec**1.5)+B_DD**2/(4*ec**2))
    exact=rule.average(2.,theta)-float(log_sigma(ec,B_DD))
    expansion_error=abs(exact/expected-1)
    with (OUT/'experiment_rows_addendum.csv').open() as f:
        experiments=list(csv.DictReader(f))
    with (OUT/'candidate_ledger.csv').open() as f:
        ledger=list(csv.DictReader(f))
    assert all(None not in row for row in experiments+ledger)
    keys={row['row_id'] for row in experiments}
    assert len(keys)==len(experiments)
    model_keys={'collision_2keV','dynamic_lowq'}
    assert all(k in keys|model_keys for row in ledger for k in row['measured_keys'].split(';'))
    assert len(ledger)==15 and len({row['candidate_id'] for row in ledger})==15
    guards=dict(cold_identity_max_abs=max(identities),equal_mass_elastic_relative_energy_error=binary,
                quadrature_48_vs_72_max_abs_log=diff,thick_48_vs_72_abs_log=thickdiff,
                lowT_expansion_relative_error=expansion_error,
                actual_source_cold_identity_max_relative=source_identity,
                candidate_rows=len(ledger),experiment_rows=len(experiments),
                population_not_assigned=True,thick_cutoff_lab_keV=.1)
    assert max(identities)<1e-13 and binary<1e-14
    assert diff<2e-6 and thickdiff<2e-5 and expansion_error<2e-5
    assert source_identity<1e-12
    result=dict(contract='private g7_gap one-off conditional estimates; no experimental likelihood or measured U',
                model='constant astrophysical S; exact nonrelativistic relative kinematics; Eq7 shift inversion',
                motion=rows,illustrative_thick=thick_rows,conditional_hot_mixtures=mixtures,
                fixed_target_energy=fixed,charge_motion=charges,reported_charge_comparisons=empirical,
                source_kinematics_motion=source_charge,
                actual_source_profile_motion=source_rows,
                overlap=overlap,source_thermal_spike=track,track_energy_budget=energy,
                oscillator_hw_PdD_eV=hw_pd,oscillator_source_status=pd['status'],
                cold_cascade_bounds=recoil,guards=guards,input_sha256=source_hashes(),
                runtime_seconds=time.monotonic()-start,numpy=np.__version__)
    (OUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'guards':guards,'runtime_seconds':result['runtime_seconds'],
                      'motion_rows':len(rows),'output':'g7_gap/results.json'},indent=2))


if __name__=='__main__':
    main()
