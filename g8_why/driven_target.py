"""Conditional steady target-motion projection and optimistic deposition budget."""
import csv
import hashlib
import json
import math
import numpy as np
from scipy.special import roots_genlaguerre, roots_legendre, gamma
from common import (ROOT, OUT, DD, QE, U, ENERGIES, MATERIALS, GX, GW,
                    response, yield_C, readout, observed, reference, write_csv)
from g4_beam.forward.consumer import csda_nodes, _csda_segments
from g4_beam.forward.profile_thermal import transport_stopping

KB = 8.617333262145e-5
AREA = 1e-4  # transported model's assumed 1 cm2 footprint, not measured spot
HW = .047607  # provisional saved PdD local oscillator, not measured dispersion


class Motion:
    def __init__(self, n=48):
        x, w = roots_genlaguerre(n, .5)
        c, a = roots_legendre(n)
        self.x, self.c = x[:, None], c[None, :]
        self.weight = w[:, None]/gamma(1.5)*a[None, :]/2

    def sigma_rate(self, lab, theta, ch):
        """<sigma(Erel)g>/v, smooth S(Erel), not a shift or sigma-only mean."""
        lab = np.atleast_1d(lab)
        if theta == 0:
            return DD.sigma(lab/2, U, ch)
        k = self.x*theta/1000
        out = np.empty(len(lab))
        for start in range(0, len(lab), 64):
            e = lab[start:start+64, None, None]/2
            er = e+.5*k-2*np.sqrt(e*.5*k)*self.c
            out[start:start+64] = np.sum(DD.sigma(er, U, ch)*np.sqrt(er/e)*self.weight, axis=(1, 2))
        return out

    def yield_C(self, m, e, theta, ch):
        lab, w = response(m, float(e))
        return float(w@self.sigma_rate(lab, theta, ch))


def zone_budget(m, e, theta_hot=100.):
    layers, sheet, state = reference(m).profiles(float(e))
    frozen = tuple(tuple(float(x) for x in layer) for layer in layers)
    length, deposited = 0., 0.
    for n, fraction in ((1, .85), (2, .15)):
        for g, w in zip(GX, GW):
            lab = 2*e/n*(1+.05/2.3548*g)
            _, _, _, depth = csda_nodes(float(lab), frozen, transport_stopping)
            parts = _csda_segments(float(lab), frozen, transport_stopping, .1)
            final = float(parts[-1][3](parts[-1][1])[0]) if parts else lab
            deposited += n*fraction*w*(lab-final)*1000
            length = max(length, depth)
    d_column = float(sheet+sum(max(0., min(zb,length)-zt)*nd for zt,zb,_,nd,*_ in layers))
    host_column = float(sum(max(0., min(zb,length)-zt)*nh for zt,zb,nh,*_ in layers))
    current = float(DD.I_cup(e))*1e-6/float(DD.collected_fraction(e, 1., 0.))
    flux_charge = current/QE/AREA
    power = flux_charge*deposited
    return dict(material=m, E_cm_keV=float(e), assumed_area_m2=AREA, current_true_A=current,
        D_number_flux_m2_s=1.15*flux_charge, terminal_zone_depth_m=length,
        fixed_D_column_m2=d_column, host_column_m2=host_column,
        deposited_eV_per_incident_charge=deposited, deposited_W_m2=power*QE,
        deposition_eV_m2_s=power, saved_conditional_temperature_K=float(state['T_K']),
        f100eV_upper_tau100ps=min(1., power*1e-10/(d_column*150.)),
        f1000eV_upper_tau100ps=min(1., power*1e-10/(d_column*1500.)),
        electron_excited_fraction_above1eV_tau1ps_upper=min(1., power*1e-12/host_column),
        electron_density_convention='conditional one active electron per host atom; all deposition to electrons',
        probability_extra_D_r10nm_tau100ps=-math.expm1(-1.15*flux_charge*math.pi*1e-16*1e-10),
        lifetime_for_mean_one_D_r10nm_s=1/(1.15*flux_charge*math.pi*1e-16))


def run():
    rule = Motion(48)
    yields, budgets = {}, []
    settings = {}
    for m in MATERIALS:
        settings[m] = np.unique(np.r_[ENERGIES, 4.5])
        for e in settings[m]:
            temp = reference(m).profiles(float(e))[2]['T_K']
            harmonic = HW/(2*np.tanh(HW/(2*KB*temp)))
            for name, theta in [('cold',0.), ('thermal',KB*temp), ('thermal_ZPM',harmonic),
                                ('hot100',100.), ('hot1000',1000.)]:
                for ch in ('p','n'):
                    yields[m,float(e),name,ch] = rule.yield_C(m,e,theta,ch)
        budgets.extend(zone_budget(m,e) for e in ENERGIES)
    budget_map = {(r['material'],r['E_cm_keV']):r for r in budgets}
    rows = []
    for m in MATERIALS:
        base = lambda e,ch: yields[m,e,'cold',ch]
        rbase = readout(m,ENERGIES,base)
        for name in ('thermal','thermal_ZPM','hot100','hot1000'):
            y = lambda e,ch: yields[m,e,name,ch]
            reported = readout(m,ENERGIES,y)
            for i,e in enumerate(ENERGIES):
                ch = 'p' if e <= 4. else 'n'
                b = budget_map[m,float(e)]
                for power in (0,1):
                    ob,_ = observed(m,ENERGIES)
                    rb = readout(m,ENERGIES,base,current_power=power)
                    q = ob[i]/rb[i]
                    hp = y(e,ch)/base(e,ch)
                    hn = y(6.5,'n')/base(6.5,'n')
                    denom = (hp-1)-q*(hn-1)
                    paired_uniform_f = (q-1)/denom if denom != 0 else float('nan')
                    # Least required occupation: leave the reference unexcited.
                    # Any positive reference heating increases this requirement.
                    f = (q-1)/(hp-1)
                    theta = 100. if name == 'hot100' else 1000. if name == 'hot1000' else None
                    tau = None if theta is None or f < 0 else f*b['fixed_D_column_m2']*1.5*theta/b['deposition_eV_m2_s']
                    rows.append(dict(material=m,E_cm_keV=float(e),state=name,current_power=power,
                        physical_p_factor=y(e,'p')/base(e,'p'),physical_n_factor=y(e,'n')/base(e,'n'),
                        normalized_factor=reported[i]/rbase[i],observed_normalized_excess=q,
                        required_hot_fraction_reference_unexcited=f if theta is not None and e!=6.5 else None,
                        paired_equal_fraction_diagnostic=paired_uniform_f if theta is not None and e!=6.5 else None,
                        fraction_in_0_1=bool(0<=f<=1) if theta is not None and e!=6.5 else None,
                        minimum_lifetime_s_all_deposition_to_D=tau if theta is not None and e!=6.5 else None))
        # Energy-budget capped hot fraction at each nominal energy; no population fitted.
        for theta,name in ((100.,'hot100'),(1000.,'hot1000')):
            fraction = {float(e):zone_budget(m,e)[f'f{int(theta)}eV_upper_tau100ps'] for e in settings[m]}
            def mix(e,ch):
                f=fraction[e]
                return (1-f)*base(e,ch)+f*yields[m,e,name,ch]
            projected = readout(m,ENERGIES,mix)
            for i,e in enumerate(ENERGIES):
                rows.append(dict(material=m,E_cm_keV=float(e),state=f'budget_upper_{name}_100ps',current_power=0,
                    physical_p_factor=mix(e,'p')/base(e,'p'),physical_n_factor=mix(e,'n')/base(e,'n'),
                    normalized_factor=projected[i]/rbase[i],observed_normalized_excess=observed(m,ENERGIES)[0][i]/rbase[i],
                    required_hot_fraction_reference_unexcited=None,paired_equal_fraction_diagnostic=None,
                    fraction_in_0_1=None,minimum_lifetime_s_all_deposition_to_D=None))
    write_csv(OUT/'driven_target_floats.csv',rows)
    write_csv(OUT/'driven_target_budget.csv',budgets)
    # One branch proof: exact cold identity and changed-consumer quadrature comparison.
    identity=max(abs(yields[m,float(e),'cold',ch]/yield_C(m,e,ch)-1)
                 for m in MATERIALS for e in ENERGIES for ch in ('p','n'))
    finer=Motion(72)
    errors=[]
    for e,theta in ((.25,100.),(.25,1000.),(6.5,1000.)):
        coarse=rule.yield_C('PdD',e,theta,'p')
        fine=finer.yield_C('PdD',e,theta,'p')
        errors.append(abs(coarse/fine-1))
    assert identity<2e-13 and max(errors)<.015
    # The 1000-eV case is deliberately broad; record its finite quadrature uncertainty.
    proof=dict(cold_shared_kernel_max_relative=identity,quadrature48_vs72_relative=errors,
        quadrature_decision_tolerance=.015,profile_unchanged=True,source_profile_sha256=hashlib.sha256((ROOT/'g4_beam/merit/profiles.json').read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256((OUT/'driven_target.py').read_bytes()).hexdigest(),
        steady_state_assumption='uniform kinetic mixture over actual traversed fixed-D column including sheet; optimistic eta=1 deposition and chosen energy lifetime')
    (OUT/'driven_target_proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof,indent=2))


if __name__=='__main__':
    run()
