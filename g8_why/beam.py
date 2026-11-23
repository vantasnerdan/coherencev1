"""Fixed-profile species, above-nominal tail and conditional shared-wake estimates."""
import hashlib
import json
import math
import numpy as np
from scipy.optimize import least_squares
from common import (ROOT, OUT, DD, QE, U, ENERGIES, MATERIALS, GX, GW,
                    response, yield_C, readout, observed, reference, write_csv)
from g4_beam.forward.consumer import forward_yields
from g4_beam.forward.profile_thermal import transport_stopping

SOURCES = {'D1_only':(1.,0.,0.), 'D2_85percent':(.15,.85,0.),
           'D2_only':(0.,1.,0.), 'D3_15percent_sensitivity':(.85,0.,.15),
           'D3_only_sensitivity':(0.,0.,1.)}
TAILS = (2., 5., 13.)  # per-D lab keV; hypotheses, not measured low-setting tail


def timing_rows():
    """Conditional 1-Angstrom separation at the target, not measured geometry."""
    result=[]
    for e in ENERGIES:
        for n in (2,3):
            lab=2*e/n
            v=math.sqrt(2*lab*1000*QE/(DD.M_D*1.66053906660e-27))
            result.append(dict(E_cm_nominal_keV=float(e),multiplicity=n,per_D_lab_keV=lab,
                assumed_D_separation_A=1.,maximum_longitudinal_arrival_gap_fs=1e-10/v*1e15,
                maximum_transverse_offset_A=1.,geometry='conditional 1A separation at target; breakup/orientation not measured'))
    return result


def run():
    grid = np.unique(np.r_[ENERGIES, 4.5])
    base, species, neutral = {}, {}, {}
    for m in MATERIALS:
        for e in grid:
            for ch in ('p','n'):
                base[m,float(e),ch]=yield_C(m,e,ch)
                for tag,f in SOURCES.items():
                    species[m,float(e),ch,tag]=yield_C(m,e,ch,fractions=f)
                for lab in TAILS:
                    neutral[m,float(e),ch,lab]=yield_C(m,e,ch,incident_lab=lab)
    rows, demands, wake, timing = [], [], [], []
    for m in MATERIALS:
        y0=lambda e,ch:base[m,e,ch]
        rb=readout(m,ENERGIES,y0)
        for tag,f in SOURCES.items():
            y=lambda e,ch:species[m,e,ch,tag]
            r=readout(m,ENERGIES,y)
            for i,e in enumerate(ENERGIES):
                rows.append(dict(material=m,E_cm_keV=float(e),hypothesis=tag,
                    physical_proton_factor=y(e,'p')/y0(e,'p'),physical_neutron_factor=y(e,'n')/y0(e,'n'),
                    normalized_factor=r[i]/rb[i],predicted_own65_R=r[i],
                    D_atoms_per_charged_C=sum((j+1)*f[j] for j in range(3))/QE))
        for tail in TAILS:
            for delta in (0,1):
                # delta1: constant neutral rate referenced to internal5kV current.
                scale=lambda e:1. if delta==0 else DD.I_INT/(DD.I_cup(e)/DD.collected_fraction(e,1.,0.))
                # A source-specific Czerski control scale is only transplanted sensitivity.
                ftest=3e-4
                y=lambda e,ch:y0(e,ch)+ftest*scale(e)*neutral[m,e,ch,tail]
                rr=readout(m,ENERGIES,y)
                for i,e in enumerate(ENERGIES):
                    rows.append(dict(material=m,E_cm_keV=float(e),hypothesis=f'neutral{tail:g}lab_f3e4_delta{delta}',
                        physical_proton_factor=y(e,'p')/y0(e,'p'),physical_neutron_factor=y(e,'n')/y0(e,'n'),
                        normalized_factor=rr[i]/rb[i],predicted_own65_R=rr[i],D_atoms_per_charged_C=(1.15+ftest*scale(e))/QE))
                for power in (0,1):
                    obs,_=observed(m,ENERGIES)
                    rbase=readout(m,ENERGIES,y0,current_power=power)
                    for i,e in enumerate(ENERGIES):
                        q=obs[i]/rbase[i]
                        ch='p' if e<=4 else 'n'
                        h=scale(e)*neutral[m,float(e),ch,tail]/y0(e,ch)
                        hr=scale(6.5)*neutral[m,6.5,'n',tail]/y0(6.5,'n')
                        denom=h-q*hr
                        fpaired=(q-1)/denom if denom!=0 else None
                        lower=(q-1)/h
                        count=None
                        if e in DD.SI1_E:
                            count=float(DD.counts_per_C_table(m)[list(DD.SI1_E).index(e)])
                        demands.append(dict(material=m,E_cm_keV=float(e),tail_lab_keV=tail,
                            neutral_rate_convention=delta,current_power=power,
                            neutral_D_per_charge_lower_if_reference_tail_zero=lower if e!=6.5 else None,
                            neutral_D_per_charge_paired_model=fpaired if e!=6.5 else None,
                            paired_nonnegative=bool(fpaired is not None and fpaired>=0) if e!=6.5 else None,
                            neutral_reactions_per_C_per_one_D_per_charge=neutral[m,float(e),ch,tail],
                            SI1_required_f_times_unknown_eta_p=None if count is None else count*DD.collected_fraction(e,1.,0.)/neutral[m,float(e),'p',tail]))
        # Conditional upper-scale test only: pretend a molecule's second encounter
        # gets a rigid positive .22eV shift at all local energies, unlike actual g6.
        boosted={}
        for e in grid:
            lab,w=response(m,float(e),(0.,1.,0.))
            for ch in ('p','n'):
                boosted[float(e),ch]=float(w@DD.sigma(lab/2,U+.00022,ch))
        y=lambda e,ch:base[m,e,ch]+.15*.5*(boosted[e,ch]-species[m,e,ch,'D2_only'])
        rw=readout(m,ENERGIES,y)
        for i,e in enumerate(ENERGIES):
            beta=.15*species[m,float(e),'p','D2_only']/base[m,float(e),'p']
            q=observed(m,ENERGIES)[0][i]/rb[i]
            wake.append(dict(material=m,E_cm_keV=float(e),D2_proton_yield_share=beta,
                both_constituents_cross_section_multiplier_needed_ref_unchanged=1+(q-1)/beta if e!=6.5 else None,
                second_constituent_multiplier_needed_ref_unchanged=1+2*(q-1)/beta if e!=6.5 else None,
                conditional_g6_022eV_second_constituent_physical_factor=y(e,'p')/base[m,float(e),'p'],
                conditional_g6_022eV_second_constituent_normalized_factor=rw[i]/rb[i]))
    fit_rows=[]
    # Eight small one-parameter fits: two energies/rate/current conventions.
    for tail in (5.,13.):
        for delta in (0,1):
            for power in (0,1):
                def evaluate(x,m,es):
                    f=float(np.exp(x[0]))
                    scale=lambda e:1. if delta==0 else DD.I_INT/(DD.I_cup(e)/DD.collected_fraction(e,1.,0.))
                    y=lambda e,ch:base[m,e,ch]+f*scale(e)*neutral[m,e,ch,tail]
                    return readout(m,es,y,current_power=power)
                def residual(x):
                    return np.concatenate([(np.log(evaluate(x,m,ENERGIES[:-1]))-np.log(observed(m,ENERGIES[:-1])[0]))/
                        observed(m,ENERGIES[:-1])[1] for m in MATERIALS])
                fit=least_squares(residual,[-6.],bounds=(-40.,15.),max_nfev=80)
                fit_rows.append(dict(tail_lab_keV=tail,neutral_rate_convention=delta,current_power=power,
                    fitted_neutral_D_per_charge=float(np.exp(fit.x[0])),chi2=float(fit.fun@fit.fun),n=16,
                    assumptions='common fraction, fixed p/n calibration; hypothetical tail at low gun settings'))
    timing=timing_rows()
    optical=[]
    for torr in (2.5e-4,3.e-4):
        column=torr*133.322368*.04/(1.380649e-23*300.)
        for xs in (1e-20,1e-19):
            p=-math.expm1(-column*xs)
            optical.append(dict(pressure_Torr=torr,path_m=.04,gas_temperature_K=300.,
                all_gas_assumed_D2=True,assumed_neutralization_cross_section_m2=xs,
                column_m2=column,conversion_probability=p,neutral_atoms_per_surviving_atomic_charge_upper=p/(1-p),
                energy_qualification='conversion cannot exceed pre-neutralization ion energy; no documented keV upstream-deceleration path at .5keV setting'))
    write_csv(OUT/'beam_floats.csv',rows)
    write_csv(OUT/'beam_neutral_requirements.csv',demands)
    write_csv(OUT/'beam_wake.csv',wake)
    write_csv(OUT/'beam_timing.csv',timing)
    write_csv(OUT/'beam_optical_depth.csv',optical)
    # Settled branch proof: atom-current multiplicity and actual neutral consumer.
    species_identity=max(abs((.85*species[m,float(e),ch,'D1_only']+.15*species[m,float(e),ch,'D2_only'])/base[m,float(e),ch]-1)
        for m in MATERIALS for e in ENERGIES for ch in ('p','n'))
    m,e='PdD',.25
    native_neutral=0.
    for g,w in zip(GX,GW):
        val=forward_yields([e],reference(m).profiles,Ue_keV=U,candidate='fast_neutrals',
            parameters={'f_n':1.,'E_n_lab_keV':13.*(1+.05/2.3548*g)},
            measurement={'channel':'p'},stopping_fn=transport_stopping,reaction='eq7',reference_E_cm_keV=None)
        native_neutral+=w*val['neutral_yield_per_C'][0]
    neutral_identity=abs(native_neutral/neutral[m,e,'p',13.]-1)
    assert species_identity<3e-13 and neutral_identity<3e-13
    proof=dict(atom_current_species_identity_max_relative=species_identity,
        neutral_changed_consumer_relative=neutral_identity,neutral_energy_distribution='5%FWHM Gaussian around declared per-D lab energy',
        source_profile_sha256=hashlib.sha256((ROOT/'g4_beam/merit/profiles.json').read_bytes()).hexdigest(),
        common_helper_sha256=hashlib.sha256((OUT/'common.py').read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256((OUT/'beam.py').read_bytes()).hexdigest())
    (OUT/'beam_proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    (OUT/'beam_tail_fits.json').write_text(json.dumps(fit_rows,indent=2)+'\n')
    print(json.dumps({'proof':proof,'common_tail_fits':fit_rows},indent=2))


if __name__=='__main__':
    run()
