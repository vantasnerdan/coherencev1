"""Conditional Pd80Rh20 H-electrosorption fit; explicitly not alloy DFT or D kinetics."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'screen'))
from screen_io import manifest,require,write_json


def fit(model,primary):
    source=manifest(primary)
    require(source['sha256']=='0ae0ed015f1dc0828c94b4de0048a75f9ce5f9d77da793363d6a94726309e2fc',
            'Primary Pd80Rh20 source edition differs from the consumed source')
    temperature=298.15
    x=np.array([.002,.4,.8]); E_RHE_V=np.array([.050,-.033,-.150])
    # Conditional reversible Nernst interpretation. Reference matches model's half-H2 at1Pa.
    mu=.5*model.KB*temperature*np.log(101325.)-E_RHE_V
    def parameters(v):
        return dict(a300=float(v[0]),b300=float(v[1]),dH=float(v[2]))
    def residual(v):
        return model.Isotherm(model.PD,temperature,alloy=parameters(v)).mu(x)-mu
    result=least_squares(residual,[.05,.62,-21500.],bounds=([.003,.401,-60000.],[.2,.95,1000.]),
                         xtol=1e-12,ftol=1e-12,gtol=1e-12,x_scale=[.05,.6,30000.])
    require(result.success and np.max(abs(residual(result.x)))<1e-8,'PdRh conditional isotherm fit failed')
    overrides=parameters(result.x);isotherm=model.Isotherm(model.PD,temperature,alloy=overrides)
    inverse=np.array([isotherm.x_of_mu(float(value)) for value in mu])
    require(np.max(abs(inverse/x-1))<1e-5,'Independent inverse does not recover source calibration points')
    model300=model.Isotherm(model.PD,300.,alloy=overrides)
    chamber=model300.x_of_mu(.5*model.KB*300*np.log(model.P_GAS))
    return dict(schema='coherence-Pd80Rh20-conditional-isotherm-v1',
        alloy_label='Pd80Rh20 H-electrosorption-fitted isotherm (D transport assumption)',
        material_overrides=overrides,source=source,
        source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC6982334/',doi='10.3390/ma13010162',
        observations=dict(T_K=temperature,E_RHE_V=E_RHE_V.tolist(),H_per_total_metal=x.tolist()),
        conditional_mu_relative_half_H2_at_1Pa_eV=mu.tolist(),
        reproduced_H_per_metal=inverse.tolist(),plateau_pressure_at298K_Pa=float(isotherm.p_pl),
        extrapolated_chamber_D_per_metal_at300K=float(chamber),
        max_mu_calibration_residual_eV=float(np.max(abs(residual(result.x)))),
        retained_estimates={key:model.PD[key] for key in ['dS','Tc','xc','w','tilt','D0','Ea','D0_hyd','Ea_hyd','Eads','s_clean','Eb','mtrap','nM','bca']},
        scope='Three parameters calibrated to three H-electrosorption points: zero holdout degrees of freedom, not predictive validation. Reversible RHE-to-H chemical-potential mapping is conditional. H-to-D isotope, kinetics, surface adsorption, implantation/density and temperature extrapolation retain Pd estimates; no PdRh DFT values supplied.',
        composition_warning='Pd80Rh20 foil data are not the ordered Pd3Rh25% screening cell. Do not combine them into a causal d-band derivative.',
        source_measurement='Total oxidised surface+absorbed H in50um foil at25C. Applied potentials after reported saturation, not direct gas-pressure measurements.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--transport',type=Path,required=True)
    p.add_argument('--primary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();sys.path.insert(0,str(a.transport.resolve().parent))
    s=importlib.util.spec_from_file_location('transport_fit_reference',a.transport);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    result=fit(m,a.primary);result['transport_source']=manifest(a.transport);result['fit_source']=manifest(Path(__file__))
    write_json(a.output,result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
