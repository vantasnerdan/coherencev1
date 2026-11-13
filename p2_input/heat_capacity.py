"""Part A: a source-bound estimate, not new heat-capacity measurements.

PdD_x reference cell: Antonov et al., PRB80,134302(2009),p1, a=4.084A.
D optical peak/shoulder: same source,p3,37/65meV.
Pd Cp=25.98 J/mol/K, rho=11995kg/m3,M=106.4g/mol at298.15K:
Michael de Podesta, Understanding the properties of matter,Table7.10,p1.
PdH Debye theta=246K, pure Pd=276K: ISSP Yamamuro2012 activity report.
Einstein and Debye harmonic heat-capacity formulas; NIST SI constants.
"""
import json
import math
from scipy.integrate import quad

KB = 1.380649e-23
EV_J = 1.602176634e-19
NA = 6.02214076e23
T = 300.0
LOADING = 0.8
N_PD = 4/(4.084e-10)**3


def einstein_D(energy_eV):
    y = energy_eV*EV_J/(KB*T)
    factor = y*y*math.exp(y)/math.expm1(y)**2
    return 3*LOADING*N_PD*KB*factor


def debye_Pd(theta_K):
    upper = theta_K/T
    integral = quad(lambda t: t**4*math.exp(t)/math.expm1(t)**2,0,upper)[0]
    return 9*N_PD*KB*(T/theta_K)**3*integral


pd_classical = 3*N_PD*KB
pd_cp_proxy = N_PD/NA*25.98
point = pd_classical+einstein_D(0.037)
lower_model = debye_Pd(276)+einstein_D(0.065)
upper_proxy = pd_cp_proxy+einstein_D(0.037)
print(json.dumps(dict(
    temperature_K=T,loading=LOADING,n_Pd_m3=N_PD,
    Einstein_temperature_K=0.037*EV_J/KB,
    shoulder_Einstein_temperature_K=0.065*EV_J/KB,
    Pd_DulongPetit_J_m3_K=pd_classical,D_Einstein_J_m3_K=einstein_D(0.037),
    Cv_point_J_m3_K=point,
    illustrative_model_range_J_m3_K=[lower_model,upper_proxy],
    warning="Range is model/proxy spread, not measured PdD Cv or a confidence interval; upper proxy uses Cp, not Cv.",
    pure_Pd_Cp_volume_J_m3_K=25.98*11995/0.1064,
    Pd_lattice_Debye_J_m3_K={str(theta):debye_Pd(theta) for theta in (246,276)},
),indent=2))
