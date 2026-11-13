"""Coherencev1 p0_scales: published reference cells, SI arithmetic, no simulation.

Sources:
NiH: Antonov et al., Europhys. Lett. 51, 140 (2000), p142:
     a=3.732 A at100 K; first optical peak89 meV at5 K, shoulder to130 meV.
PdD: Antonov et al., Phys. Rev. B80,134302 (2009), pp1,3:
     a=4.084 A (quoted stoichiometric PdD reference), peak37 meV, shoulder65 meV.
Host optical fits: Rakic et al., Appl. Opt.37,5271 (1998), Tables1-3.
Constants: NIST2022 CODATA complete listing. Bare ionic charge=+e.

Fixed a and phonon energies across x and T are SCENARIOS, not measured a(x,T).
All frequencies denoted omega are angular; E=hbar*omega in eV.
All wavelengths are VACUUM REFERENCES, not metal propagation/domain predictions.
"""
import json
from math import pi, sqrt

E_CHARGE = 1.602176634e-19  # C; also J/eV
PLANCK = 6.62607015e-34  # J s
C_LIGHT = 299792458.0  # m/s
KB = 1.380649e-23  # J/K
EPS0 = 8.8541878188e-12  # F/m
HBAR = PLANCK / (2 * pi)
HC_EV_UM = PLANCK * C_LIGHT / E_CHARGE * 1e6
# Independent NIST inverse-metre/eV conversion, rounded only as tabulated.
NIST_HC_EV_UM = 1.239841984
MATERIALS = [
    dict(material="Ni/H", a_A=3.732, mass_kg=1.67262192595e-27,
         vibration_eV=0.089, band_peak_to_shoulder_eV=[0.089, 0.130],
         host_optical_wp_eV=15.92, host_f0=[0.083, 0.096],
         host_electron_damping_eV=[0.022, 0.048]),
    dict(material="Pd/D", a_A=4.084, mass_kg=3.3435837768e-27,
         vibration_eV=0.037, band_peak_to_shoulder_eV=[0.037, 0.065],
         host_optical_wp_eV=9.72, host_f0=[0.330, 0.330],
         host_electron_damping_eV=[0.008, 0.009]),
]


def main():
    rows = []
    for material in MATERIALS:
        n_metal = 4 / (material["a_A"] * 1e-10) ** 3
        evib = material["vibration_eV"]
        for loading in (0.6, 0.8, 1.0):
            n_ion = loading * n_metal
            omega = sqrt(n_ion * E_CHARGE**2 / (EPS0 * material["mass_kg"]))
            energy = HBAR * omega / E_CHARGE
            wavelength = 2 * pi * C_LIGHT / omega * 1e6
            # E*lambda has units eV*um; compare once below with NIST.
            rows.append(dict(
                material=material["material"], loading=loading,
                a_A=material["a_A"], n_metal_m3=n_metal, n_ion_m3=n_ion,
                ion_omega_rad_s=omega, ion_plasma_eV=energy,
                ion_plasma_vacuum_wavelength_um=wavelength,
                vibration_eV=evib,
                host_optical_wp_eV=material["host_optical_wp_eV"],
                host_intraband_plasma_eV=[material["host_optical_wp_eV"] * sqrt(f)
                                         for f in material["host_f0"]],
                kBT_300_eV=KB * 300 / E_CHARGE,
                kBT_600_eV=KB * 600 / E_CHARGE,
                illustrative_Ec_Egap_eV=evib,
                illustrative_rho=energy / evib,
                illustrative_Dicke_g_f1=energy / (2 * evib),
                illustrative_diamagnetic_d=(energy / evib) ** 2 / 4,
                illustrative_detuning=0.0,
                illustrative_temperature_300=KB * 300 / E_CHARGE / evib,
                illustrative_temperature_600=KB * 600 / E_CHARGE / evib,
                illustrative_mode_vacuum_wavelength_um=HC_EV_UM / evib,
            ))
    unit_reference = rows[0]
    print(json.dumps(dict(
        rows=rows, host_reference_parameters=MATERIALS,
        parametric_mapping=dict(rho="Epi/Ec", r="Egap/Ec", detuning="r-1",
                                Dicke_g="rho*sqrt(f*r)/2", TRK_d="rho**2/4",
                                temperature="kBT/Egap",
                                vacuum_wavelength_um="1.239841984/Ec[eV]",
                                losses="kappa_field/omega_c, gamma_Bloch/omega_c",
                                loss_guess_each=[1e-3, 1.0]),
        unit_crosscheck=dict(
            calculated_eV_um=(unit_reference["ion_plasma_eV"] *
                              unit_reference["ion_plasma_vacuum_wavelength_um"]),
            NIST_eV_um=NIST_HC_EV_UM),
    ), indent=2))


if __name__ == "__main__":
    main()
