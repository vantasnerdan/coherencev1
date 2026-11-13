# g1_rank/pumps.py
"""Physical pump feasibility, no arbitrary stimulated-gain cross section.
Reuses hybrid-normalized rank.py and domain.py. Run: python pumps.py > pumps.json
Decision: chemical arrival, coherent resonant/parametric light, or pulsed current?
NumPy and standard library; negligible compute, one thread.
"""
import json
import numpy as np
from domain import DEFAULT, HBAR, HC, MH, KB, impedance, parameters
from rank import MATS, gain_branch

CHARGE = 1.602176634e-19
EPS0 = 8.8541878188e-12
C = 299792458.
KB_J = KB*CHARGE
MH_KG = 1.67262192595e-27


def optical(m):
    g = gain_branch(m)
    z = parameters(m, DEFAULT | {'force':0.})
    nu = g['mode_frequency_THz']*2*np.pi*1e12*HBAR/m['E']
    qd = np.pi/2
    kq = z['onsite']-z['Jperx']*z['x']+4*DEFAULT['gradient']*np.sin(qd/2)**2
    D = impedance(nu,kq,z['g'],DEFAULT['photon_ratio'],
                  DEFAULT['electron_strength'],z['we'],DEFAULT['gamma'])
    target = .1
    force = target*abs(D+.75*DEFAULT['beta']*target*target)
    l0 = HC/np.sqrt(m['mass']*m['E'])
    Eeff = force*m['E']/l0*1e10  # V/m for effective charge +e
    Iref = .5*EPS0*C*Eeff*Eeff
    duration = 3*g['after_pump_amplitude_decay_ps']*1e-12
    Omega = g['mode_frequency_THz']*2*np.pi*1e12
    kappa = g['passive_energy_decay_s_inv']
    metric = g['hybrid_action_metric']
    # Modulation of the ACTUAL potential kq, not an independently scaled spring.
    delta_k_th = 2*kappa*nu*metric/(m['E']/HBAR)
    hbare = delta_k_th/kq
    return dict(material=m['name'], carrier_THz=g['mode_frequency_THz'],
        coherent_resonant_onset_field_V_m=0., target_displacement_A=target*l0,
        target_internal_effective_field_V_m=float(Eeff),
        incident_field_scaling='divide by abs(Z_mode/e * internal_field_transmission)',
        target_unit_coupling_intensity_W_m2=float(Iref),
        target_pulse_duration_ps=float(duration/1e-12),
        target_unit_coupling_fluence_mJ_cm2=float(Iref*duration*.1),
        parametric_pump_THz=2*g['mode_frequency_THz'],
        parametric_fractional_effective_stiffness_threshold=float(2*kappa/Omega),
        parametric_fractional_bare_stiffness_threshold=float(hbare),
        parametric_reduced_delta_k_threshold=float(delta_k_th),
        reference_strain_derivative_dlnK_deps=5.,
        parametric_strain_threshold_for_derivative5=float(hbare/5),
        parametric_threshold_field_V_m=None,
        missing_field_conversion='mode symmetry and Raman/IR/strain derivative mapping pump into stiffness modulation',
        parametric_phase='frequency locked to half pump; two phase states, not a free-phase laser linewidth')


def chemical_pd():
    # PdH estimated by mass-rescaling the measured PdD mode (not a new measurement).
    pd = MATS[2]
    gd = gain_branch(pd)
    scale = np.sqrt(pd['mass']/MH)
    kappa = gd['passive_energy_decay_s_inv']*scale
    eph = gd['mode_frequency_THz']*2*np.pi*1e12*HBAR*scale
    nm = 4/pd['a']**3*1e30
    thickness = {'monolayer':pd['exclusion_A']*1e-10,
                 'gain_length_patch':gd['gain_length_nm']*1e-9}
    atom_flux_per_Pa = 2/np.sqrt(2*np.pi*(2*MH_KG)*KB_J*300.)
    rows = []
    for label,L in thickness.items():
        # Direct population inversion in a recycled two-bound-level idealization:
        # rate of hot arrivals per site must exceed kappa/2.
        max_hot_probability = min(1., .1/eph)
        flux = .5*kappa*nm*L/max_hot_probability
        energy_floor = .5*kappa*nm*L*eph/.1
        rows.append(dict(layer=label, thickness_nm=float(L/1e-9),
            H_arrival_flux_inversion_lower_bound_m2_s=float(flux),
            chemical_power_energy_only_flux_bound_m2_s=float(energy_floor),
            ideal_H2_pressure_at300K_bar=float(flux/atom_flux_per_Pa/1e5),
            chemical_power_W_m2=float(flux*.1*CHARGE)))
    return dict(material='matched PdH from PdD force-constant reference',
        isotope_mode_eV=float(eph), passive_energy_decay_s_inv=float(kappa),
        solution_heat_eV_per_H=.1,
        atomic_H_extra_energy_eV_per_H=2.24,
        energy_reference='0.1eV heat of solution is already relative to half H2; do not add H2 dissociation again. Prepared atomic H can release an additional2.24eV/H; preparing it paid that energy.',
        chemical_inversion_status='optimistic hot-arrival/recycled finite-level bound; thermalizing a harmonic ladder alone never inverts it',
        sticking_and_hot_branching_scaling='actual threshold divides by hot-state branching and sticking probabilities',
        rows=rows)


if __name__ == '__main__':
    gs = [gain_branch(m) for m in MATS]
    current = [dict(material=g['material'],
        existing_optic_drift_floor_A_m2=g['current_kinematic_A_m2'],
        requested_pulse_range_A_m2=[1e12,1e13],
        largest_pulse_over_drift_floor=1e13/g['current_kinematic_A_m2'],
        result='below necessary finite-q optic drift condition; no new steady-current gain model') for g in gs]
    print(json.dumps(dict(
        choice='parametric optical stiffness modulation for the self-sustained branch; resonant mid-IR response calibrates coupling and loss first',
        chemical=chemical_pd(), optical=[optical(m) for m in MATS],
        pulsed_or_constriction_current=current,
        corrected_hybrid_gain=gs,
        source_benchmark=dict(material='LiNbO3, not metal hydride',
            carrier_THz=20., duration_fs=150., peak_field_MV_cm=30.,
            max_fluence_mJ_cm2=200.,
            url='https://doi.org/10.1038/s41567-022-01512-3'),
        missing_physical_inputs='loaded-hydride mode effective charge/field penetration, stiffness-modulation derivative, hot-arrival branching, and actual linewidth; none is replaced by a geometric gain cross section'),
        indent=2, allow_nan=False))

