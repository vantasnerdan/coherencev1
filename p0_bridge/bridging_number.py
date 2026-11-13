"""p0_bridge: conditional requirements from reported heat, not a LENR validation.
Claims: Fleischmann, Pons & Hawkins (1989), J. Electroanal. Chem.261,301-308,
p304 Table1: 0.095 and21.4 W/cm3. Read original table, including length footnote.
Densities: Coherencev1 p0_scales.
Ordering energy: p0_domain result, 1-10 meV/site with the diamagnetic term off.
Q is assumed fully thermalized; no reaction, domain geometry or lifetime assumed.
C is the fraction of occupied H/D sites in the REFERENCE volume that are coherent.
"""
import json

EV_J = 1.602176634e-19  # exact SI J/eV, NIST
DENSITIES = [["Ni/H",0.6,4.617281326359789e+28],["Ni/H",0.8,6.156375101813054e+28],["Ni/H",1,7.695468877266316e+28],["Pd/D",0.6,3.5233358392018157e+28],["Pd/D",0.8,4.6977811189357545e+28],["Pd/D",1,5.872226398669693e+28]]  # [material, loading, occupied-site m^-3]
CLAIMS = [("FPH1989_p304_1mm_8mAcm2", 0.095),
          ("FPH1989_p304_4mm_512mAcm2", 21.4)]  # W per full-electrode cm3
Q_MeV = (1.0, 25.0)
ORDERING_ENERGY_meV = (1.0, 10.0)


def required_rate(power_W_cm3, q_MeV, n_m3, coherent_fraction=1.0):
    power_W_m3 = power_W_cm3 * 1e6
    q_J = q_MeV * 1e6 * EV_J
    return power_W_m3 / (q_J * n_m3 * coherent_fraction)


def main():
    all_n = [row[2] for row in DENSITIES]
    pd_n = [row[2] for row in DENSITIES if row[0] == "Pd/D"]
    rows = []
    for claim, power in CLAIMS:
        for q in Q_MeV:
            rows.append(dict(
                claim=claim, reported_W_cm3=power, Q_MeV=q,
                event_rate_m3_s=power * 1e6 / (q * 1e6 * EV_J),
                C1_rate_full_density_range_s1=[required_rate(power,q,max(all_n)),
                                              required_rate(power,q,min(all_n))],
                C1_rate_PdD_density_range_s1=[required_rate(power,q,max(pd_n)),
                                             required_rate(power,q,min(pd_n))],
            ))
    ratios = [dict(Q_MeV=q, ordering_energy_meV=list(ORDERING_ENERGY_meV),
                   Q_over_abs_du_equals_energy_equivalent_sites=[
                       q * 1e6 / (max(ORDERING_ENERGY_meV) * 1e-3),
                       q * 1e6 / (min(ORDERING_ENERGY_meV) * 1e-3)]) for q in Q_MeV]
    reference_n = next(row[2] for row in DENSITIES
                       if row[0] == "Pd/D" and row[1] == 1.0)
    reference_rate = required_rate(21.4,25.0,reference_n)
    print(json.dumps(dict(
        inputs=dict(densities=DENSITIES,claims_W_cm3=CLAIMS,C=1.0,
                    Q_MeV=Q_MeV,assumption="all Q becomes measured heat"),
        rows=rows,energy_ratios=ratios,
        reference=dict(reported_W_cm3=21.4,Q_MeV=25.0,n_m3=reference_n,
                       C1_rate_s1=reference_rate),
        scaling="every listed per-coherent-site rate scales as 1/C",
        energy_ledger_reference_W_cm3=(reference_rate * reference_n *
                                      (25.0 * 1e6 * EV_J) / 1e6),
    ),indent=2))


if __name__ == "__main__":
    main()
