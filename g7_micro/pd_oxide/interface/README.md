# Pd interface readout

Actual radii (A): [0.1, 0.2, 0.4]; missing: []. Only actual common-R saved buried-site numbers are compared; no new bare calculation, baseline replay, interpolation or far stitching.

D force: 0.18832834085029737 eV/A (partly-relaxed); nearest O/Pd: 1.5242634583865406/1.5560270889516299 A. All 23 non-D atoms are frozen; no full-slab stationarity or force-to-energy error bound.

Assumptions: PBE finite differences h=0.20 A, Gamma-only, 0.10 eV Fermi width, fixed coated step-5 Pd/O frame, bare q+1 probe with no added electron and the precommitted ray.

curve.json retains original quartet energies; internal = free - negative TS. contact.json reports signed W0 = QHe + Bnative - Bsource + regular images. Contact is compared to V(0.1)-C/0.1, NOT the positive lowering; its signed difference, physical electronic/nuclear decomposition and scale are reported without an equality tolerance.

Vacancy: Reused original vacancy: exact selected nonD frame match; already charged.

Original state energies and integer CPU accounting are retained, not recertified. Machine paths are redacted to basenames with existing hashes. No GPWs, logs, source trees or support programs are bundled.

Stage charge: 171.4532733333 Native core-min, 5 preparation core-min and 0.9327458337 measured collector core-min = 177.3860191670 of 270; 200 reserve remains held. This includes the initial capped geometry attempt, which is outside the completed-case unit sum. The later geometry driver's actual exit-code finalizer failure is retained in the unit accounting; its already-published complete frame was used without repeating SCF.

D-only interface geometry in the exact frozen coated step-5 host; frozen-host forces are not slab stationarity. Finite-cell q+1/no-added-electron internal BO near curve and signed regular contact only. No geometry-error bound, global minimum, site stability, isolated limit, full h/k/cell convergence, Ueq, action/Eq7 or experimental phase claim.
