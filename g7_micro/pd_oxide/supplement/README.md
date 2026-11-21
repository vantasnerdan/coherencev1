# Compact Pd/PdO geometry supplement

Base capsule: private draft453374618, archive774098833, SHA25697d55c55be5bb9151b5310a8ff1b0b480933279b425deba8bd4191e731895bb2; integrated by Axiom on coherencev1 main as a31c95e. The base is unchanged.

## Actual base D sites

Signed heights use the mean actual z of the original top substrate Pd layer (indices0–4) and of the added4Pd+4O sheet (indices15–22). Negative height is below the reference plane. Distances use minimum images in the saved periodic cell.

- **Bare Pd:** D z=9.4466586004A, height relative to top-metal mean=-2.0261052585A; no oxide sheet or O atoms (oxide-relative height/nearest O not applicable); nearest Pd=1.9466586004A.
- **Coated Pd:** D z=9.5589173814A, height relative to top-metal mean=-1.9036039085A and to oxide-sheet mean=-4.4056927608A; nearest O=4.1898073130A, nearest Pd=1.8770967666A. The target is buried in the metal beneath the oxide, **not** directly at the oxide/metal plane.

The saved coordinates are in geometry-coordinate-ledger.json; unrounded metrics and original-geometry byte/SHA bindings are in site-geometry.json. Independent NumPy minimum-image and ASE find_mic results agree within1.8e-15A; no Native calculation was performed for these geometry lines.

Both geometry-error bounds remain null. Bare free fmax=0.03603609470523719eV/A; coated originalstep5 free fmax=0.42174246367039653eV/A is partly relaxed. Paired k and h checks were omitted under entry3065; no equilibrium, site-stability or experimental claim follows.

## Second-frame R0.1 attempt: unfinished

Original accepted step3 was chosen rather than step6 because it is more separated from base step5: D displacement0.07236025933371798A and added-sheet RMS displacement0.05900090459337368A (maximum0.07962165997627553A). Its free fmax is0.621541089029348eV/A; this is not force-qualified.

The occupied host completed84 SCF iterations at215 electrons; its actual unit used1008485201000ns. The vacant host reached its1079s manager limit, using1078449841000ns, without a completed state JSON, GPW, boundary or SCF-complete marker. Both probe states were never started. The four-state quartet is incomplete: **no second-frame lowering or two-frame spread exists**; both are explicitly null in second-frame-attempt.json. No second-frame contact was attempted.

Actual cost34.7822507 of the80core-min ceiling is retained. A conservative forecast from measured host/cold-probe costs exceeds the remaining45.2177493core-min and presses03:00; this is an allocation decision, not proof that every recovery route is impossible. Stop accepted: no retry, checkpoint import, extra probe, h/k or other Native work in this stage. An operational timeout is not proof of scientific nonconvergence or a physical no-go. The completed base capsule remains unchanged; part(b) stays open for the lead's next allocation decision.
