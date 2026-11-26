# Pd monovacancy with six D: relaxed 27-site cell

This capsule provides the force-qualified neutral PBE structure of one Pd monovacancy with six D in a fixed 3×3×3 primitive fcc repeat (26 Pd plus 6 D). The maximum unconstrained residual force is 0.0451522 eV/Å, below the 0.05 eV/Å target. All six D remain in the registered vacancy cage. The 15 periodic D–D distances range from 1.959332 to 3.289906 Å; the complete labelled distances, atom order, positions and forces are in `native/reference.json`. D labels are zero-based in their saved order. The nearest D–D distance to a different periodic vacancy image is 5.587228 Å; this geometric clearance does not by itself bound the finite-cell error.

This endpoint identifies one relaxed basin of the chosen cell and seed. It does not establish a global minimum, thermal occupancy, dynamic screening, or a large-cell resting separation. No vacancy pair cut is triggered by the small-cell result alone.

`native/saved-start/step-0002.json` binds the exact position used to start this relaxation; the electronic state and optimizer history were recomputed. `sha256sum -c SHA256SUMS` verifies the member bytes. `replay/pd_multivac.py` is the exact producer; literature and PAW source hashes are in the native records and the external binaries are not redistributed.
