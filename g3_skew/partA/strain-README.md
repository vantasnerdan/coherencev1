# PdD strain/shear collective-stiffness law

PW500/PBE/PAW10, Gamma12, Fermi0.1eV, nonmagnetic, central D displacement +/-0.01A. Reference a0=4.084A. Metals are fixed during displacement; the biaxial table relaxes only c_z at fixed reference fractional sites. These are collective Gamma Born tensors, not local Kii.

## Measured shear families

The actual source table gives the full lab/material tensor, polar R/U and every force column. Lowest eigenvalues at gamma100=2/4/6/8 percent are0.154652354891/0.118712648942/0.078071212258/0.034843198653eV/A2; gamma110=8percent gives0.029644101604.

## Relaxed films

|in-plane strain|c_z/a0|Kxx|Kyy|Kzz (eV/A2)|
|---|---|---|---|---|
|-0.04|1.1025199840929252|0.8684012758994606|0.8683085708314837|-0.34101169339463017|
|-0.02|1.066469433759545|0.39489337024636234|0.3949375505469876|-0.2895514625078555|
|0|1.032461532734177|0.08389307398256204|0.08473493251258385|-0.18066595820475415|
|0.02|1|-0.1269342929521637|-0.12671698078200555|0.057068663903223765|
|0.04|0.9711585943664419|-0.28059003847026015|-0.28060599476345316|0.42063156495110676|

At zero imposed biaxial strain the relaxed c_z/a0=1.032461532734177; its across-layer bare Gamma curvature is -0.18066595820475415, distinct from the cubic reference +0.18105969055575036. Registered-centre Born curvature is not a quantum-stability calculation.

## Callable and achieved fit

```python
import json
from pathlib import Path
from strain_tensor_law import stiffness_tensor, calibration_mask
model = json.loads(Path("strain-law.json").read_text())
K_material = stiffness_tensor(U_minus_I, model)
K_lab = R @ K_material @ R.T
```

The selected cubic-covariant quadratic closure has12 coefficients. Maximum training-entry error0.043972615764eV/A2; independent-entry RMS0.008057042725. Largest identified leave-one-case-out error0.160595968499; removing the single110 case leaves cross-strain coefficients unidentified. These are fit discrepancies, not confidence intervals. The raw families are the primary result; the polynomial is an approximate mixed-strain closure.

Calibration-box normal range[-0.04,0.102519984093], maximum absolute Biot shear0.039968038349. Being inside the box does not certify arbitrary mixed strain. Outside it the function refuses by default; extrapolate=True explicitly requests an extrapolation. There is no eigenvalue clipping. Both transferred EMT cases have0 registered D sites inside this box. The actual-map consumer verifies F=R*(E_law+I)/(a_Pd/a0), and rejects accidental use of the geometry-only strain reference.

Adding a Gamma contrast to a local well assumes unchanged sum(J); absolute Gamma must not replace local Kii or be subtracted again. The EMT geometry has no hydrogen chemical pressure.

## Reproduction and provenance

`python strain_table.py; python fit_strain_law.py; python plot_strain_table.py` reproduces the saved-data fit and actual native-family figure. NumPy is required; plotting also requires Matplotlib. GPAW/ASE and matching PAW10 setups are needed only to regenerate SCFs with film_relax.py/tensor_case.py. Each producer refuses to overwrite existing case directories.

The packet includes all raw central-force JSON, constrained-relaxation stress/energy points, case inputs, execution-profile metadata and source-byte bindings. GPW checkpoints and private Native command logs are omitted. g3-cost-ledger.json binds all35 completed timed GPAW processes, total200.0055core-min. Covariance was exercised over all48 signed cubic permutations, with a cross-shear-index mutation detected; the actual calibrated consumer also refused out-of-box original EMT strains.
