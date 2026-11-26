# G13 Pass 2 — live whole-Pd-foil draft

This folder is the only changed surface; accepted base `9c26f8cb688918ae18b63dbd187fdf99793c1e90`. No Pass 1 replay, Ti, electronic calculation, strain-to-D calibration, publication or commit.

Governing steering is 4097, 4103 **as corrected by 4104**, and structured-beam addition 4109. A rogue hole means a transient signed trough, **not damage**. The fixed-source draft is retained as a baseline only; its conclusions are not substituted for the corrected live arms.

## Keep / drop / backwards dependencies

**Keep:** a 250-μm-thick, 21.269-mm-square free foil with the same area and mass as the actual 24-mm disk, centrally driven over 1 cm². Vector Ritz ∫BᵀCB propagation retains cubic anisotropy, explicit rotated grains, both faces and lateral-edge reflections. New launches depend predictably on the field already present before each physical-time bin and carry an actual distribution of bounded marked energies, not a fixed mean packet. Spatially peaked/hot-spot and temporally correlated burst/ripple beams are compared at matched charge against uniform Poisson. Separate impact-overlap, off-impact directional-convergence and pre-existing signed peak/hole impingement statistics retain their own denominators. Signed holes are scored whether struck or not; actual quadratic strain energy remains separate.

**Drop:** atomistic cascades, Ti, persistent damage/dose, new electronic/force or pair calculations, unrepresented short wavelengths, measured grain-map claims and calibrated strain-to-D transfer. Pure-Pd constants retain their secondary-transcribed near-300-K provenance; the loaded population profile is not a measured PdD tensor. The square/disk/flange exterior, band-launch law and material loss remain declared parameters. The live classes do not use the retired roughness interpretation.

**Work backwards:** chart current × accepted host stopping × unknown band-launch fraction → signed source → actual finite variational propagation → region/depth/time strain-energy tails → acoustic storage per local D → separately conditional local attenuation-to-D conversion and holding time → retained 0.6423 meV/D *reference*, not an instability threshold for x<1. Source concentration, attenuation, nonlinear strength and holding time are separate axes; a joint frontier must state its normalization, not claim a universal cheapest combination.

**First decision:** how the three built event populations and symmetric positive/negative tails change with live-state launch feedback, structured beam, lifetime and propagation range. “Common” means 50% with a named impact or receiving-window denominator. No class is dismissed before propagation. A/C overlap is counted without duplicating their union; B is an off-impact receiving event, not added to an impact rate.

## Run

From repository root, with NumPy/SciPy available:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python g13_rogue/pass2/coarse.py --live --exposure 0.0001 --realizations 2 --output g13_rogue/pass2/output
```

The corrected live outputs are under `output/live/`; equations and event/beam semantics are in [live_method.md](live_method.md). `output/` without `live/`, [fixed_source_readout.md](fixed_source_readout.md), and [method.md](method.md) retain the already-computed fixed-source comparison and conservative nonlinear witness, not a replacement for the three live classes. The exact inherited scientific inputs are in [inputs_ledger.md](inputs_ledger.md).

## Fixed-source baseline, not the corrected live result

The retained [fixed_source_readout.md](fixed_source_readout.md) records a conservative localized nonlinear model and the separate direct attenuation/holding map. Those results are conditional on its fixed source; they are not claims about the later live launch, rogue-hole or clumped-beam populations. The corrected numerical readout is in [readout.md](readout.md).

The propagated source is conditional on the named launch spectrum/fraction and loss, not a measured coherent yield. Six 300-μs records provide finite tails, not a fitted universal tail law. Source/component selection, mounted-disk exterior, loaded-tensor applicability and D storage/transfer remain explicit.

The same CLI can consume the saved duration tables without propagating again:

```sh
python g13_rogue/pass2/coarse.py --event-requirements
python g13_rogue/pass2/coarse.py --package
```

`--package` archives only this new folder, excluding the archive itself and external `delivery.json`; its manifest binds each included file. A draft-release byte readback is required separately and conveys no publication or merge permission.

## Corrected live result

[readout.md](readout.md) reports 48 live cases. In the canonical structured-beam case, A/C impact fractions are **93.43% / 1.064%** and off-impact B occupancy **40.28%**; the peaked/steady B case reaches **56%**. The primary A/C states are averaged over 2.5-mm source neighborhoods and frozen before each **50-ns bin**, so within-bin ordering is unresolved. Positive/negative strain and displacement holes, independent hole events, physical class heights, clumping controls and μs-scale plate/edge timing are tabulated separately.

Equal-charge and energy-accounting consumers passed; the impact-height consumer selected the identical parent counts and source work. This is a finite 5–250-kHz **modal-band** model with a declared square/free-edge/mount exterior, not a measurement of acoustic launch yield, loaded elasticity, below-band global plate motion or D-mode conversion. Persistent damage is omitted by the correction, not repurposed as class C.
