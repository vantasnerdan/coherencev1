# Electronic heating changes both folded PdD stiffness and coupling

In ideal rocksalt PdD, heating the electrons raises both the diagonal stiffness of the four-atom projected model and the coupling between its two deuterium sublattices. The latter subtracts from the diagonal change for collective motion. The diagonal includes interactions with same-phase periodic images: it is not the isolated on-site stiffness. These are fixed-volume, fixed-Pd electronic free-energy curvatures, not lattice-temperature measurements.

| Electronic T (K) | ΔU (eV/D site) | ΔKΓ (eV/Å²) | X[100] Δfolded diagonal | X[100] Δcross coupling | X[010] Δfolded diagonal | X[010] Δcross coupling |
|---:|---:|---:|---:|---:|---:|---:|
| 3,000 | 0.073680 | +0.07645 | +0.12158 | +0.04512 | +0.12466 | +0.04821 |
| 10,000 | 0.830291 | +0.40254 | +0.66406 | +0.26152 | +0.62747 | +0.22493 |

The table uses the finer equivalent 16³ primitive mesh; all contrasts are relative to that mesh’s own 300 K row. At 3,000 K, the 12³ mesh gives Δcross coupling +0.03789 eV/Å² longitudinally and +0.08381 transversely; at 10,000 K, +0.25387 and +0.26024. The positive changes agree in sign, but their polarization difference is not mesh-stable. Do not fold them into a scalar change. Absolute curvatures, crossed force projections, matched excitation and both mesh differences are in `folded-response.json` and `thermal-changes.json`. All six folded Γ points reproduce their accepted primitive counterparts to 0.000185 eV/Å²; the largest force/energy curvature difference is 0.000991 eV/Å². The transverse total-force residual reaches 6.30×10⁻⁵ eV/Å and is retained, not subtracted. The two meshes are sensitivities, not continuum error bounds.

Under the consumer’s scalar five-shell representation only, the Γ/X difference measures the combined change in shells 1 and 3. It does not identify either individually or disentangle the same-phase shells from the bare on-site term. A redistribution between shells 1 and 3 can leave both observed modes unchanged while changing the small-q curvature. The modes therefore do not determine a unique wave speed. Use the projected contrast directly or treat any assignment to individual shells as a separate model assumption.

The two new primitive temperatures give ΔKΓ +0.00142/+0.00516 eV/Å² at 1,500 K and +0.01599/+0.02234 at 2,000 K (12³/16³), all relative to 300 K. Both 1,000 K samples are lower than the 300 K reference; both 1,500 K samples are higher. **1,000–1,500 K is the provisional sampled return bracket**, not a converged root or the temperature of minimum stiffness: neither endpoint’s smaller shift clears the full intermesh spread. Stiffening is numerically resolved at 2,000 K in this two-mesh comparison. `crossover.json` retains the individual ΔU, signs and mesh differences. Nothing here resolves the previous sign disagreement at the much lower modeled track energy.

## Replay from this directory

```sh
sha256sum -c SHA256SUMS
mkdir replay
python read_folded.py --native-root native-fold --item1 sources/item1-response.json --output replay/folded-response.json
python dispersion_changes.py replay/folded-response.json replay/thermal-changes.json
python verify_folded.py native-fold replay/folded-response.json replay/thermal-changes.json
python read_crossover.py --native-root native-primitive --item1 sources/item1-response.json --output replay/crossover.json
for f in folded-response.json thermal-changes.json crossover.json; do cmp "$f" "replay/$f"; done
```

`DISPERSION.md` gives the registered geometry and reference. `inputs_ledger.md` contains proposed field-consumer rows. The stored native coordinates, forces, electronic occupation metadata and energies reproduce the results without rerunning the electronic calculations. No vacancy, partial-loading or nonthermal calculations were made in this follow-up.
