# Result-producing commands

Run from the Coherencev1 integration root. These are portable reproduction commands, not historical argv receipts. Existing scientific outputs are immutable; use a fresh output directory. Optional Native/proof sections are not permission to spend a new allocation. No accepted proof was rerun during sealing.

## Pure saved-output pipeline

```bash
D=g4_beam/screen/dynamic
PYTHON=${PYTHON:-python}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p scratch
OUT=$(mktemp -d scratch/g5-dynamic-readout.XXXXXX)

"$PYTHON" -m g4_beam.screen.dynamic.pd18_dataset \
  --directory "$D/data/atomic" >"$OUT/atomic-admission.json"

"$PYTHON" -m g4_beam.screen.dynamic.dynamic_results \
  --moments "$D/data/moments/result.extract.json" >"$OUT/speed10.json"

"$PYTHON" -m g4_beam.screen.dynamic.static_results \
  --registry "$D/data/static/registry.json" \
  --output "$OUT/static-profiles.json" --csv "$OUT/static-profiles.csv"

"$PYTHON" -m g4_beam.screen.dynamic.beam_screen_results \
  --moments "$D/data/moments/result.extract.json" --output-dir "$OUT/beam"

"$PYTHON" -m g4_beam.screen.dynamic.fit_reference_results \
  --cause-fits "$D/data/static/benchmarks/fits_transported.json" \
  --beam-result "$OUT/beam/comparison.json" --speed-result "$OUT/speed10.json" \
  --bulk-screening "$D/data/static/benchmarks/bulk_screening.csv" \
  --output "$OUT/fit-reference.json" --csv "$OUT/fit-reference.csv"

"$PYTHON" -m g4_beam.screen.dynamic.screening_report_results \
  --static "$OUT/static-profiles.json" --fit-reference "$OUT/fit-reference.json" \
  --speed-result "$OUT/speed10.json" \
  --output "$OUT/screening-report.json" --csv "$OUT/screening-report.csv"
```

## Native reproduction — separate allocation required

Requires compatible recorded GPAW/ASE dependencies, source/setup authentication and an explicitly bounded job allocation. Preserve the original flat Native helper import convention. The supplied specifications preserve the actual physical configurations; their paths and source bindings identify the bundled copies. The old two-point wrapper is replaced by the qualified general family driver, not a second physics implementation.

```bash
D=g4_beam/screen/dynamic
export PYTHONPATH="$PWD:$PWD/$D:$PWD/$D/data/static/native-source"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# Retained-wavefunction velocity moments; no SCF.
"$PYTHON" -m g4_beam.screen.dynamic.postprocess_tenvalence \
  --input "$D/data/moments/reproduction-input.json" --output "$OUT/moments"

# Exact admitted atomic recipe, including the high-energy sampled diagnostic.
"$PYTHON" -m g4_beam.screen.dynamic.generate_pd18 \
  --output "$OUT/generated-pd18" --radius-Bohr 1.6 --scan-max-Hartree 12 \
  --projectors '4s,5s,10.0s,4p,5p,10.0p,4d,d,10.0d,0.0f,5.0f,10.0f,G'

# One command per supplied actual response-bank configuration.
for SPEC in "$D"/data/native-specs/*.json; do
  NAME=$(basename "$SPEC" .json)
  "$PYTHON" "$D/run_static_family.py" --input "$SPEC" --output "$OUT/native-$NAME"
done
```

## Optional source-bound verification reproduction

Original accepted numerical results and original source identities are already bundled as scientific JSON. These commands would be NEW executions of identified import-only adapters, not reruns attributed to historical source bytes. NumPy/SciPy, SymPy/mpmath and pytest are required as appropriate. Quadrature qualification additionally imports GPAW but performs no SCF or response calculation.

Supply the exact original PDFs as `YAN_PDF` and `HUKE_PDF`:
- Yan: https://arxiv.org/pdf/1104.1273v1 —443238B, SHA256 `44df983b77c96798d9527af1f84163acd1015f780898f1f255c0b4b2d590a6f8`.
- Huke: https://arxiv.org/pdf/0805.4538v1 —653081B, SHA256 `7a428728bc7d86c12bc8aa9bb659c57bb74d9070ec87d91ba0dbb0aea74cd33b`.

The exact consumed Yan text, Huke rendered page and GPAW Drude source are bundled; the hard source-binding checks remain intact.

```bash
"$PYTHON" -m g4_beam.screen.dynamic.proof.static.verify \
  --yan-text "$D/data/primary/yan-paw-response.txt" \
  --yan-pdf "${YAN_PDF:?Supply the exact original Yan v1 PDF}" \
  >"$OUT/static-proof.json"

"$PYTHON" -m g4_beam.screen.dynamic.proof.static.qualify_quadrature \
  --bank "$D/data/static/banks/qrepair-20261118t055431z" \
  --setups "$D/data/atomic" --output "$OUT/quadrature-proof.json"

"$PYTHON" -m g4_beam.screen.dynamic.proof.verify \
  --huke-pdf "${HUKE_PDF:?Supply the exact original Huke v1 PDF}" \
  --huke-page "$D/data/speed-proof-sources/huke-page10.png" \
  --drude-source "$D/data/speed-proof-sources/chi0_drude.py" \
  >"$OUT/speed-proof.json"

"$PYTHON" -m pytest -q "$D/proof/speed_edges.py"
```
