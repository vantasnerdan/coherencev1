# Compact Pd/PdO saved-output delivery

**PARTLY RELAXED GEOMETRY: authority-labelled/conditional calibration only, NOT force-qualified.**

**GEOMETRY-CONVERGENCE ERROR IS UNBOUNDED: geometry_error_bound_eV=null. Residual forces and h/k tests do not bound pair-lowering/contrast error.**

Paired h/k endpoint criterion: NOT QUALIFIED. Geometry force criterion: NOT FORCE-QUALIFIED. Missing case files: plain-h, oxygen-h, plain-k, oxygen-k.

| Geometry | Status | Actual free fmax (eV/A) | Actual raw fmax (eV/A) |
| --- | --- | ---: | ---: |
| plain | force-qualified | 0.036036094705237193 | 0.22639077804335472 |
| oxygen | partly-relaxed | 0.42174246367039653 | 0.42174246367039653 |

| BO role | Ray status | Final coated admitted |
| --- | --- | --- |
| oxygen | precommitted shared-open | true |
| plain | precommitted ray; final coated shared-open admission pending | false |

| BO role | Actual completed R (A) | Missing commissioned R (A) |
| --- | --- | --- |
| oxygen | [0.1, 0.2, 0.4, 0.8, 1.2, 1.6] | [] |
| plain | [0.1, 0.2, 0.4, 0.8, 1.2, 1.6] | [] |

Actual common positive-R support: [0.1, 0.2, 0.4, 0.8, 1.2, 1.6] A. Only these exact measured radii are contrasted; no interpolation or invented far points. R0.1 is selected by value, never by execution ordinal.

The `*_BO_pair_common.json` files retain only that measured overlap; the ordinary per-variant curves retain every actual completed point. `common-support.json` records the paired contrast and endpoint support. Signed contact remains a separate sidecar, not a fabricated positive-R point.

Recorded BO stop(s): 2026-11-21T03:00:00Z; not an authorization assertion. `--allow-partial` permits missing case files only, never missing geometry/cap/site evidence or false force labels.

Bare 15Pd+D versus coated 19Pd+4O+D in the same fixed cell; the added sheet is 4Pd+4O, not oxygen insertion at fixed total Pd. Only bottom substrate indices 10–14 are fixed. D and every other atom are free. The actual relaxed target position, migration and raw fixed-atom forces are reported, not an inferred retained-octahedral site.

Entry_2980/2904 retains strict <=0.30 eV/A partial eligibility. The distinct entry_2994/3001 route permits the specifically authorized ORIGINAL step5 after the actual startup-failed resume: **partly relaxed, fmax 0.42 eV/A, no geometry error bound**. Its original Native state, specification, accepted step and process clock remain unchanged. The actual failed attempt is separate closed evidence and cost: zero new resumed SCFs/steps, no qualified checkpoint import or geometry advance. A startup failure is not misrepresented as a new cap. Any second-frame spread is sensitivity only, never a rigorous geometry-error bound.

Commissioned six base radii: 1.6, 1.2, 0.8, 0.4, 0.2, 0.1 A; only the actual completed sets above are delivered. Near branches execute 0.1, 0.2, 0.4 first; optional far branches start independent cold sequences and reuse original fixed host/contact records, not recomputed hosts. Parallel container metadata is a checked derived view, not an executed monolithic specification; each Raw state retains its actual branch input, original bytes and execution/sequence indices. Paired Gamma versus 2x2x1 and h0.20 versus h0.16 checks apply only at R0.1/contact. Each contrast-change tolerance is max(5 eV, 10% of the absolute base coated-minus-bare contrast). Endpoint h/k qualification is separate from geometry force qualification; no full scientific qualification is claimed. Slab-thickness, cell/vacuum and full-radius h/k controls were omitted; no corresponding uncertainty bound or isolated-limit claim.

Under entry_2971 the immutable precommit selects on force-qualified bare plus the actual accepted coated boundary: vacuum-side normal first if its paired clearance is >=1 A, otherwise the most-open common declared candidate on those CURRENT frames. The coated boundary is a geometric witness only: not final, not force-qualified, not partly-relaxed BO admission and not BO-eligible. Bare data retain final-coated admission pending. Actual final coated data require strict force-qualified, explicit entry_2904 partial eligibility, or distinct exact entry_2994/3001 labelled-original admission and continuous clearance >1e-6 A along the EXACT same precommitted direction; their label is **precommitted shared-open**, never maximum-common-final. The 1 A threshold is normal preference only, not a final-admission cutoff. A blocked final ray is refused without changing the bare ray. No global or final-frame optimum is claimed. Geometry-only records have no selected common ray. Positive-R data are not spliced to signed regular contact. No Ueq, action/Eq7, or experimental phase claim.

Private inputs were checked by reader_pd.replay with the accepted Native-state, epoch, contact-operator, Green and energy helpers. This small package contains no GPWs, raw logs or original private source trees. Provider basenames and byte/SHA bindings are retained. Run `python3 export_pd.py verify DIRECTORY` for archive hashes, numerical replay from the saved energy ledger, and comparisons. That small-package check is not a fresh raw-log/checkpoint/operator-source/cap audit; use reader_pd.py on the original private case roots for that.
