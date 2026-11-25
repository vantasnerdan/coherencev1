# Lead's note on the PdD refit through 0.48 A

Superseded by the refit through 0.8 A: see [bulk-refit-080](../bulk-refit-080/LEAD.md).

The well of a deuteron in PdD was refitted with first-principles points out to 0.48 A from the site (twice the earlier reach), with no quartic form assumed inside that range. The deuteron comes out slightly more spread than before: mean-square displacement 0.0806 to 0.0828 A^2 at 300 K (plus 2.8 percent) and 0.0635 to 0.0641 A^2 at 77 K (plus 0.9 percent). The share of its probability beyond 0.48 A is 3.3 percent at 300 K. The earlier quartic well was close; the picture in the parent page stands.

What is still extrapolated: everything beyond 0.48 A, and directions between the three computed rays. The two continuations tried differ by more than the numerical error. Points to about 0.8 A, NiH, the pair cut and the surface sites are still to come.

Archive SHA-256 52c9b8fc0973d274244949d00ac1eadddac9997b153f8359cfaf023e8e48d9aa; all file sums check. The lead reran `solve.py --npar 10 --order 48 --tag main` (NumPy and SciPy only): the mean-square displacements and kinetic energies agree with the delivered ones to 5 parts in 1e10. The coupled reweighting was not rerun.

Not committed, for size: the spectrum and full three-dimensional distribution arrays (16 `.npz` files over 1 MB, 84 MB together). `solve.py` regenerates them. One interpreter path was replaced by `python` in `replay.md`, `coupled_sensitivity.py`, `coupled_sensitivity.json` and `results/report_table.json`, so those four files and the 16 omitted ones do not match `SHA256SUMS`.
