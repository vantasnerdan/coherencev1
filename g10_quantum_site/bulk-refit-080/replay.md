# Replay and exact scope

This is one unpublished local G10 bulk update. No push, release/publication mutation, pair/contact or G11 work is authorized or included. Governing plan `coherencev1-g10-bulk080-20261125`, template revision `14c5716f18cd5b1b8ee77624079ed31a0879c64388a625f4db88162dda02faf4`. Accepted .48 evidence comes from `ed06115`; workspace source base is `7349c26e626d41262f9ad66ab94e9e266373425c`.

## Exact input

Giuliano's authenticated Fleet `entry_3684`, GitHub asset **789418823**, supplies the complete native/signed bulk input: archive **94,239 bytes**, SHA-256 `8085c5623dfbbf137a097249e6637c42182565dcb3760ee544e36df1f90553f4`; signed JSON SHA-256 `18fa6232525243391b84f2fb495cfe9947f7d41eef270931fb3dc8b99857e28d`. All 20 extracted members had matched provider hashes before fitting. Exact extracted members, not a redundant nested archive, are included at `g10_quantum_site/.input-080/g10_quantum_site/bulk-wide-wells/`. Its README is authoritative for sign, original atom order and native/derived semantics. The provider references identify adopted zeros/checkpoints; checkpoints are not redistributed.

Each row moves only original atom 25 in the frozen 64-ion cell. `native_force_parallel_eV_A` is along the **fixed signed axis**, so outward slope is −sign(displacement) times that force. Energy is the electronic `delta_energy_free_eV`, not a nuclear free energy. Full forces remain in original atom order. Separately native [110] signs are averaged; inversion-derived negatives carry their parent label and count as zero additional calculations. Old near PdD [100]/[111] E/F and near K/F0 are retained from the pinned accepted input. No absolute electronic energy is compared between materials.

## Prerequisites without duplicating old data

The repository alone is **not sufficient for full report replay**. Main's accepted-tree readback at `4f1fd691af7d797b4bea7dfedc8103154d0b36cf` matched 26 of the 33 exact consumed prerequisite records. Six accepted .48 main distribution arrays were deliberately omitted from that tree; acquire those bytes from the **already retrieved** prior capsule, not by rerunning old spectra. The seventh byte difference is the reporting-only portable script edit described below. The scientific result and existing numerical validation are unchanged.

Prepare a **separate replay tree** from accepted current main, retaining its portable scripts, and extract this new .80 capsule there. Do not overlay the authoritative current checkout or the read-only prior delivery. `prerequisites.json` records each consumed source/result/path by repository-relative name, bytes and SHA-256, plus the accepted script equivalent and exact prior-archive acquisition map. No old arrays are copied into this new archive.

### Six-array, no-regeneration overlay

Use prior asset **789149873**, `g10-pdd-refit-048.tar.gz`, SHA-256 **`52c9b8fc0973d274244949d00ac1eadddac9997b153f8359cfaf023e8e48d9aa`**, from the [authenticated asset endpoint](https://api.github.com/repos/vantasnerdan/coherencev1/releases/assets/789149873). The prior capsule was already retrieved; point the following portable `python` command at that local file and the separate replay tree. It extracts **only**:

```text
g10_quantum_site/bulk-refit-048/results/D_77_old_main_distributions.npz
g10_quantum_site/bulk-refit-048/results/D_77_quadratic_main_distributions.npz
g10_quantum_site/bulk-refit-048/results/D_77_quartic_main_distributions.npz
g10_quantum_site/bulk-refit-048/results/D_300_old_main_distributions.npz
g10_quantum_site/bulk-refit-048/results/D_300_quadratic_main_distributions.npz
g10_quantum_site/bulk-refit-048/results/D_300_quartic_main_distributions.npz
```

```sh
python - /path/to/g10-pdd-refit-048.tar.gz /path/to/separate-replay-tree <<'PY'
import hashlib
import json
from pathlib import Path
import sys
import tarfile

archive = Path(sys.argv[1]).expanduser()
root = Path(sys.argv[2]).expanduser().resolve()
manifest = json.loads((root / "g10_quantum_site/bulk-refit-080/prerequisites.json").read_text())
overlay = manifest["prior_archive_overlay"]
with archive.open("rb") as source:
    assert hashlib.file_digest(source, "sha256").hexdigest() == overlay["archive_sha256"]
wanted = set(overlay["members"])
expected = {entry["path"]: entry for entry in manifest["external_prerequisites"]}
seen = set()
with tarfile.open(archive, "r:gz") as source:
    for member in source:
        if member.name not in wanted:
            continue
        assert member.isfile() and member.name not in seen
        stream = source.extractfile(member)
        assert stream is not None
        with stream:
            raw = stream.read()
        record = expected[member.name]
        assert len(raw) == record["bytes"]
        assert hashlib.sha256(raw).hexdigest() == record["sha256"]
        target = root / member.name
        assert target.resolve().is_relative_to(root)
        if target.exists():
            with target.open("rb") as existing:
                assert hashlib.file_digest(existing, "sha256").hexdigest() == record["sha256"]
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as output:
                output.write(raw)
        seen.add(member.name)
assert seen == wanted
print("Six exact prior arrays available; no script overwritten and no numerical regeneration.")
PY
```

### Consumed and accepted portable script provenance

Keep current main's **portable** `g10_quantum_site/bulk-refit-048/coupled_sensitivity.py`. Do **not** extract the prior script over it. The consumed pre-portability script was **16,751 bytes**, SHA-256 `5a7b84c13be1171adf5603a64012754bafbef1897744c4ae20849a42a79c00ee`; those exact bytes remain available as that member in prior asset 789149873 if strict source-epoch reproduction is needed in another isolated tree.

Main observed exactly one reporting-string substitution, `(an interpreter path under the author home directory)` → `python`, with callbacks/operator unchanged. Applying that declared single byte substitution to the consumed source gives the accepted portable identity: **16,709 bytes**, SHA-256 `10988fafc1061660f234101526288d885d5705af3bcc66caf1f40b9f45ddb476`. Both identities and the default **keep accepted portable script** selection are explicit in `prerequisites.json`. The consumed hash is provenance, not a demand to overwrite the accepted portability edit. No numerical recalculation or old/current mutation is required.

Large generated arrays (>1 MB) are omitted from this small capsule and listed separately in `prerequisites.json` as **regenerable generated outputs**; run the relevant command below if you need eigenvectors or arbitrary-vector densities. They are not required to read the supplied common tables, directional curves, figures or numerical comparisons. No source implementation is omitted by this rule.

## Machine-independent commands

Python 3.12 with NumPy, SciPy and Matplotlib is sufficient. Select your local environment at invocation; replay code and commands use portable `python`, not an author's home path. No network or installation action is performed by the scripts.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
D=g10_quantum_site/bulk-refit-080
python -B "$D/solve.py" --variant wide_quadratic
python -B "$D/solve.py" --variant wide_quadratic --check
python -B "$D/solve.py" --variant wide_quartic
python -B "$D/solve.py" --variant wide_quartic --check
python -B "$D/coupled.py"
# NiH replay commands and chosen angular models: nih/README.md
python -B "$D/deliver.py"
python -B "$D/verify.py"
python -B "$D/package.py"
```

At most two 1-thread heavy fits concurrently. The .80 PdD replay imports the exact accepted `bulk-refit-048/solve.py`; it does not copy the spectral/Fourier implementation. `well.py` rebinds only the radial callback in the accepted analytic cubic angular operator. Main uses 20 HO levels per axis, potential quadrature 48, position pitch .020 Å, momentum pitch .60 Å⁻¹ and angular GL24×48. The bounded comparison uses 24 levels, quadrature 64, pitches .016 Å/.48 Å⁻¹ and GL32×64. Old/main/check results are read, never solved again. Temperature changes canonical weights, not electronic surface or geometry.

## Equations and interpretation

Isolated operator: H=−ℏ²∇²/(2m)+V(x), positions in Å and energies in eV. Genuine canonical mixtures of position/Fourier eigenfunctions produce density and radial measure; momentum is k=p/ℏ in Å⁻¹. All parity sectors are retained through representatives with multiplicities 1/3/3/1. Equal coordinate variances follow cubic symmetry; angular line anisotropy does not break it. The signed [111] table is ∫ρ(sn+a u+b v) da db, where n=+[111]/√3 and u,v span its transverse plane. It is not ρ(sn) and does not specify an impact trajectory or a slab.

For R=.80, t=max(r−R,0), primary exterior is E_R+s_R t+K t²/2. Alternative is E_R+s_R t+(K+3a_ray R²)t²/2+a_ray R t³+a_ray t⁴/4. Both confine and coincide through R. Hermite interpolation exactly matches E and slope at all retained ray knots; analytical segment derivative extrema are checked for monotonicity. Cubic nonnegative weights w=(1−4q+9p,4q−36p,27p), q=sum_{i<j} z_i z_j, p=z_1z_2z_3, z_i=x_i²/r², reuse the accepted angular closure. The exact old quartic is used at r≤.03 Å. The analytic force includes angular derivatives. This preserves near K, F0 and angular closure but cannot validate unmeasured off-ray physics.

Tagged finite-P consumer changes **one** on-site term in the old occupied host, with ratio weight exp[−β/P Σ_b(Vnew−Vold)]. It uses all archived beads and their actual neighbour field, not a mean field or independently redrawn ensemble. Nine contiguous 89-frame blocks yield ratio influence SE and paired old/.48/.80 differences; ESS is sampled participation, not autocorrelation correction or an unseen-overlap bound. At P=32/128, no continuum-P extrapolation is claimed. The signed [111] displaced-path characteristic uses the already accepted callback at separation 0/.25/.50 Å (conjugate to k), four cuts and antithetic signs. Three characteristic samples do not determine a momentum density or kinetic energy. No Fourier inversion is attempted for coupled paths.

`results/observables.json` carries common observable fields; support/model/numerical metadata is separate. Null means unsupported, never a fitted or truncated zero. Empty sampled >.80 bins are unobserved rare tails. Isolated tail comparisons are empirical resolution checks on the declared finite model, not rigorous tail enclosures, DFT errors or angular uncertainty estimates. The paired exterior band is a chosen model sensitivity, not a probability/confidence interval and not guaranteed to bracket later data.

`verify.py` is one settled changed-consumer check: source identity, joined table readback, distribution normalizations, signed-marginal meaning/symmetry, anchored forces and imported coupled witness evidence. It does not repeat the spectra, NiH fitting, old sampling, or independent physical validation. The parent owns final archive readback and recipient delivery.
