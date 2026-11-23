# Scientific toolset: lattice, few-body, and open dynamics

Installed and exercised on Linux x86-64, Python 3.12 and GNU Fortran 13. These are small reference cases, not foil predictions.

| Layer | Installed tools and observed case | What the case establishes |
|---|---|---|
| Quantum lattice | [i-PI](https://docs.ipi-code.org/getting-started.html) 3.3.1 + ASE 3.29.0, analytic force driver; one D, 16 beads, 300 K, archived cubic `g2_refit` well | Axis mean-square displacement **0.022155(505) Å²** (block SE), versus **0.022732 Å²** from direct diagonalization of the same well. The observable includes bead spread, not just the centroid. |
| Anharmonic free energy | [python-sscha and CellConstructor](https://sscha.eu/download/) 1.6.2; translationally invariant dimer with the same D relative mass and well; 1,024 configurations per force gradient | SSCHA stochastic force-gradient zero at **1.6039 eV/Å²**, versus analytic Gaussian self-consistency **1.6186 eV/Å²**. This curvature is not a discrete vibrational level. i-PI samples equilibrium; neither case implies real-time nuclear motion. |
| Nuclear structure | [NuHamil](https://arxiv.org/abs/2302.07962) `eb65165`, intrinsic Jacobi no-core shell model; bundled bare N³LO EM500 NN interaction + Coulomb, no 3N or SRG | Ground energies in MeV: D **−2.224348** (`Nmax=100`, 30 MeV) vs **−2.224575** deuteron binding; ³H **−7.695892** and ³He **−6.947934** (`16`, 30 MeV); ⁴He **−0.648207** (`2`, 18 MeV) vs **−0.6527** in [published Fig. 2(e)](https://doi.org/10.1016/j.physletb.2024.138749), a 4.5-keV difference, not exact parity. The ⁴He (`8`, 28 MeV) example gives **−17.414016**; basis and frequency differ, and these are not converged binding predictions. |
| Relativistic coupling | NumPy 2.5.5 reduced positive-energy [Dirac operator, Eq. 27](https://doi.org/10.1088/1361-6455/acf3be); two-state harmonic reference | Signed off-diagonal matrix element **−1.45949×10⁻⁴ MeV**, matching the analytic example. This is a test input, not a material coupling. NuHamil ⁴He export gives a five-component Jacobi state and matching Hamiltonian; $\lVert H\psi-E\psi\rVert=2.6\times10^{-9}$ MeV. |
| Collective dynamics | [QuTiP](https://qutip.readthedocs.io/en/latest/guide/dynamics/dynamics-master.html) 5.3.1 and [quimb](https://quimb.readthedocs.io/) 1.16.0; four driven, exchanging spins with local amplitude damping | Against a dense matrix-exponential reference, maximum density-matrix element difference **8.2×10⁻¹²** (QuTiP) and **1.6×10⁻⁵** (quimb local-gate split, max bond 16). The tensor part applies scheduled local Lindblad and exchange gates to the vectorized density operator; it is not a built-in general MPO Lindblad solver. |

**Selection.** i-PI uses a cheap force provider without a DFT call per bead; SSCHA gives a separate Gaussian free-energy reference. NuHamil includes the interaction and Jacobi bound-state solver, avoiding a separate shell-model interaction conversion or a valence-only space. QuTiP gives an independent dense reference for quimb's programmable local tensor contractions. The sizes above are observed cases, not scaling promises: keeping a dense density matrix costs the square of the Hilbert dimension; tensor evolution depends on operator entanglement and mode connectivity.

## Install and run

From the repository root. Prerequisites on Ubuntu 24.04: `python3-venv python3-dev gfortran libopenblas-dev libgsl-dev libhdf5-dev zlib1g-dev pkg-config make gcc`. Keep the environment outside the checkout; set `TOOL_ENV` to its desired path.

```sh
: "${TOOL_ENV:?set TOOL_ENV to the installation directory}"
python3 -m venv "$TOOL_ENV/venv"
. "$TOOL_ENV/venv/bin/activate"
python -m pip install -r tools_marvin/requirements.txt
python -m pip install --no-build-isolation cellconstructor==1.6.2 python-sscha==1.6.2

git clone --depth 1 https://github.com/Takayuki-Miyagi/NuHamil-public.git "$TOOL_ENV/NuHamil"
ref=eb65165f5ecf11fde6e188a65385c98e7526ccb6
if [ "$(git -C "$TOOL_ENV/NuHamil" rev-parse HEAD)" != "$ref" ]; then
  git -C "$TOOL_ENV/NuHamil" fetch --depth 1 origin "$ref"
  git -C "$TOOL_ENV/NuHamil" switch --detach "$ref"
fi
git -C "$TOOL_ENV/NuHamil" submodule update --init --recursive --depth 1
mkdir -p "$TOOL_ENV/bin"
(
  cd "$TOOL_ENV/NuHamil"
  make dirs INSTLDIR="$TOOL_ENV/bin"
  make 'FFLAGS=-O3 -fopenmp -ffree-line-length-none -DVERSION=\"$(VERSION)\"' \
    "LFLAGS=$(pkg-config --cflags --libs-only-L hdf5) -lgsl -lz -lhdf5_fortran -lhdf5 -lm -ldl -lblas -llapack"
  install -m 0755 exe/NuHamil.exe "$TOOL_ENV/bin/NuHamil.exe"
)

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python tools_marvin/lattice.py ipi
python tools_marvin/lattice.py sscha
python tools_marvin/nuclear.py --exe "$TOOL_ENV/bin/NuHamil.exe" \
  --interaction "$TOOL_ENV/NuHamil/input_nn_files/N3LO_EM500_kmax8_N100_Jmax8.bin"
python tools_marvin/boost.py
python tools_marvin/collective.py
```

`TOOL_ENV` must be an absolute path so the nested build and simulation directories refer to the same installation. The `ipi` distribution is not the older `i-PI` distribution. `six` is explicitly installed for the CellConstructor runtime. The NN input SHA-256 is `0f3af5611bdc0d187120ebf3e2cb1312d094b6b738c5a1593a4914344ac0f64a`. For the demonstrated native eigenpair export, add `--export-dir NEW_DIR he4-n2` to the nuclear command; it creates an NPZ outside the transient solver scratch.

## First proposed project uses

1. **Lattice.** Compare D distribution and anisotropy at the PdD bulk octahedral, adsorbed, first-subsurface, and vacancy-trapped sites. Start with a fitted, held-out-tested local force surface and exact levels on the same input, then use i-PI for equilibrium spread and SSCHA for selected anharmonic mode curvatures. The electronic-structure side supplies relaxed cells/slabs, masses, site definitions, displacements, energies and forces across each well and between sites; validate a learned potential before extrapolating to NiH/TiD or larger cells.
2. **Nuclear.** Use bound-state exports to check spin, parity and intrinsic-momentum conventions, then form one allowed opposite-parity transition with a consistently projected spin–momentum operator and one specified vibration. `projected_boost` consumes `[3,n,n]` intrinsic $\sum_j\boldsymbol\sigma_j\times(\boldsymbol\pi_j c)$ matrix elements in the same eigenbasis, energies, and the lattice $\langle f|\mathbf P/(Mc)|i\rangle$. It computes the term linear in that lattice matrix element, not a decay width. The bound-state export does not supply the spin–momentum transition matrix or a scattering-state normalization; D–D entrance and open final channels need a consistent continuum treatment. Supply mass-normalized mode eigenvectors and a time-dependent nuclear center-of-mass momentum projection, not a scalar density multiplier.
3. **Collective.** Start with explicit nuclear levels, one candidate receiver-mode set, and separate population/dephasing and gamma jumps. Compare energy reaching the modeled modes to energy going into the gamma channel; keep the **24 MeV** redistribution question as an input-dependent observable, not an assumed cascade. Use dense QuTiP for the small graph and the same local gates in quimb when a larger graph has manageable bond dimension. Supply mode energies/eigenvectors, transition amplitudes with phases, linewidths with their definitions, bath spectra/temperatures, and radiative branching. The four-spin case establishes none of those material inputs or a microsecond prediction.
