"""Reproduce G11's intrinsic NN witness and explicitly normalized DD model.

Run from pinned repository root: python -m g11_cmns.coupling.replay --nuhamil PATH
"""
import argparse
import ctypes as ct
from ctypes.util import find_library
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.special import loggamma, roots_legendre, spherical_jn

from tools_marvin.boost import projected_boost


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NUHAMIL_PIN = "eb65165f5ecf11fde6e188a65385c98e7526ccb6"
NN_FILE = "N3LO_EM500_kmax8_N100_Jmax8.bin"
NN_SHA256 = "0f3af5611bdc0d187120ebf3e2cb1312d094b6b738c5a1593a4914344ac0f64a"
M_N = 938.92  # MeV
M_4 = 4 * M_N
M_D = 1875.61294257  # MeV
HBARC = 197.3269804  # MeV fm
ALPHA = 1 / 137.035999177
Q_DD = 23.84653  # MeV, threshold relative to alpha ground state
CASES = ((4, True), (6, True), (6, False))


def sha(path):
    with path.open("rb") as fp:
        return hashlib.file_digest(fp, "sha256").hexdigest()


def build(src, dst):
    objects = sorted(p for p in (src / "obj").glob("*.o") if p.name != "NuHamilMain.o")
    if not objects or not (src / "mod" / "abodymanager.mod").is_file():
        raise ValueError("Build the pinned NuHamil checkout as in tools_marvin/README.md first")
    hdf5 = subprocess.check_output(["pkg-config", "--libs-only-L", "hdf5"], text=True).split()
    subprocess.run(["gfortran", "-O2", "-fopenmp", "-ffree-line-length-none", "-I", src / "mod",
                    HERE / "intrinsic.F90", *objects, *hdf5, "-lgsl", "-lz", "-lhdf5_fortran",
                    "-lhdf5", "-lm", "-ldl", "-lblas", "-llapack", "-o", dst], check=True)


def as_pair(z):
    return [float(z.real), float(z.imag)]


def intrinsic_case(exe, folder, nn, nmax, coulomb, momentum):
    folder.mkdir()
    path = folder / "input.dat"
    path.write_text(f'''&input
rank=4, hw=28, renorm='bare', Operators='hamil'
NNInt='N3LO_EM500', input_nn_file='{nn}'
bra=0,1,0,0, ket=0,1,0,0, trans2lab=.false.
NN_only=.true., genuine_3bf=.false., coul=.{str(coulomb).lower()}., pn_same_mass=.true.
N2max={nmax}, NAmax={nmax}, N3max={nmax}, jmax3={2*nmax+1}, ramp='flat{nmax}'
/
''')
    proc = subprocess.run([exe, path.name], cwd=folder, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    if proc.returncode or not (folder / "g11_matrices.bin").is_file():
        raise RuntimeError(proc.stdout[-5000:])
    with (folder / "g11_matrices.bin").open("rb") as fp:
        n, m = map(int, np.fromfile(fp, dtype="=i4", count=2))
        mu, pair, pair_ref, herm, closure, closure_ref = np.fromfile(fp, dtype="=f8", count=6)
        e, o = (np.fromfile(fp, dtype="=f8", count=k) for k in (n, m))

        def fmatrix(a, b):
            return np.fromfile(fp, dtype="=f8", count=a*b).reshape((a, b), order="F")

        even, odd, up, down = fmatrix(n, n), fmatrix(m, m), fmatrix(m, n), fmatrix(n, m)
        assert not fp.read(1)
    anchors_e, anchors_o = np.argmax(abs(even), axis=0), np.argmax(abs(odd), axis=0)
    even *= np.sign(even[anchors_e, np.arange(n)])
    odd *= np.sign(odd[anchors_o, np.arange(m)])
    residuals, orthogonality = [], []
    for v, energies, name in ((even, e, "even"), (odd, o, "odd")):
        h = np.fromfile(folder / f"Hamil_{name}.op", dtype="=f8").reshape(v.shape, order="F")
        residuals.append(float(np.max(np.linalg.norm(h @ v - v*energies, axis=0))))
        orthogonality.append(float(np.max(abs(v.T @ v - np.eye(v.shape[0])))))
    assert max(residuals) < 1e-9 and max(orthogonality) < 1e-10
    # Actual nuclear eigenbasis, not HO basis indices; direct sum ordered [even, odd].
    tensor = np.zeros((3, n+m, n+m))
    tensor[2, n:, :n] = odd.T @ up @ even
    tensor[2, :n, n:] = even.T @ down @ odd
    energies = np.r_[e, o]
    C = projected_boost(energies, tensor, [0, 0, 1/M_4], M_N)  # Pc = 1 MeV; result divided by 1 MeV
    assert np.count_nonzero(C[:n, :n]) == 0 and np.count_nonzero(C[n:, n:]) == 0
    assert np.max(abs(C-C.T.conj())) < 1e-14 and abs(tensor[2, n, 0]) > 1.
    pc0 = momentum["optical_37meV_two_mass_reference"]["one_cell_one_moving_D_abs_Pc_eV"]
    pc_driven = momentum["pumped_PdD_archived"]["spatial_rms_of_site_peak_Pc_eV_from_mass_times_velocity"]
    return {
        "Nmax": nmax, "hw_MeV": 28, "Coulomb_fixed_T": coulomb,
        "dimensions_even_odd": [n, m], "phase": "largest absolute native eigenvector component chosen positive (zero-based anchor row); cross-truncation signs need not be gauge aligned",
        "E_absolute_even_MeV_first_3": e[:3].tolist(), "E_absolute_odd_MeV_first_3": o[:3].tolist(),
        "phase_anchor_basis_rows_even_first_3": anchors_e[:3].tolist(), "phase_anchor_basis_rows_odd_first_3": anchors_o[:3].tolist(),
        "O_z_odd_bra_even_ket_MeV_first_3_by_3": tensor[2, n:n+3, :3].tolist(),
        "C_z_odd_bra_even_ket_eV_per_eV_first_3_by_3": [[as_pair(z) for z in row] for row in C[n:n+3, :3]],
        "same_parity_blocks_nonzero_elements": int(np.count_nonzero(C[:n, :n]) + np.count_nonzero(C[n:, n:])),
        "O_z_alpha_to_all_retained_odd_norm_MeV": float(np.linalg.norm(tensor[2, n:, 0])),
        "H_1phonon_odd0_alpha_one_cell_one_D_eV": as_pair(C[n, 0]*1j*pc0),
        "H_drive_peak_abs_odd0_alpha_one_moving_D_eV": float(abs(C[n, 0])*pc_driven),
        "max_eigen_residual_MeV_even_odd": residuals, "max_orthogonality_error_even_odd": orthogonality,
        "O_z_reverse_transpose_max_error_MeV": herm,
        "m_pair_reduced_c2_MeV": mu,
        "stretched_pair_O_z_MeV": pair, "stretched_pair_analytic_MeV": pair_ref,
        "intrinsic_0s4_O_z_squared_norm_MeV2": closure, "intrinsic_0s4_analytic_MeV2": closure_ref,
    }


class GslResult(ct.Structure):
    _fields_ = [("value", ct.c_double), ("absolute_error", ct.c_double)]


def coulomb_regular():
    gsl = ct.CDLL(find_library("gsl"))
    gsl.gsl_set_error_handler_off.argtypes = []
    gsl.gsl_set_error_handler_off.restype = ct.c_void_p
    gsl.gsl_set_error_handler_off()
    fg = gsl.gsl_sf_coulomb_wave_FG_e
    fg.argtypes = [ct.c_double]*3 + [ct.c_int] + [ct.POINTER(GslResult)]*4 + [ct.POINTER(ct.c_double)]*2
    fg.restype = ct.c_int

    def regular(eta, x, ell):
        f, df, g, dg = GslResult(), GslResult(), GslResult(), GslResult()
        ef, eg = ct.c_double(), ct.c_double()
        status = fg(eta, x, ell, 0, f, df, g, dg, ef, eg)
        if status:
            raise ArithmeticError(f"GSL Coulomb F status {status} for eta={eta}, kr={x}, l={ell}")
        return float(f.value*np.exp(ef.value))
    return regular


def coulomb_origin_series(eta, x, ell):
    """Independent small-kr Frobenius control for the charged regular wave."""
    logc = ell*np.log(2.)-np.pi*eta/2+loggamma(ell+1.+1j*eta).real-loggamma(2*ell+2.)
    prev, term, total = 0., 1., 1.
    for n in range(1, 24):
        prev, term = term, (2*eta*x*term - x*x*prev)/(n*(n+2*ell+1))
        total += term
    return np.exp(logc)*x**(ell+1)*total


def cluster_sensitivity(momentum):
    """Frozen spin-1 DD clusters, Gaussian 0+ surrogate and exact repulsive Coulomb 1- waves."""
    F = coulomb_regular()
    assert abs(F(0., .73, 1)-.73*spherical_jn(1, .73)) < 1e-13
    mu = M_D/2
    radii = (1.2, 2.0)
    bands = ((237.5, 262.5), (2375., 2625.))  # eV; one normalized top-hat packet per band
    pc0 = momentum["optical_37meV_two_mass_reference"]["one_cell_one_moving_D_abs_Pc_eV"]
    pc_drive = momentum["pumped_PdD_archived"]["spatial_rms_of_site_peak_Pc_eV_from_mass_times_velocity"]

    op = np.zeros((3, 2, 2))
    op[2, 1, 0] = op[2, 0, 1] = 1.
    def integrate(nr, ne):
        x, w = roots_legendre(nr)
        r, wr = 12*(x+1), 12*w  # radial cutoff 24 fm, 12 times largest b
        # Reduced radial u=rR is normalized in dr, and (d/dr-1/r)u= -r*u/b².
        u = np.stack([2/np.pi**.25/b**1.5 * r*np.exp(-r*r/(2*b*b)) for b in radii])
        du = -u*r/np.square(radii)[:, None]
        assert np.max(abs((u*u)@wr-1)) < 1e-12
        # Rank-one separable s-wave potential -|g><g|/<u|g>, g=(H_C+Q)u,
        # makes the prescribed u an eigenfunction at -Q in the same cluster Hamiltonian.
        # It has no p-wave action; changing b is model sensitivity, not a wave-function fit.
        b = np.asarray(radii)[:, None]
        hu = (-HBARC**2/(2*mu)*(r*r/b**4-3/b**2)+ALPHA*HBARC/r)*u
        separable_denominator = ((hu+Q_DD*u)*u)@wr
        assert min(separable_denominator) > 0

        def at(eV):
            E = eV*1e-6
            k = np.sqrt(2*mu*E)/HBARC
            eta = ALPHA*mu/(HBARC*k)
            sigma = loggamma(2.+1j*eta).imag
            # u from <E|E'>=delta(E-E') with E in eV, not in MeV; factor 1e-3.
            scale = 1e-3*np.sqrt(2*mu/(np.pi*HBARC**2*k))
            u_p = scale*np.array([F(eta, k*rr, 1) for rr in r])
            me = 4/3*HBARC*(du @ (u_p*wr))
            # alpha bra, p-wave ket. Apply Eq.27 at each continuum energy, then combine amplitudes.
            at7 = [F(eta, k*7, l) for l in (0, 1)]
            error = max(abs(f/coulomb_origin_series(eta, k*7, l)-1) for l, f in enumerate(at7))
            assert error < 1e-7, error
            C_per_E = projected_boost([-Q_DD, E], op, [0, 0, 1/M_4], M_N)[0, 1]*me
            assert np.all(me < 0) and np.all(C_per_E.imag > 0)
            return C_per_E, complex(np.exp(1j*sigma)), me, float(eta), [scale*f for f in at7], error

        en, ew = roots_legendre(ne)
        out = []
        for low, high in bands:
            width = high-low
            mid = (low+high)/2
            c, phase, me, eta, entrance, origin_error = at(mid)
            weights = ew*np.sqrt(width)/2  # multiply dE, then divide by sqrt(width)
            values = [at(mid+x*width/2) for x in en]
            standing = weights @ np.stack([v[0] for v in values])
            incoming = weights @ np.stack([v[0]*v[1] for v in values])
            for j, radius in enumerate(radii):
                out.append({
                    "gaussian_b_fm": radius, "mid_Ecm_eV": mid, "packet_Ecm_eV": [low, high],
                    "Coulomb_eta": eta, "Coulomb_incoming_phase_e_i_sigma1": as_pair(phase),
                    "u_s_u_p_at_7fm_per_sqrt_eV_fm": entrance,
                    "Coulomb_F_relative_difference_from_origin_series_at_7fm": origin_error,
                    "separable_form_factor_denominator_MeV": float(separable_denominator[j]),
                    "O_z_standing_per_energy_MeV_per_sqrt_eV": float(me[j]),
                    "C_standing_per_energy_per_sqrt_eV": as_pair(c[j]),
                    "C_incoming_per_energy_per_sqrt_eV": as_pair(c[j]*phase),
                    "C_standing_packet_eV_per_eV": as_pair(standing[j]),
                    "C_incoming_packet_eV_per_eV": as_pair(incoming[j]),
                    "H_incoming_packet_1phonon_one_cell_one_D_eV": as_pair(incoming[j]*1j*pc0),
                    "H_incoming_packet_drive_peak_abs_one_moving_D_eV": float(abs(incoming[j])*pc_drive),
                })
        return out

    lo, hi = integrate(72, 12), integrate(144, 24)
    for a, b in zip(lo, hi, strict=True):
        key = "C_incoming_packet_eV_per_eV"
        x, y = complex(*a[key]), complex(*b[key])
        b["relative_packet_quadrature_difference"] = abs(x-y)/abs(y)
        assert b["relative_packet_quadrature_difference"] < 1e-5
    return hi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nuhamil", type=Path, required=True, help="Pinned built NuHamil checkout with obj/ and mod/")
    args = parser.parse_args()
    src = args.nuhamil.resolve(strict=True)
    rev = subprocess.check_output(["git", "-C", src, "rev-parse", "HEAD"], text=True).strip()
    if rev != NUHAMIL_PIN:
        raise ValueError(f"NuHamil HEAD {rev} does not match {NUHAMIL_PIN}")
    nn = src / "input_nn_files" / NN_FILE
    if sha(nn) != NN_SHA256:
        raise ValueError("NN input SHA-256 mismatch")
    p = HERE / "collective_input.json"
    contract = json.loads(p.read_text())
    momentum = contract["momentum"]
    with tempfile.TemporaryDirectory(prefix="g11-intrinsic-") as temp:
        work = Path(temp)
        exe = work / "intrinsic"
        build(src, exe)
        rows = [intrinsic_case(exe, work / f"n{n}-c{coul}", nn, n, coul, momentum) for n, coul in CASES]
    print(json.dumps({
        "schema": "coherencev1.g11.coupling.evidence.v1", "NuHamil_revision": rev,
        "interaction_file": NN_FILE, "interaction_sha256": NN_SHA256,
        "collective_input_sha256": sha(p), "boost_source_sha256": sha(ROOT / "tools_marvin/boost.py"),
        "intrinsic_status": "4-nucleon Jacobi finite HO NN eigenvectors, positive-energy Eq.27 witness; no computed microscopic DD continuum or pole normalization",
        "intrinsic_witness": rows,
        "cluster_status": "single-partition pure Coulomb sensitivity, frozen pointlike spin-1 D, Gaussian surrogate alpha; incoming convention adds e^(i sigma_1); packet not a material g or cross section",
        "cluster_DD_s_wave_direct_eq27_C": [0., 0.],
        "cluster_coulomb_dd_sensitivity": cluster_sensitivity(momentum),
    }, indent=2))


if __name__ == "__main__":
    main()
