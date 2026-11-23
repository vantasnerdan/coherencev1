"""NuHamil Jacobi NCSM with the bundled bare Entem–Machleidt NN force."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np


CASES = {
    "d": (2, "1,1,0,0", 30, 100),
    "t": (3, "1,1,1,1", 30, 16),
    "he3": (3, "1,1,1,-1", 30, 16),
    "he4-n2": (4, "0,1,0,0", 18, 2),
    "he4": (4, "0,1,0,0", 28, 8),
}
VALUE = r"(-?\d+\.\d{8})"
ENERGY = {2: rf"^\s*{VALUE}\s*$", 3: rf"^\s*Energy\s*\n\s*{VALUE}\s*$", 4: rf"^\s*Energy:\s*{VALUE}\s*$"}


def run_case(exe, interaction, label, export_dir):
    rank, sector, hw, nmax = CASES[label]
    # NuHamil uses twice J,T,Tz for A>=3 and actual J,T,Tz for A=2.
    with tempfile.TemporaryDirectory(prefix="nuhamil-") as tmp:
        dir_ = Path(tmp)
        if rank == 2:
            truncation = f"N2max={nmax}"
        else:
            truncation = f'NAmax={nmax}, N3max={nmax}, jmax3={2*nmax+1}, ramp="flat{nmax}"'
        export = 'fn_fewbody_wave_function="case"' if export_dir is not None and rank >= 3 else ""
        (dir_ / "input.dat").write_text(f'''&input
 rank={rank}, hw={hw}, renorm="bare", Operators="hamil"
 NNInt="N3LO_EM500", input_nn_file="{interaction}"
 bra={sector}, ket={sector}, trans2lab=.false.
 NN_only=.true., genuine_3bf=.false., coul=.true.
 {truncation}
 {export}
/
''')
        result = subprocess.run([exe, "input.dat"], cwd=dir_, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        matches = re.findall(ENERGY[rank], result.stdout, flags=re.MULTILINE)
        if result.returncode or len(matches) != 1:
            raise RuntimeError(f"{label}: expected one ground-state energy\n{result.stdout[-6000:]}")
        version = re.search(r"^# Version is (.*)$", result.stdout, re.MULTILINE)
        if version is None:
            raise RuntimeError("NuHamil did not report its build revision")
        energy = float(matches[0])
        row = {"nucleus": label, "NuHamil": version[1],
               "hw_MeV": hw, "Nmax": nmax, "ground_MeV": energy}
        if export:
            wf, = dir_.glob("wf_case*.wav")
            h_file, = dir_.glob("Hamil_case*.op")
            vector = np.fromfile(wf, dtype=np.float64)
            h = np.fromfile(h_file, dtype=np.float64).reshape((vector.size, vector.size), order="F")
            residual = float(np.linalg.norm(h@vector-energy*vector))
            assert np.isclose(np.linalg.norm(vector),1.) and residual < 5e-8, residual
            assert export_dir is not None
            with (export_dir/f"{label}.npz").open("xb") as f:
                np.savez(f, ground_state=vector, hamiltonian_MeV=h, energy_MeV=energy,
                         hw_MeV=hw, nmax=nmax, sector=sector, revision=version[1])
            row.update(jacobi_dimension=vector.size, export_residual_MeV=residual)
        print(json.dumps(row), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--interaction", type=Path, required=True)
    p.add_argument("--export-dir", type=Path, help="New directory for A>=3 bound eigenpairs")
    p.add_argument("case", nargs="*", default=list(CASES))
    args = p.parse_args()
    source = args.interaction.resolve(strict=True)
    with source.open("rb") as f:
        print(json.dumps({"interaction_sha256": hashlib.file_digest(f, "sha256").hexdigest()}), flush=True)
    if args.export_dir is not None:
        args.export_dir.mkdir()
    for label in args.case:
        run_case(args.exe.resolve(strict=True), source, label, args.export_dir)


if __name__ == "__main__":
    main()
