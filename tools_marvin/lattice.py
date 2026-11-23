"""Exercise i-PI and SSCHA on the same archived cubic isotope well."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from xml.sax.saxutils import escape

from ase.units import Bohr, Hartree, kB
import numpy as np
from scipy.optimize import root_scalar


__DRIVER_CLASS__ = "CubicDriver"
__DRIVER_NAME__ = "cubic-well"
HERE = Path(__file__).resolve().parent
TEMPERATURE = 300.


def well(x, k, a, b):
    squared = x*x
    total = squared.sum(axis=-1, keepdims=True)
    energy = (0.5*k*total[...,0] + (a-b)/4 * (squared**2).sum(axis=-1) + b/4*total[...,0]**2)
    force = -(k + (a-b)*squared + b*total) * x
    return energy, force


class CubicDriver:
    def __init__(self, k, a, b, **_kwargs):
        self.k, self.a, self.b = k, a, b
        self.virial = np.zeros((3,3))

    def __call__(self, _cell, positions):
        energies, forces = well(np.asarray(positions)*Bohr, self.k, self.a, self.b)
        if isinstance(positions, list):
            return [(float(e.sum()/Hartree), f*(Bohr/Hartree), self.virial, "{}")
                    for e,f in zip(energies,forces)]
        return float(energies.sum()/Hartree), forces*(Bohr/Hartree), self.virial, "{}"


def input_well():
    data = json.loads((HERE.parent / "g2_refit/computed_inputs.json").read_text())
    d = next(row for row in data["materials"] if row["name"] == "PdD")
    a = d["K4_100_eV_A4"]
    b = (3*d["K4_111_eV_A4"]-a)/2
    return d["K2_eV_A2"], a, b, d["mass_eVc2"]


def ipi_case(k, a, b, mass):
    beads, steps, stride = 16, 20000, 10
    # Mass is supplied explicitly in electron masses, not inferred from the XYZ label.
    mass_au = mass / 510998.95069
    executable = Path(sys.executable).with_name("i-pi")
    with tempfile.TemporaryDirectory(prefix="cubic-pimd-") as tmp:
        work = Path(tmp)
        (work/"init.xyz").write_text("1\n# CELL(abcABC): 20 20 20 90 90 90 positions{angstrom} cell{angstrom}\nD 0 0 0\n")
        (work/"input.xml").write_text(f'''<simulation verbosity="low">
  <output prefix="well">
    <properties stride="{stride}" filename="out">[step, r_gyration{{angstrom}}, atom_x{{angstrom}}(0)]</properties>
  </output>
  <total_steps>{steps}</total_steps>
  <prng><seed>731</seed></prng>
  <ffdirect name="cubic">
    <pes>custom</pes><pes_path>{escape(str(Path(__file__).resolve()))}</pes_path>
    <batch_size>{beads}</batch_size><parameters>{{ k:{k}, a:{a}, b:{b} }}</parameters>
  </ffdirect>
  <system>
    <initialize nbeads="{beads}">
      <file mode="xyz">init.xyz</file>
      <masses mode="manual">[{mass_au}]</masses>
      <velocities mode="thermal" units="kelvin">{TEMPERATURE}</velocities>
    </initialize>
    <forces><force forcefield="cubic"/></forces>
    <ensemble><temperature units="kelvin">{TEMPERATURE}</temperature></ensemble>
    <motion mode="dynamics"><fixcom>False</fixcom><dynamics mode="nvt">
      <timestep units="femtosecond">0.5</timestep>
      <thermostat mode="pile_l"><tau units="femtosecond">20</tau></thermostat>
    </dynamics></motion>
  </system>
</simulation>''')
        run = subprocess.run([executable, "input.xml"], cwd=work, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, check=False)
        out = work/"well.out"
        if run.returncode or not out.is_file():
            raise RuntimeError(run.stdout[-6000:])
        rows = np.loadtxt(out)
    assert len(rows) == steps//stride+1 and rows.shape[1] == 5, rows.shape
    # Bead-mean r^2 = centroid r^2 + radius_of_gyration^2.
    samples = ((rows[-1600:,1:]**2).sum(axis=1))/3
    block = samples.reshape(16,-1).mean(axis=1)
    moment = float(block.mean())
    sem = float(block.std(ddof=1) / np.sqrt(len(block)))
    data = json.loads((HERE.parent/"g2_refit/anharmonic.json").read_text())
    spectrum = next(row["spectrum"] for row in data["rows"] if row["material"] == "PdD")
    reference = spectrum["local_mean_square_axis_A2_T300"]
    assert abs(moment-reference) < 0.002, (moment, reference, sem)
    print(json.dumps({"case": "i-PI isolated D", "beads": beads, "temperature_K": TEMPERATURE,
                      "mean_square_axis_A2": moment, "block_standard_error_A2": sem,
                      "reference_mean_square_axis_A2": reference}, indent=2))


def sscha_case(k, a, b, mass):
    from ase import Atoms
    from ase.calculators.calculator import Calculator, all_changes
    from cellconstructor import Phonons, Structure, Units
    from sscha.Ensemble import Ensemble

    # Relative mass of two equal 2*m particles is m; translations stay zero.
    sep = np.array([2.,0.,0.])
    atoms = Atoms("HH", positions=[[0,0,0],sep], cell=[20.,20.,20.],
                 pbc=True, masses=np.full(2, 2*mass/931494102.42))
    structure = Structure.Structure()
    structure.generate_from_ase_atoms(atoms)
    projector = np.kron([[1.,-1.],[-1.,1.]],np.eye(3))

    class CalculatorWell(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            assert self.atoms is not None
            vector = self.atoms.positions[1]-self.atoms.positions[0]-sep
            energy, force = well(vector,k,a,b)
            self.results = {"energy": float(energy), "forces": np.array([-force,force])}

    def stochastic_gradient(trial):
        np.random.seed(587)
        dyn = Phonons.Phonons(structure, nqirr=1, force_real=True)
        dyn.dynmats[0] = trial * projector / (Units.RY_TO_EV * Units.A_TO_BOHR**2)
        ensemble = Ensemble(dyn, TEMPERATURE, supercell=(1,1,1))
        ensemble.generate(1024, evenodd=True)
        ensemble.compute_ensemble(CalculatorWell(), compute_stress=False, verbose=False)
        grad = ensemble.get_preconditioned_gradient(verbose=False)[0].real
        return float(np.sum(projector * grad))

    def exact_gradient(trial):
        h2m = 1973.269804**2/mass  # (hbar*c)^2/mc^2, eV Å^2
        quantum = np.sqrt(h2m*trial)
        spread = h2m/(2*quantum*np.tanh(quantum/(2*kB*TEMPERATURE)))
        return trial-k-(3*a+2*b)*spread

    analytic = root_scalar(exact_gradient, bracket=(k, 2*k)).root
    fitted = root_scalar(stochastic_gradient, bracket=(k, 2*k), xtol=1e-4).root
    assert abs(fitted-analytic) < 0.10, (fitted, analytic)
    print(json.dumps({"case": "SSCHA dimer, relative mode", "configurations_per_gradient": 1024,
                      "temperature_K": TEMPERATURE, "spring_eV_A2": fitted,
                      "analytic_SCHA_spring_eV_A2": analytic}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("ipi", "sscha"))
    args = parser.parse_args()
    {"ipi": ipi_case, "sscha": sscha_case}[args.case](*input_well())


if __name__ == "__main__":
    main()
