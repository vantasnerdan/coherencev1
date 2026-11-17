"""ASE Vibrations force recording with unchanged raw clamp reactions."""
from pathlib import Path
from ase.io import write
from ase.vibrations import Vibrations
from surface_io import manifest


class RecordedVibrations(Vibrations):
    def __init__(self, *args, bank_directory, **kwargs):
        super().__init__(*args, **kwargs)
        self.bank_directory = Path(bank_directory)
        self.bank_directory.mkdir(parents=True, exist_ok=True)

    def calculate(self, atoms, disp):
        forces = atoms.get_forces(apply_constraint=False)
        energy = float(atoms.get_potential_energy())
        checkpoint = self.bank_directory / f'{disp.name}.gpw'
        geometry = self.bank_directory / f'{disp.name}.traj'
        atoms.calc.write(str(checkpoint), mode='')
        write(geometry, atoms)
        return dict(forces=forces, energy_eV=energy,
                    positions_A=atoms.positions.copy(),
                    gpw=manifest(checkpoint), geometry=manifest(geometry))
