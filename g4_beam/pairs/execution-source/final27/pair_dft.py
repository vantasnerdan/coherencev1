"""SCF discriminator on the same neutral Pd10/H1 PBE surface."""
from math import ceil
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import serial_comm


def calculator(atoms, settings, log):
    electrons = sum(10 if symbol == 'Pd' else 1 for symbol in atoms.symbols)
    return GPAW(legacy_gpaw=False, communicator=serial_comm,
                mode=PW(settings['ecut_eV']), xc='PBE',
                setups={'Pd': '10', 'H': 'paw'},
                kpts={'size': settings['kmesh'], 'gamma': True},
                spinpol=True, charge=0,
                nbands=ceil(electrons / 2) + max(12, ceil(electrons * 0.08)),
                occupations=FermiDirac(settings['occupation_width_eV']),
                symmetry={'point_group': True, 'time_reversal': True},
                convergence={'energy': 1e-7, 'density': 1e-6, 'eigenstates': 1e-9},
                maxiter=150, parallel={'domain': 1, 'band': 1, 'kpt': 1},
                txt=str(log))
