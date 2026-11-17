"""Qualified serial new-PW GPAW, matched slab energetics."""
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import serial_comm


def calculator(atoms, log, cutoff=300., kmesh=(6, 3, 1)):
    isolated = all(s == 'H' for s in atoms.symbols)
    atom = isolated and len(atoms) == 1
    electrons = sum(10 if s == 'Pd' else 1 for s in atoms.symbols)
    return GPAW(legacy_gpaw=False, communicator=serial_comm,
                mode=PW(cutoff), xc='PBE', setups={'Pd': '10'},
                kpts={'size': (1, 1, 1) if isolated else kmesh, 'gamma': True},
                spinpol=atom, nbands=4 if isolated else (electrons+1)//2 + 12,
                occupations=FermiDirac(0.01 if isolated else 0.05),
                poissonsolver={} if isolated else {'dipolelayer': True},
                convergence={'energy': 1e-5, 'density': 1e-4, 'eigenstates': 1e-7},
                maxiter=100, parallel={'domain': 1, 'band': 1, 'kpt': 1}, txt=str(log))
