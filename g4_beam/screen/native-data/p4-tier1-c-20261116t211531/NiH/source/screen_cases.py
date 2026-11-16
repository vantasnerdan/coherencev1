"""Explicit ordered bulk representatives; seed dimensions are not DFT results."""
import numpy as np
from ase import Atoms
from ase.build import bulk

PD_REFERENCE_A = 3.890
OCTAHEDRAL = np.array([[.5, .5, .5], [.5, 0, 0], [0, .5, 0], [0, 0, .5]])
SETUPS = {'Pd': ('10', 'Pd.10.PBE.gz', 10),
          'Rh': ('9', 'Rh.9.PBE.gz', 9),
          'Ag': ('11', 'Ag.11.PBE.gz', 11),
          'Ni': ('10', 'Ni.10.PBE.gz', 10),
          'Ti': ('paw', 'Ti.PBE.gz', 12),
          'H': ('paw', 'H.PBE.gz', 1)}


def geometry(atoms):
    return dict(symbols=atoms.get_chemical_symbols(),
                positions_A=atoms.positions.tolist(),
                cell_A=atoms.cell.array.tolist(),
                initial_magnetic_moments_muB=atoms.get_initial_magnetic_moments().tolist())


def case_record(name, atoms, relaxation, note, magnetic=False):
    symbols = atoms.get_chemical_symbols()
    if magnetic:
        atoms.set_initial_magnetic_moments([1.0 if s == 'Ni' else 0.0 for s in symbols])
    electrons = sum(SETUPS[s][2] for s in symbols)
    return dict(case=name, geometry=geometry(atoms),
                model='neutral collinear PBE / explicit valence PAW; static charge RPA',
                magnetic=magnetic, valence_electrons=electrons,
                setups={s: SETUPS[s][0] for s in sorted(set(symbols))},
                setup_files=[SETUPS[s][1] for s in sorted(set(symbols))],
                relaxation=dict(mode=relaxation, max_steps=60,
                                force_tolerance_eV_A=.02,
                                stress_tolerance_eV_A3=.0005),
                metal_indices=[a for a, s in enumerate(symbols) if s != 'H'],
                D_indices=[a for a, s in enumerate(symbols) if s == 'H'],
                note=note,
                isotope_convention='H nuclear charge/electrons represent D in Born-Oppenheimer electronic quantities; no isotope kinetic prefactor inferred')


def bulk_cases():
    cases = []
    for percent in [-2, -1, 0, 1, 2]:
        atoms = bulk('Pd', 'fcc', a=PD_REFERENCE_A*(1+percent/100))
        cases.append(case_record(
            f'Pd_hydro_{percent:+d}pct', atoms, 'fixed',
            'Fixed hydrostatic strain relative to empirical a=3.890 A, not a redefined PBE equilibrium; constrained nonmagnetic branch'))
    for percent in [-2, 2]:
        atoms = bulk('Pd', 'fcc', a=PD_REFERENCE_A)
        atoms.set_cell(atoms.cell.array@np.diag([1+percent/100, 1+percent/100, 1]),
                       scale_atoms=True)
        cases.append(case_record(
            f'Pd_biax001_{percent:+d}pct', atoms, 'normal_001',
            'Cartesian (001) in-plane strain relative to a=3.890 A; relax normal zz and internal positions; constrained nonmagnetic branch'))
    for solute in ['Rh', 'Ag']:
        atoms = bulk('Pd', 'fcc', a=PD_REFERENCE_A, cubic=True)
        atoms[0].symbol = solute
        cases.append(case_record(
            'Pd3'+solute+'_L12', atoms, 'hydrostatic',
            'Ordered L12 corner solute/three Pd face centres; cubic volume and internal coordinates relaxed; a=3.890 A is only the initial seed'))
    for numerator in [1, 2, 3]:
        x = numerator/4
        a = PD_REFERENCE_A+.2*x
        atoms = bulk('Pd', 'fcc', a=a, cubic=True)
        sites = OCTAHEDRAL[:1] if numerator == 1 else (
            OCTAHEDRAL[:2] if numerator == 2 else OCTAHEDRAL[1:])
        atoms += Atoms('H'*numerator, positions=sites*a)
        cases.append(case_record(
            f'PdD_x{x:g}', atoms, 'cell',
            'Ordered octahedral occupancy, not a disordered configurational average; D and full cell relaxed; linear seed interpolation is not an empirical loading law'))
    cases.append(case_record(
        'PdD_x1', bulk('PdH', 'rocksalt', a=4.09), 'hydrostatic',
        'Full octahedral rocksalt occupancy; cubic volume and D positions relaxed; a=4.09 A is only an initial seed'))
    cases.append(case_record(
        'Ni', bulk('Ni', 'fcc', a=3.52), 'hydrostatic',
        'Ferromagnetic collinear initialization; relaxed cubic volume; a=3.52 A is only an initial seed', True))
    cases.append(case_record(
        'NiH', bulk('NiH', 'rocksalt', a=3.75), 'hydrostatic',
        'Collinear magnetic initialization may relax to zero moment; cubic volume and H positions relaxed; a=3.75 A is only an initial seed', True))
    cases.append(case_record(
        'Ti_hcp', bulk('Ti', 'hcp', a=2.951, c=4.686), 'hexagonal',
        'hcp representative; a,c and internal positions relaxed; supplied dimensions are initial seeds'))
    cases.append(case_record(
        'TiD2_fluorite', bulk('TiH2', 'fluorite', a=4.44), 'hydrostatic',
        'Stoichiometric fluorite tetrahedral D sites, not the transport TiD1.85 estimate; cubic volume and D positions relaxed'))
    return cases
