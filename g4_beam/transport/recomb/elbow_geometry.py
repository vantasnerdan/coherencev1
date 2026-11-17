"""Two explicit Pd111 approach charts, neither a global reaction path."""
import numpy as np
from ase import Atoms


def elbow_pair(clean, separation_A, height_A, site='bridge'):
    """Keep the saved clean host and the declared molecular approach chart."""
    if (len(clean) != 6 or set(clean.numbers) != {46}
            or tuple(clean.pbc) != (True, True, False)):
        raise ValueError('This chart requires the declared three-layer 1x2 Pd111 host')
    if separation_A <= 0 or height_A <= 0:
        raise ValueError('Positive bond separation and surface height are required')
    primitive_a = np.asarray(clean.cell[0])
    primitive_b = np.asarray(clean.cell[1]) / 2
    # ASE's top-layer origin is Pd4. Bridge spans fcc/hcp hollows;
    # top follows a nearest-neighbour Pd axis toward bridge-site fragments.
    origin = clean.positions[4].copy()
    if site == 'bridge':
        direction = primitive_a + primitive_b
        center = origin + direction / 2
    elif site == 'top':
        direction = primitive_b.copy()
        center = origin.copy()
    else:
        raise ValueError('Unknown molecular approach chart')
    direction /= np.linalg.norm(direction)
    center[2] = float(np.max(clean.positions[:, 2])) + height_A
    positions = np.array([center - separation_A * direction / 2,
                          center + separation_A * direction / 2])
    atoms = clean.copy()
    atoms.calc = None
    atoms += Atoms('HH', positions=positions, masses=[2.014101778] * 2)
    return atoms
