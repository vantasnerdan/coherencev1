"""One requested primitive mesh/cutoff pair; no fit to the measured frequency."""
import argparse
import json
from pathlib import Path
import traceback

import numpy as np
from scipy.constants import e, hbar, m_p, physical_constants
from primitive import calculate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metal", choices=("Ni", "Pd"))
    parser.add_argument("kmesh", type=int, choices=(6, 8, 12, 16, 20))
    parser.add_argument("--cutoff", type=int, choices=(500, 600), default=500)
    args = parser.parse_args()
    name = f"k{args.kmesh:02d}_pw{args.cutoff}"
    result = calculate(args.metal, name, np.eye(3), kmesh=args.kmesh, cutoff_eV=args.cutoff)
    K = result["K_Gamma_eV_A2"]
    mass = m_p if args.metal == "Ni" else physical_constants["deuteron mass"][0]
    magnitude = float(1e3*hbar*np.sqrt(abs(K)*e/1e-20/mass)/e)
    result.update(harmonic_quantum_meV=magnitude if K > 0 else None,
                  imaginary_magnitude_meV=magnitude if K < 0 else None,
                  measured_peak_meV=89.0 if args.metal == "Ni" else 37.0,
                  mass_kg=mass)
    path = Path("primitive_"+args.metal)/name/"result.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result, indent=2, allow_nan=False), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
