"""One central-strain case at the accepted primitive mesh; zero strain is reused."""
import argparse
import json
import traceback

import numpy as np
from primitive import calculate

CASES = ("hydro_minus", "hydro_plus", "parallel_minus", "parallel_plus", "perpendicular_minus", "perpendicular_plus")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metal", choices=("Ni", "Pd"))
    parser.add_argument("kmesh", type=int, choices=(6, 8, 12, 16, 20))
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--cutoff", type=int, choices=(500, 600), default=500)
    args = parser.parse_args()
    eps = -0.01 if args.case.endswith("minus") else 0.01
    deformation = np.eye(3)*(1+eps) if args.case.startswith("hydro") else np.diag([1+eps, 1., 1.])
    axis = 1 if args.case.startswith("perpendicular") else 0
    name = f"k{args.kmesh:02d}_pw{args.cutoff}_{args.case}"
    result = calculate(args.metal, name, deformation, kmesh=args.kmesh, cutoff_eV=args.cutoff, axis=axis)
    print(json.dumps(result, indent=2, allow_nan=False), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
