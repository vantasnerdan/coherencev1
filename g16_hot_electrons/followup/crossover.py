"""Run only the two commissioned primitive crossover temperatures."""
import argparse
from pathlib import Path

from native import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paw', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--temperature-K', type=int, choices=(1500, 2000), required=True)
    parser.add_argument('--mesh', type=int, choices=(12, 16), required=True)
    args = parser.parse_args()
    parent = args.root / f'k{args.mesh}' / str(args.temperature_K)
    for folder, direction, values in [('zero', 'zero', [0.]), ('100-fd', '100', [-0.01, 0.01])]:
        run(argparse.Namespace(paw=args.paw, output=parent / folder, mesh=args.mesh, temperature_K=args.temperature_K,
                               direction=direction, displacements=values))


if __name__ == '__main__':
    main()
