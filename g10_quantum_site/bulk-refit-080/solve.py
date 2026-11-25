"""Replay .80 PdD using the accepted .48 eigenstate/Fourier implementation."""
import argparse
import json
import well

solver = well.load_module('pdd048_solver', well.PRIOR / 'solve.py')
solver.OUT = well.HERE / 'results'

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variant', choices=well.VARIANTS, required=True)
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    solver.OUT.mkdir(exist_ok=True)
    tag = 'check' if args.check else 'main'
    npar, order = (12,64) if args.check else (10,48)
    s = solver.solve('D', args.variant, npar, order)
    solver.save_spectrum(s, solver.OUT / f'D_{args.variant}_{tag}_spectrum.npz')
    rows = solver.analyze(s, tag, full=not args.check)
    (solver.OUT / f'summary_{args.variant}_{tag}.json').write_text(json.dumps(rows, indent=2)+'\n')

if __name__ == '__main__':
    main()
