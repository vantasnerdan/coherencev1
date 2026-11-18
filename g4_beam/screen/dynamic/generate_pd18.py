"""Portable Pd18 generator adapter; source identities are in ledger.json."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np
import gpaw
from gpaw.atom import generator2
from g4_beam.screen.dynamic.pd18_dataset import binding, phase_diagnostics, read_identity




def scan_logarithmic_derivatives(gen, output, scan_max_Hartree):
    energies = np.linspace(-5., scan_max_Hartree, int(round((scan_max_Hartree+5.)/.05))+1)
    curves = dict(energy_Hartree=energies)
    for ell in range(4):
        ae = gen.aea.logarithmic_derivative(ell, energies, 3.)
        ps = gen.logarithmic_derivative(ell, energies, 3.)
        curves[f'ae_l{ell}'] = ae
        curves[f'pseudo_l{ell}'] = ps
    target = output/'logarithmic-derivatives.npz'
    np.savez_compressed(target, **curves)
    return phase_diagnostics(target)


def generate(output, radius_Bohr, scan_max_Hartree, projectors):
    parser = argparse.ArgumentParser()
    generator2.CLICommand.add_arguments(parser)
    recipe = ['Pd', '-f', 'PBE', '-e', '18', '-t', '18', '-w']
    if radius_Bohr is not None:
        recipe.extend(['-r', str(radius_Bohr)])
    if projectors is not None:
        recipe.extend(['-P', projectors])
    args = parser.parse_args(recipe)
    kwargs = generator2.get_parameters(args.symbol, args)
    old_cwd = Path.cwd()
    os.chdir(output)
    try:
        with (output/'generation.txt').open('x') as log:
            gen = generator2.generate(**dict(kwargs, output=log))
            if not gen.check_all():
                raise RuntimeError('Official default occupied bound-state check failed; no --no-check override')
            setup = gen.make_paw_setup(args.tag)
            setup.generatordata = ',\n    '.join(
                f'{key}={value!r}' for key, value in kwargs.items() if value is not None)
            setup.write_xml()
            diagnostic = scan_logarithmic_derivatives(gen, output, scan_max_Hartree)
        xml = output/'Pd.18.PBE'
        identity = read_identity(xml)
        return dict(schema='pd18-generator-candidate-v1',
                    status='generated_candidate_requiring_root_diagnostic_admission',
                    generator_recipe_argv=recipe, generator_kwargs=kwargs,
                    gpaw_version=gpaw.__version__, generator_source=binding(generator2.__file__),
                    producer_source=binding(__file__), default_occupied_bound_state_check=True,
                    default_check_limitation='Scalar-relativistic extra unoccupied-state errors are printed but not rejected by check_all; no blanket ghost-free claim.',
                    identity=identity, XML=binding(xml),
                    logarithmic_derivative_diagnostics=diagnostic,
                    generation_log=binding(output/'generation.txt'),
                    qualified_for_solid=False,
                    qualification='The integrator must inspect all atomic diagnostics before any dependent solid/response wave.')
    finally:
        os.chdir(old_cwd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--radius-Bohr', type=float)
    parser.add_argument('--scan-max-Hartree', type=float, default=10.)
    parser.add_argument('--projectors')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    status = 'failed'
    print('PD18_ATOMIC_START', datetime.now(timezone.utc).isoformat(), flush=True)
    try:
        result = generate(output, args.radius_Bohr, args.scan_max_Hartree, args.projectors)
        (output/'candidate.json').write_text(json.dumps(result, indent=2)+'\n')
        status = 'candidate_complete'
        print('PD18_ATOMIC_CANDIDATE_COMPLETE', flush=True)
    except Exception:
        traceback.print_exc()
        raise
    finally:
        terminal = dict(status=status, elapsed_wall_s=time.monotonic()-started,
                        process_cpu_s=time.process_time(), argv=sys.argv,
                        finished_at=datetime.now(timezone.utc).isoformat())
        (output/'terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')


if __name__ == '__main__':
    main()
