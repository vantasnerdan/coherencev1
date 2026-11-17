"""Reanalyze authentic portable route banks without launching Native work."""
import argparse
from pathlib import Path
from elbow_path_results import analyze_path, bound_json
from surface_io import manifest, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Route reanalysis output is immutable')
    path = args.manifest.resolve()
    base = path.parent
    bundle = bound_json(dict(path=str(path), sha256=args.sha256))
    if bundle['schema'] != 'coherence-portable-elbow-inputs-v1':
        raise ValueError('Explicit portable elbow input schema required')
    def absolute(reference):
        return {**reference, 'path': str(base/Path(reference['path']))}
    segments = [{**segment, 'directory': str(base/segment['directory']),
                 'reference_base': str(base)} for segment in bundle['segments']]
    report = analyze_path(segments, absolute(bundle['target_route']),
        absolute(bundle['clean_reference']), absolute(bundle['gas_reference']))
    report.update(portable_input_manifest=manifest(path), Native_launched=False,
        peak_k_slab_convergence_qualified=False,
        transport_ledger_replacement_authorized=False,
        interpretation='Restricted-route reanalysis; not a global activation barrier')
    save(args.output, report)
    print('PORTABLE_ELBOW_REANALYSIS', report['complete_route'],
          report.get('sampled_peak_excess_eV'), flush=True)


if __name__ == '__main__':
    main()
