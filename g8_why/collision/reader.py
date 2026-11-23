"""Verify the collision namespace and replay saved model evidence; no Native code.

Run: python collision/reader.py G8_WHY_DIRECTORY
Requires NumPy/SciPy. Other mechanism pages in the shared directory are outside
this capsule; every collision/ file and collision.md must be manifested.
"""
import argparse
from fractions import Fraction
import hashlib
import importlib.abc
import json
from pathlib import Path, PurePosixPath
import sys

sys.dont_write_bytecode = True


class NoNativeImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'gpaw', 'ase'}:
            raise ImportError('Native library forbidden in saved-output capsule: '+fullname)
        return None


sys.meta_path.insert(0, NoNativeImports())
import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def inventory(root):
    manifest = json.loads((root/'collision/manifest.json').read_bytes())
    require(manifest['schema'] == 'g8-collision-capsule-v1', 'Wrong collision capsule schema')
    require(manifest['authority_entry'] == 'entry_3259' and not manifest['publishing'], 'Wrong delivery authority')
    actual = {'collision.md'}|{p.relative_to(root).as_posix() for p in (root/'collision').rglob('*') if p.is_file()}
    require(actual == set(manifest['files'])|{'collision/manifest.json'}, 'Missing or unmanifested collision file')
    for name, item in manifest['files'].items():
        path = PurePosixPath(name)
        require(not path.is_absolute() and '..' not in path.parts and
                (name == 'collision.md' or name.startswith('collision/')), 'Unsafe capsule member')
        source = root/name
        require(not source.is_symlink() and source.is_file(), 'Link or missing capsule file')
        require(digest(source.read_bytes()) == item['EXPORTED'], 'EXPORTED byte mismatch: '+name)
    return manifest


def exact_monotonicity(common, report):
    domain = Fraction(16)
    for channel, coefficients in common.DD.BH.items():
        a = [Fraction(str(v)) for v in coefficients]
        lower = a[0]+sum(min(v,Fraction(0))*domain**i for i,v in enumerate(a[1:],1))
        slope = a[1]+sum(min(i*v,Fraction(0))*domain**(i-1) for i,v in enumerate(a[2:],2))
        certificate = report['symbolic']['BH_lower_bounds'][channel]
        require(lower == Fraction(certificate['S_lower_keV_mb']) and lower > 0, 'S coefficient certificate changed')
        require(slope == Fraction(certificate['S_prime_lower_mb']) and slope > 0, 'S derivative certificate changed')
    gap = Fraction(str(common.DD.BG))**2-domain-Fraction(str(common.U))
    require(gap == Fraction(report['symbolic']['BG_squared_minus_D_minus_U_keV']) and gap > 0,
            'Gamow monotonic domain changed')


def replay_envelope(root, common, report):
    differences = []
    with np.load(root/'collision/evidence/envelope-nodes.npz',allow_pickle=False) as nodes:
        for row in report['results']:
            energy = row['CM_keV']
            label = f"{row['material']}_{energy}_{int(row['nonnegative_common_removal'])}"
            lab,path = nodes[label+'_lab'],nodes[label+'_path']
            ref_lab,ref_path = nodes[label+'_ref_lab'],nodes[label+'_ref_path']
            cap = (np.sqrt(lab.max())+np.sqrt(2*lab.max()+.2))**2
            removal = lab.max()+.1 if row['nonnegative_common_removal'] else 0.
            lower_energy = (ref_lab-cap)/2-removal
            upper_rate = common.DD.sigma(lab+.1,common.U,'p')*np.sqrt(2*(lab+.1)/lab)
            lower_rate = np.zeros_like(ref_lab)
            positive = lower_energy > 0
            lower_rate[positive] = .5*common.DD.sigma(lower_energy[positive],common.U,'n')*np.sqrt(2*lower_energy[positive]/ref_lab[positive])
            require(np.allclose(upper_rate,nodes[label+'_upper_rate'],rtol=1e-11,atol=0), 'Saved numerator rate changed')
            require(np.allclose(lower_rate,nodes[label+'_lower_rate'],rtol=1e-11,atol=0), 'Saved reference rate changed')
            upper,lower = float(path@upper_rate),float(ref_path@lower_rate)
            require(upper > 0 and lower > 0, 'Nonpositive reported reference')
            def rates(e, channel):
                if e == energy and channel == 'p':
                    return upper
                if e == 6.5 and channel == 'n':
                    return lower
                raise ValueError('Unexpected readout energy/channel')
            for power in (0,1):
                value = float(common.readout(row['material'],[energy],rates,current_power=power)[0])
                difference = abs(value/row['reported_upper'][str(power)]-1)
                require(difference < 1e-10, 'Saved normalized readout changed')
                differences.append(difference)
    return dict(reported_values=len(differences),max_relative_difference=max(differences))


def replay_pilot(root, report):
    with np.load(root/'collision/evidence/joint-spectrum.npz',allow_pickle=False) as joint:
        probability,kinetic = joint['joint_probabilities'],joint['kinetic_eV']
        require(np.all(probability >= 0) and abs(probability.sum()-1) < 1e-10, 'Joint probability changed')
        marginal = probability.sum(axis=0)
        mean,second = float(marginal@kinetic),float(marginal@(kinetic*kinetic))
    expected = report['original_pilot']['saved_output_consumer']['thermal300K']
    require(abs(mean-expected['mean_kinetic_eV']) < 1e-11 and
            abs(second-expected['second_kinetic_moment_eV2']) < 1e-11, 'Finite-local pilot moments changed')
    return dict(mean_kinetic_eV=mean,second_moment_eV2=second,
                scope='Finite-local pilot only; physical Coulomb tail remains an explicit closed omission')


def replay(root):
    manifest = inventory(root)
    sys.path.insert(0,str(root/'collision/source'))
    from yield_kernel import source_common
    common = source_common(root/'collision/model')
    envelope = json.loads((root/'collision/evidence/envelope.json').read_bytes())
    pilot = json.loads((root/'collision/evidence/pilot.json').read_bytes())
    exact_monotonicity(common,envelope)
    return dict(schema='g8-collision-portable-consumer-v1',passed=True,
        manifested_files=len(manifest['files']),envelope=replay_envelope(root,common,envelope),
        pilot=replay_pilot(root,pilot),omissions=manifest['closed_omissions'],
        native_imports_forbidden=['gpaw','ase'],new_native_invocations=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(replay(args.root.resolve()),indent=2,sort_keys=True))
    except (ValueError,KeyError,OSError,AssertionError) as error:
        print('REFUSED: '+str(error),file=sys.stderr)
        sys.exit(2)
