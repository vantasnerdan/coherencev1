"""Finite saved-output intake; ORIGINAL bindings never refer to redacted bytes."""
import hashlib
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def load_json(path):
    return json.loads(Path(path).read_bytes())


class Capsule:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.manifest = load_json(self.root / 'manifest.json')
        require(self.manifest['schema'] == 'pdd-collision-capsule-v1', 'Wrong capsule schema')
        self.entries = self.manifest['files']
        require(len(self.entries) == len(set(self.entries)), 'Duplicate manifest paths')
        actual = {str(p.relative_to(self.root)) for p in self.root.rglob('*') if p.is_file()}
        require(actual == set(self.entries) | {'manifest.json'}, 'Unmanifested or missing capsule files')
        for name, entry in self.entries.items():
            path = self.root / name
            require(path.resolve().is_relative_to(self.root) and not path.is_symlink(), 'Unsafe capsule member')
            require(digest(path.read_bytes()) == entry['EXPORTED'], 'EXPORTED byte integrity failed: ' + name)

    def json(self, name):
        return load_json(self.root / name)

    def original(self, name):
        return self.entries[name]['ORIGINAL']

    def bind(self, binding, name):
        original = self.original(name)
        require(binding['sha256'] == original['sha256'] and binding['bytes'] == original['bytes'],
                'ORIGINAL source binding mismatch: ' + name)

    def bound_file(self, binding, directory):
        name = directory + '/' + binding['path']
        self.bind(binding, name)
        return name
