"""Atomic, non-overwriting numeric artifacts and raw-byte provenance."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def manifest(path, name=None):
    path = Path(path)
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            size += len(chunk)
            digest.update(chunk)
    return dict(path=name or path.name, bytes=size, sha256=digest.hexdigest())


def publish(path, writer):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    require(not path.exists(), f'Output already exists: {path.name}')
    fd, temporary = tempfile.mkstemp(prefix='.'+path.name+'.', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            writer(stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return manifest(path)


def write_json(path, value):
    data = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()
    return publish(path, lambda stream: stream.write(data))


def write_npz(path, **arrays):
    return publish(path, lambda stream: np.savez_compressed(stream, **arrays))
