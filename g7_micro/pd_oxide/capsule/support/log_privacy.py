"""Axiom-approved GPAW log metadata transforms; never scientific JSON/source."""
import re

USER = re.compile(rb'(?im)^[ \t]*User:[^\r\n]*(?:\r?\n|$)')
CWD = re.compile(rb'(?im)^[ \t]*CWD:[^\r\n]*(?:\r?\n|$)')
LIBRARIES = re.compile(rb'(?m)^Paths:[ \t]*\r?\n(?:[ \t]+[^\r\n]*(?:\r?\n|$))+')
DATASET = re.compile(rb'(?im)^([ \t]*file:[ \t]+)(/[^\r\n \t]+)')
CHECKPOINT = re.compile(rb'''(?im)^((?:Writing to|Reading from)[ \t]+)(/[^\s"'()]+\.gpw)(?=[ \t]|$)''')
LOG_TRANSFORM_COUNTS = frozenset({'gpaw_user_headers_removed', 'gpaw_cwd_headers_removed',
                                  'gpaw_library_path_blocks_removed', 'dataset_path_prefixes_removed',
                                  'gpw_path_prefixes_removed'})


def allowed_log_transform(changes):
    return (type(changes) is list and len(changes) == 1 and type(changes[0]) is dict and
            bool(changes[0]) and set(changes[0]) <= LOG_TRANSFORM_COUNTS and
            all(type(value) is int and value > 0 for value in changes[0].values()))


def sanitize_log(data):
    """Remove identity/location metadata, preserving dataset/checkpoint basenames."""
    data.decode('utf-8')
    counts = {}
    for name, pattern in (('gpaw_user_headers_removed', USER),
                          ('gpaw_cwd_headers_removed', CWD),
                          ('gpaw_library_path_blocks_removed', LIBRARIES)):
        data, count = pattern.subn(b'', data)
        if count:
            counts[name] = count
    for name, pattern in (('dataset_path_prefixes_removed', DATASET),
                          ('gpw_path_prefixes_removed', CHECKPOINT)):
        data, count = pattern.subn(lambda match: match[1]+match[2].rsplit(b'/', 1)[-1], data)
        if count:
            counts[name] = count
    return data, ([counts] if counts else [])
