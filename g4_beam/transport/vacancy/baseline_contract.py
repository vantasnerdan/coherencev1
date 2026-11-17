"""Extract literal transport parameters and exact consuming expressions from source bytes."""
import ast
import hashlib
from pathlib import Path


def extract(path):
    data = Path(path).read_bytes()
    text = data.decode()
    tree = ast.parse(text)
    parameters = None
    kb = None
    for node in tree.body:
        if isinstance(node, ast.Assign):
            names = [target.id for target in node.targets if isinstance(target, ast.Name)]
            if 'PD' in names:
                parameters = {kw.arg: ast.literal_eval(kw.value) for kw in node.value.keywords}
            if 'KB' in names:
                kb = ast.literal_eval(node.value)
    if parameters is None or kb is None:
        raise ValueError('Original PD/KB literal definitions missing')
    expressions = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name in ('v_eq', 'mu0', 'mu', 'foil_temperature'):
            expressions[node.name] = ast.get_source_segment(text, node)
    for line in text.splitlines():
        if any(token in line for token in ('theta_eq =', 'cap = m["mtrap"]', 'ktr =', 'R = ktr')):
            expressions.setdefault('trap_lines', []).append(line.strip())
    return dict(source=dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest()),
                PD=parameters, KB_eV_K=kb, expressions=expressions,
                units=dict(Eb='eV per D transferred from dilute solution',
                           mtrap='D per missing metal atom', Sf_v='k_B',
                           mu='eV per D relative to half D2 gas at1Pa'),
                interpretation='Independent equivalent traps: xt_eq=mtrap*v*logistic((mu-mu0+Eb)/kT). No ionic vibrational terms in trap model; phi_v uses the same independent-trap partition function.')
