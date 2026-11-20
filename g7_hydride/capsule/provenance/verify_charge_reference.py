"""Check only common charged-cell gauge/background cancellation, not image convergence."""
import hashlib
import json
from pathlib import Path
import sympy as sp


def verify():
    q, other_q = sp.symbols('q other_q', real=True)
    k, phi, constant = sp.symbols('k phi constant', real=True)
    volume = sp.symbols('volume', positive=True)
    correction = lambda charge: k*charge**2/volume**sp.Rational(1, 3)+phi*charge+constant
    # Probe occupied/vacant have the same declared q; both host references have q=0.
    common = sp.expand(correction(q)-correction(q)-correction(0)+correction(0))
    wrong = sp.expand(correction(q)-correction(other_q)-correction(0)+correction(0))
    assert common == 0 and all(common.subs(q, charge) == 0 for charge in (0, 1))
    wrong_branch = sp.simplify(wrong.subs({q: 1, other_q: 0}))
    assert wrong_branch != 0 and sp.diff(wrong_branch, k) != 0
    # Atomic energy-reference shifts are linear in nucleus count, including empty G.
    assert 2-1-1+0 == 0
    source = Path(__file__).read_bytes()
    return dict(schema='hydride-charge-reference-check-v1', passed=True,
                source=dict(bytes=len(source), sha256=hashlib.sha256(source).hexdigest()),
                matched_charge_expression=str(common),
                wrong_charge_expression=str(wrong_branch), wrong_charge_sensitivity_nonzero=True,
                gas_atomic_reference_count_coefficients=[2, -1, -1, 0],
                scope='An identical same-cell q-squared background/gauge term cancels in the quartet. This does NOT cancel geometry-dependent electronic/multipole images; actual F two-cell checks are still required.',
                empty_reference='With zero nuclei/electrons the electronic/ionic/atomic-reference sums are empty and the adopted Hamiltonian zero is explicit; no failed Native state is substituted.')


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
