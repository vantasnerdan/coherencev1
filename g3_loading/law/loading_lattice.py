"""PdD beta-phase lattice estimates at 77 K, from Schirber & Morosin's abstract.

Only the publisher's original abstract HTML is in hand, not the full 1975 paper.
No individual measured points or numerical uncertainty are available here.
The log-linear interpolation is OUR construction from the abstract's rounded
slope and author-extrapolated endpoint, not a reported exact measurement.
This is not a room-temperature or zero-temperature equation of state.
"""

import json
from hashlib import sha256
from math import exp, isfinite
from pathlib import Path

SOURCE = Path(__file__).with_name("loading-lattice-source.json")
SOURCE_SHA256 = "027f5d80d8cc30026b23c0279f0d84f15ca14208edb82ad38109a228127456ee"


def a(x, *, isotope="D", temperature_K=77.0, phase="beta"):
    """Return a provenance-bearing estimate in angstrom; refuse unsupported scope.

    Study domain: 0.75 <= x <= 1.0; measured range remains 0.8 < x < 0.98.
    Campaign-authorized x=0.75 extension is an extrapolation; uncertainty is unknown.
    """
    if isotope != "D" or temperature_K != 77.0 or phase != "beta":
        raise ValueError("Source identity requires isotope D, beta phase, and 77 K.")
    x = float(x)
    if not isfinite(x) or not 0.75 <= x <= 1.0:
        raise ValueError("Unsupported loading: campaign range is 0.75 <= x <= 1.")
    raw = SOURCE.read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError("Numerical calibration facts do not match the pinned source.")
    calibration = json.loads(raw)
    return {
        "x": x,
        "a_angstrom": calibration["author_extrapolated_a_x1_A"] * exp(calibration["dln_a_dx"] * (x - 1.0)),
        "status": "author_extrapolated" if x == 1.0 else
                  "model_interpolated" if 0.8 < x < 0.98 else "model_extrapolated",
        "measured_point": False,
        "isotope": isotope,
        "phase": phase,
        "temperature_K": temperature_K,
        "source_uncertainty_angstrom": None,
        "uncertainty_status": "not reported in accessible abstract; not zero",
        "construction": "author x=1 endpoint" if x == 1.0 else "a(x)=4.084*exp(0.044*(x-1)); constant log-slope ansatz",
        "source": "10.1103/PhysRevB.12.117, publisher abstract (article pp.117-118)",
        "source_sha256": calibration["source"]["original_html_sha256"],
        "calibration_facts_sha256": SOURCE_SHA256,
        "source_gate": "primary abstract; campaign-authorized extrapolation; pointwise uncertainty unavailable",
    }


if __name__ == "__main__":
    import json

    outcomes = []
    for x in (1.0, 0.875, 0.75):
        try:
            outcomes.append(a(x))
        except ValueError as exc:
            outcomes.append({"x": x, "status": "refused", "reason": str(exc)})
    for identity in ({"isotope": "H"}, {"temperature_K": 300.0}):
        try:
            outcomes.append(a(0.875, **identity))
        except ValueError as exc:
            outcomes.append({"x": 0.875, "identity": identity, "status": "refused", "reason": str(exc)})
    print(json.dumps(outcomes, indent=2, allow_nan=False))
