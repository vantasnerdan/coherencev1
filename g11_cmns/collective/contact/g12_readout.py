"""Conditional proton readout for the fixed G12 PdD, 0.25-keV CM row."""

from math import isfinite, log10

# Frozen Part 1 proton-only, allowed-spin s-wave coefficient and own n65
# readout; see ../part1/rates.csv and the derivation in readout_contract.md.
_PROTON_A3_PER_S = 58122036.7098898
_OWN_N65_READOUT_PER_C = 434741.4188244709
_CUP_FOCUSING_FRACTION = 1.0
_CUP_CURRENT_A = 38e-6
_REF_CUP_CURRENT_A = 238e-6
TRUE_CURRENT_A = _CUP_CURRENT_A / _CUP_FOCUSING_FRACTION  # chart model, not a run measurement

_ELEMENTARY_CHARGE_C = 1.602176634e-19
_D_PER_BEAM_CHARGE = 0.85 + 2 * 0.15
_LOG_PROTON_A3_PER_S = log10(_PROTON_A3_PER_S)
_LOG_TRUE_CURRENT_A = log10(TRUE_CURRENT_A)
_LOG_D_PER_TRUE_C = log10(_D_PER_BEAM_CHARGE / _ELEMENTARY_CHARGE_C)
_LOG_K0_DENOMINATOR = log10(_CUP_FOCUSING_FRACTION * _OWN_N65_READOUT_PER_C)
_LOG_K1_FACTOR = log10(_REF_CUP_CURRENT_A / _CUP_CURRENT_A)


def _readout(log_per_s: float, log_per_c: float) -> dict[str, float]:
    log_delta_r_k0 = log_per_c - _LOG_K0_DENOMINATOR
    return {
        "log10_emitted_p_per_s": log_per_s,
        "log10_emitted_p_per_true_C": log_per_c,
        "log10_delta_R_k0": log_delta_r_k0,
        "log10_delta_R_k1": log_delta_r_k0 + _LOG_K1_FACTOR,
    }


def rate_readout(log10_contact_A3: float) -> dict[str, float]:
    """Convert positive contact in Å⁻³ for one continuously active normalized pair."""
    if not isfinite(log10_contact_A3):
        raise ValueError("expected finite log10 of positive contact in Å⁻³")
    log_per_s = log10_contact_A3 + _LOG_PROTON_A3_PER_S
    return _readout(log_per_s, log_per_s - _LOG_TRUE_CURRENT_A)


def event_readout(log10_excess_Xi_A3_s: float) -> dict[str, float]:
    """Convert positive excess in Å⁻³ s, supplied after one reference subtraction.

    The beam readout is conditional on one such preparation per incident D,
    with the same supplied excess for each beam component; no capture or
    preparation probability is inferred. Per-second output uses model current.
    """
    if not isfinite(log10_excess_Xi_A3_s):
        raise ValueError("expected finite log10 of positive integrated excess in Å⁻³ s")
    log_per_c = log10_excess_Xi_A3_s + _LOG_PROTON_A3_PER_S + _LOG_D_PER_TRUE_C
    return _readout(log_per_c + _LOG_TRUE_CURRENT_A, log_per_c)
