"""Read-only lookup in the accepted PdD radial spectra and conditional size scenarios."""

from bisect import bisect_left
from dataclasses import dataclass
from fractions import Fraction
import json
from math import ceil, fsum, isfinite, log10
from pathlib import Path
from runpy import run_path


_CONTACT = Path(__file__).resolve().parent.parent / "contact"
_rate_readout = run_path(str(_CONTACT / "g12_readout.py"))["rate_readout"]
_TARGETS_EV = (100, 1000)
_ENERGIES_EV = (0.0006, 0.0048, 0.025)
_SIZES = (32, 256, 2048)


@dataclass(frozen=True)
class ContactSpectrum:
    join_A: float
    include_flagged: bool
    ground_energy_eV: float
    excitation_eV: tuple[float, ...]
    log10_contact: tuple[float, ...]
    contact_nondecreasing: bool


def load_contact_spectra(path: Path | str = _CONTACT / "results.json") -> tuple[ContactSpectrum, ...]:
    """Read native retained spectra only, not thermal observables or refined-model metadata."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    spectra = []
    for variant in data["variants"]:
        model, spectrum = variant["model"], variant["spectrum"]
        energy = spectrum["energy_eV"]
        contact = tuple(spectrum["log10_contact_A_minus3"])
        if not energy or len(energy) != len(contact) or not all(map(isfinite, (*energy, *contact))):
            raise ValueError("paired finite archived energy/contact arrays required")
        excitation = tuple(e - energy[0] for e in energy)
        if any(hi <= lo for lo, hi in zip(excitation, excitation[1:])):
            raise ValueError("archived excitation energies must strictly increase")
        spectra.append(ContactSpectrum(
            model["join_A"], model["include_flagged"], energy[0], excitation, contact,
            all(b >= a for a, b in zip(contact, contact[1:])),
        ))
    if not spectra:
        raise ValueError("no retained contact variants")
    return tuple(spectra)


def _peak_is_supported(row: dict) -> bool:
    exit_time = row["first_site_or_bond_exit_fs"]
    return row["candidate_supported"] is True and (exit_time is None or row["peak_time_fs"] < exit_time)


def _level(spectrum: ContactSpectrum, index: int) -> dict:
    contact = spectrum.log10_contact[index]
    return {
        "level_index": index,
        "excitation_eV": spectrum.excitation_eV[index],
        "log10_contact_A_minus3": contact,
        **_rate_readout(contact),
    }


def _span(variants: list[dict], side: str) -> dict:
    levels = [variant[side] for variant in variants]
    keys = ("log10_contact_A_minus3", "log10_emitted_p_per_s", "log10_emitted_p_per_true_C",
            "log10_delta_R_k0", "log10_delta_R_k1")
    return {key: [min(level[key] for level in levels), max(level[key] for level in levels)]
            for key in keys}


def conditional_readout(row: dict, spectra: tuple[ContactSpectrum, ...]) -> dict:
    """Assign the classical peak K_rel as excitation above each model's own ground.

    The endpoints are two stationary normalized l=0 eigenlevel readouts, not a
    population or a probability bound. Missing support is unresolved, not zero.
    """
    delta = row["peak_relative_kinetic_eV"]
    result = {
        "assigned_delta_E_eV": delta,
        "status": "unresolved_lattice_support",
        "variants": [],
        "across_join_endpoint_span": None,
    }
    if not _peak_is_supported(row):
        return result
    if not isfinite(delta) or delta < 0:
        raise ValueError("expected finite nonnegative supported pair kinetic energy")
    variants = result["variants"]
    for spectrum in spectra:
        variant = {
            "join_A": spectrum.join_A,
            "include_flagged": spectrum.include_flagged,
            "ground_energy_eV": spectrum.ground_energy_eV,
            "max_excitation_eV": spectrum.excitation_eV[-1],
            "contact_nondecreasing": spectrum.contact_nondecreasing,
            "status": "unresolved_contact_support",
            "lower": None,
            "upper": None,
        }
        variants.append(variant)
        j = bisect_left(spectrum.excitation_eV, delta)
        if j == len(spectrum.excitation_eV):
            continue
        i = j if spectrum.excitation_eV[j] == delta else j - 1
        variant.update(status="supported", lower=_level(spectrum, i), upper=_level(spectrum, j))
    if all(variant["status"] == "supported" for variant in variants):
        result.update(status="supported", across_join_endpoint_span={
            side: _span(variants, side) for side in ("lower", "upper")
        })
    else:
        result["status"] = "unresolved_contact_support"
    return result


def energy_budget_floors() -> dict:
    """Perfect-concentration floors for the campaign's exact decimal launch energies."""
    return {
        "N_at_fixed_energy_per_D": [
            {"target_eV": t, "energy_per_D_eV": e, "n_D_min": ceil(t / Fraction(str(e)))}
            for t in _TARGETS_EV for e in _ENERGIES_EV
        ],
        "energy_per_D_at_fixed_N": [
            {"target_eV": t, "n_D": n, "energy_per_D_eV_min": t / n}
            for t in _TARGETS_EV for n in _SIZES
        ],
    }


def conditional_size_scaling(rows: list[dict]) -> dict:
    """Compare one pattern and scope at fixed energy and identical selection/loss rules.

    Supply only the distributed-launch series, not initial pair-allocation
    controls. An unsupported size is never fitted or counted as zero. The
    linear-efficiency scenario uses the largest supported size independently
    of any fit; the power law requires growth at all three comparable sizes.
    """
    if not rows:
        raise ValueError("supply a single size series")
    group = {(r["pattern"], r["scope"], r["energy_per_D_eV"]) for r in rows}
    if len(group) != 1 or len({r["n_D"] for r in rows}) != len(rows):
        raise ValueError("size scaling requires a unique row per size for one pattern, scope and energy")
    pattern, scope, e = next(iter(group))
    supported = sorted((r for r in rows if _peak_is_supported(r)), key=lambda r: r["n_D"])
    result = {
        "pattern": pattern,
        "scope": scope,
        "energy_per_D_eV": e,
        "supported_n_D": [r["n_D"] for r in supported],
        "saturation_limit": "unresolved",
        "constant_efficiency_scenario": None,
        "power_law_scenario": {"status": "insufficient_supported_sizes"},
    }
    if not supported:
        return result
    n, peak = supported[-1]["n_D"], supported[-1]["peak_relative_kinetic_eV"]
    efficiency = peak / (n * e)
    if 0 < efficiency <= 1:
        result["constant_efficiency_scenario"] = {
            "interpretation": "assumed constant, not demonstrated at larger N",
            "anchor_n_D": n,
            "peak_energy_fraction": efficiency,
            "targets": [
                {"target_eV": t, "n_D_if_constant_efficiency": t / (efficiency * e)}
                for t in _TARGETS_EV
            ],
        }
    if len(supported) < 3:
        return result
    peaks = [r["peak_relative_kinetic_eV"] for r in supported]
    if any(not 0 < p <= r["n_D"] * e for p, r in zip(peaks, supported)):
        result["power_law_scenario"] = {"status": "nonpositive_peak_or_energy_budget_violation"}
        return result
    x = [log10(r["n_D"]) for r in supported]
    log_e = log10(e)
    y = [log10(p) - log_e for p in peaks]
    local = [(v - u) / (b - a) for a, b, u, v in zip(x, x[1:], y, y[1:])]
    if min(local) <= 0:
        result["power_law_scenario"] = {"status": "nonincreasing_step", "local_exponents": local}
        return result
    x0, y0 = fsum(x) / len(x), fsum(y) / len(y)
    dx = [a - x0 for a in x]
    alpha = fsum(a * (b - y0) for a, b in zip(dx, y)) / fsum(a * a for a in dx)
    targets = []
    for t in _TARGETS_EV:
        log_gain = log10(t) - log_e
        log_n = x0 + (log_gain - y0) / alpha
        budget_violation = log_n < log_gain
        targets.append({
            "target_eV": t,
            "status": "violates_total_energy_budget" if budget_violation else "conditional_only",
            "log10_n_D_if_power_law": log_n,
        })
    result["power_law_scenario"] = {
        "status": "finite_size_fit_not_a_prediction",
        "alpha": alpha,
        "local_exponents": local,
        "pivot_log10_n_D": x0,
        "pivot_log10_pair_gain": y0,
        "max_abs_log10_residual": max(abs(b - y0 - alpha * a) for a, b in zip(dx, y)),
        "targets": targets,
    }
    return result
