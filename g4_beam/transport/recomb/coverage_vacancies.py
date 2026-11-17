"""Hold a saved finite-time vacancy profile fixed in the bulk-fed steady closure.

Cell-average remapping conserves integral(v dz); v is vacancies/Pd, not D/Pd.
The primary source, condition, energy, run time and temperature remain bound.
"""
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import numpy as np


@dataclass(frozen=True)
class FrozenVacancies:
    edges_m: np.ndarray
    fraction_per_Pd: np.ndarray
    provenance: dict

    def on_grid(self, target_edges_m):
        target = np.asarray(target_edges_m, dtype=float)
        if (target.ndim != 1 or len(target) < 2 or np.any(~np.isfinite(target))
                or np.any(np.diff(target) <= 0)
                or abs(target[0] - self.edges_m[0]) > 1e-15
                or abs(target[-1] - self.edges_m[-1]) > 1e-15):
            raise ValueError('Target grid must cover the same complete metal foil')
        integral = np.concatenate(([0.], np.cumsum(
            self.fraction_per_Pd * np.diff(self.edges_m))))
        cumulative = np.interp(target, self.edges_m, integral)
        remapped = np.diff(cumulative) / np.diff(target)
        return remapped


def read_run_end(path, expected_sha256, *, condition, energy_cm_keV,
                 temperature_K, source_revision):
    path = Path(path)
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != expected_sha256:
        raise ValueError('Run-end vacancy primary-source identity changed')
    data = json.loads(raw)
    matches = [(i, row) for i, row in enumerate(data['records'])
               if row['case']['condition'] == condition
               and row['case']['E_cm_keV'] == energy_cm_keV
               and abs(row['case']['T_K'] - temperature_K) < 1e-10]
    if len(matches) != 1:
        raise ValueError('Exactly one accepted run-end operating point required')
    index, row = matches[0]
    if row['state']['d'] != 0:
        raise ValueError('This frozen-v closure requires the declared bare Pd foil')
    layers = np.asarray(row['layers'], dtype=float)
    vacancy = np.asarray(row['state']['v'], dtype=float)
    edges = np.concatenate((layers[:, 0], [layers[-1, 1]]))
    if (len(vacancy) != len(layers) or np.any(~np.isfinite(vacancy))
            or np.any(np.diff(edges) <= 0)
            or not np.allclose(layers[:-1, 1], layers[1:, 0], atol=1e-15, rtol=0)
            or not np.allclose((edges[:-1] + edges[1:]) / 2,
                               row['state']['zc'], atol=1e-15, rtol=0)):
        raise ValueError('Run-end cell geometry or vacancy fractions are inconsistent')
    edges.setflags(write=False)
    vacancy.setflags(write=False)
    return FrozenVacancies(edges, vacancy, dict(
        primary_path=str(path), primary_bytes=len(raw), primary_sha256=digest,
        source_revision=source_revision, record_index=index, case=row['case'],
        raw_minimum_v_per_Pd=float(vacancy.min()),
        raw_negative_cells=int(np.sum(vacancy < 0)),
        model_vacancy_projection='Original transport.py rhs uses max(raw_v,0); raw saved profile is retained',
        ensemble='Solute/trap/sheet steady state with run-end v(z) held fixed; no vacancy thermal equilibration'))
