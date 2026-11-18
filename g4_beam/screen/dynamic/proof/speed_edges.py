"""Consumer-visible branch, shape and normalization edges for issue 243."""

from __future__ import annotations

import numpy as np
import pytest

from g4_beam.screen.dynamic._upstream.electronic_speed import (
    fermi_rms_velocity_au,
    speed_factor,
)


def test_continuous_endpoints_and_scalar_array_contract():
    near_one = np.array([np.nextafter(1., 0.), 1., np.nextafter(1., 2.)])
    np.testing.assert_allclose(speed_factor(near_one), 1/np.sqrt(2), rtol=1e-14, atol=0)
    assert speed_factor(0) == 1.0
    assert speed_factor(np.nextafter(0., 1.)) == 1.0
    assert isinstance(speed_factor(np.array(0.5)), float)
    values = np.array([[0., 0.5], [1., 2.]])
    assert speed_factor(values).shape == values.shape
    np.testing.assert_array_equal(speed_factor(values),
                                   np.array([[speed_factor(x) for x in row] for row in values]))
    assert speed_factor(np.empty((0, 2))).shape == (0, 2)


def test_high_speed_stays_positive_without_cancellation_clipping():
    ratios = np.array([1e8, 1e12, 1e16, 1e100, np.finfo(float).max])
    with np.errstate(invalid="raise", divide="raise", over="raise"):
        factors = speed_factor(ratios)
    assert np.all(np.isfinite(factors)) and np.all(factors > 0)
    np.testing.assert_allclose(ratios*factors, 1/np.sqrt(3), rtol=3e-14, atol=0)
    at_two = np.array([np.nextafter(2., 0.), 2., np.nextafter(2., 3.)])
    np.testing.assert_allclose(speed_factor(at_two),
                               np.sqrt(0.5-3*np.log(3)/8), rtol=3e-15, atol=0)


@pytest.mark.parametrize("ratio", [-1., -np.nextafter(0., 1.), np.inf, np.nan,
                                  [0., np.inf], 0.5+0j])
def test_speed_rejects_outside_domain(ratio):
    with pytest.raises(ValueError):
        speed_factor(ratio)


def test_velocity_trace_volume_and_offdiagonal_rotation():
    # Synthetic moment tensor, not material or free-electron-density data.
    moment = np.array([[4., 1., 0.], [1., 3., -0.5], [0., -0.5, 2.]])
    volume, dos = 13., 3.
    plasma = 4*np.pi/volume*moment
    expected = np.sqrt(3.)
    assert fermi_rms_velocity_au(plasma.astype(complex), dos, volume) == pytest.approx(expected)
    rotation = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    assert fermi_rms_velocity_au(rotation@plasma@rotation.T, dos, volume) == pytest.approx(expected)
    # Replicating a cell multiplies DOS and volume together, not omega_p^2.
    assert fermi_rms_velocity_au(plasma, 5*dos, 5*volume) == pytest.approx(expected)
    assert fermi_rms_velocity_au(np.zeros((3, 3)), dos, volume) == 0.


@pytest.mark.parametrize("plasma,dos,volume", [
    (np.eye(2), 1., 1.),
    (np.diag([1., 1., -0.1]), 1., 1.),
    (np.array([[1., 1., 0.], [0., 1., 0.], [0., 0., 1.]]), 1., 1.),
    (np.full((3, 3), np.nan), 1., 1.),
    (np.eye(3), 0., 1.),
    (np.eye(3), np.inf, 1.),
    (np.eye(3), 1., -1.),
    (np.eye(3), [1.], 1.),
    (np.eye(3), 1.+0j, 1.),
])
def test_velocity_rejects_nonplasma_or_mismatched_scalar_domains(plasma, dos, volume):
    with pytest.raises(ValueError):
        fermi_rms_velocity_au(plasma, dos, volume)
