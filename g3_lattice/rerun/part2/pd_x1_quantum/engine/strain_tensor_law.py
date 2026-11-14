"""Cubic-covariant PdD Gamma-stiffness closure in the reference material frame.

Input is the symmetric Biot strain U-I, not engineering shear. Calibration and
measured residuals live in strain-law.json; no eigenvalue clipping is performed.
"""
import numpy as np

NAMES = (
    "normal_parallel", "normal_transverse", "shear",
    "normal_parallel_squared", "normal_parallel_times_transverse_sum",
    "normal_transverse_squares", "normal_transverse_product",
    "shear_squares_touching_axis", "shear_square_opposite_axis",
    "shear_times_pair_normal_sum", "shear_times_opposite_normal",
    "cross_shear_product",
)


def feature_tensor(strain):
    """Return the complete linear/quadratic cubic basis, (...,12,3,3)."""
    E = np.asarray(strain)
    b = np.zeros(E.shape[:-2]+(12, 3, 3), dtype=E.dtype)
    for i in range(3):
        j, k = [axis for axis in range(3) if axis != i]
        ei, ej, ek = E[..., i, i], E[..., j, j], E[..., k, k]
        b[..., 0, i, i] = ei
        b[..., 1, i, i] = ej+ek
        b[..., 3, i, i] = ei*ei
        b[..., 4, i, i] = ei*(ej+ek)
        b[..., 5, i, i] = ej*ej+ek*ek
        b[..., 6, i, i] = ej*ek
        b[..., 7, i, i] = E[..., i, j]**2+E[..., i, k]**2
        b[..., 8, i, i] = E[..., j, k]**2
    for i, j in ((0, 1), (0, 2), (1, 2)):
        k = 3-i-j
        s = E[..., i, j]
        b[..., 2, i, j] = b[..., 2, j, i] = s
        b[..., 9, i, j] = b[..., 9, j, i] = (E[..., i, i]+E[..., j, j])*s
        b[..., 10, i, j] = b[..., 10, j, i] = E[..., k, k]*s
        b[..., 11, i, j] = b[..., 11, j, i] = E[..., i, k]*E[..., j, k]
    return b


def calibration_mask(strain, model):
    """Bounding-box eligibility, not certification of all mixed strain states."""
    E = np.asarray(strain)
    diagonal = np.diagonal(E, axis1=-2, axis2=-1)
    off = np.stack([E[..., 0, 1], E[..., 0, 2], E[..., 1, 2]], axis=-1)
    bounds = model["calibration_box"]
    return ((diagonal >= bounds["normal_min"]-1e-12).all(axis=-1)
            & (diagonal <= bounds["normal_max"]+1e-12).all(axis=-1)
            & (np.abs(off) <= bounds["shear_abs_max"]+1e-12).all(axis=-1))


def stiffness_tensor(strain, model, *, extrapolate=False):
    """Return material-frame collective K_Gamma in eV/A2, scalar or batch input.

    Lab-frame result is R @ K @ R.T. Applying a Gamma contrast to local Kii
    assumes unchanged sum(J); absolute Gamma must not replace local Kii.
    """
    if np.iscomplexobj(strain):
        raise ValueError("Biot strain must be real")
    E = np.asarray(strain, dtype=float)
    if E.shape[-2:] != (3, 3) or not np.isfinite(E).all():
        raise ValueError("Expected finite (...,3,3) strain")
    if np.max(np.abs(E-np.swapaxes(E, -1, -2))) > 1e-10:
        raise ValueError("Biot strain must be symmetric")
    if not extrapolate and not np.all(calibration_mask(E, model)):
        raise ValueError("Outside the calibration box; explicit extrapolate=True required")
    coefficients = np.asarray(model["coefficients_eV_A2"], dtype=float)
    if coefficients.ndim != 1 or len(coefficients) not in (3, 12):
        raise ValueError("Unsupported cubic polynomial basis")
    basis = feature_tensor(E)[..., :len(coefficients), :, :]
    return model["K0_eV_A2"]*np.eye(3)+np.einsum("...kij,k->...ij", basis, coefficients)
