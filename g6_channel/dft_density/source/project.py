"""Project the canonical real trigonometric interpolant, not a modulo-N diagonal."""
import numpy as np


def split_nyquist(density):
    if density.ndim != 3 or len(set(density.shape)) != 1 or density.shape[0] % 2:
        raise ValueError('Expected an even cubic Native density grid')
    if not np.isfinite(density).all():
        raise ValueError('Nonfinite density')
    n = density.shape[0]
    coefficients = np.fft.fftshift(np.fft.fftn(density)) / density.size
    for axis in range(3):
        shape = list(coefficients.shape)
        shape[axis] += 1
        expanded = np.zeros(shape, dtype=complex)
        middle = [slice(None)] * 3
        middle[axis] = slice(1, n)
        expanded[tuple(middle)] = coefficients[tuple(middle)]
        first = [slice(None)] * 3
        last = [slice(None)] * 3
        first[axis], last[axis] = 0, n
        expanded[tuple(first)] = coefficients[tuple(first)] / 2
        expanded[tuple(last)] = coefficients[tuple(first)] / 2
        coefficients = expanded
    return coefficients, np.arange(-n // 2, n // 2 + 1)


def project(density, grid=128):
    coefficients, frequencies = split_nyquist(density)
    if grid % 2 or 4 * max(frequencies) >= grid:
        raise ValueError('Output grid would alias projected modes')
    projected = np.zeros((grid, grid), complex)
    without_nyquist = np.zeros_like(projected)
    supported = np.zeros((grid, grid), bool)
    bound = int(max(frequencies))
    for ix, h in enumerate(frequencies):
        for iz, l in enumerate(frequencies):
            index = (int(2 * h) % grid, int(l) % grid)
            value = coefficients[ix, bound - h, iz]
            projected[index] += value
            supported[index] = True
            if abs(h) < bound and abs(l) < bound:
                without_nyquist[index] += value
    image = np.fft.ifft2(projected) * grid**2
    if np.max(np.abs(image.imag)) > 1e-11:
        raise ValueError('Non-Hermitian projected density')
    if abs(float(image.real.mean()) - float(density.mean())) > 1e-12:
        raise ValueError('Projection changed the density zero mode')
    return image.real, projected, supported, without_nyquist


def evaluate(coefficients, fractions):
    """Evaluate exactly at transverse fractional points; no nearest-pixel sampling."""
    q = np.fft.fftfreq(coefficients.shape[0]) * coefficients.shape[0]
    phase_x = np.exp(2j * np.pi * np.asarray(fractions)[:, 0, None] * q)
    phase_y = np.exp(2j * np.pi * np.asarray(fractions)[:, 1, None] * q)
    values = np.einsum('ij,pi,pj->p', coefficients, phase_x, phase_y)
    if np.max(np.abs(values.imag)) > 1e-11:
        raise ValueError('Nonreal point evaluation')
    return values.real
