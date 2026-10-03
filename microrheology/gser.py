"""MSD → viscoelastic moduli via the generalized Stokes–Einstein relation.

Implements T. Mason's algebraic approximation to the GSER: the
mean-squared displacement's local log-log slope and curvature give
``|G*(ω)|``, and a second-order analytic continuation splits it into
``G′(ω)`` and ``G″(ω)``.

References:

- Mason, T. G. & Weitz, D. A., Phys. Rev. Lett. **74**, 1250–1253
  (1995). https://doi.org/10.1103/PhysRevLett.74.1250
- Mason, T. G., Rheol. Acta **39**, 371–378 (2000).
  https://doi.org/10.1007/s003970000094
- Mason, T. G. et al., Phys. Rev. Lett. **79**, 3284–3287 (1997)
  (second-order corrections).

This is an independent implementation of the published method.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy import constants
from scipy.special import gamma

__all__ = ["msd_to_moduli"]


def _local_log_derivatives(x, y, width):
    """Local slope and curvature of log(y) vs log(x).

    At each point, fits a quadratic in log-log space weighted by a
    Gaussian of the given width; returns the smoothed value and its
    first two logarithmic derivatives.
    """
    lx = np.log(x)
    ly = np.log(y)
    smooth = np.empty_like(y)
    slope = np.empty_like(y)
    curve = np.empty_like(y)
    for i, u in enumerate(lx):
        w = np.exp(-((lx - u) ** 2) / (2 * width ** 2))
        sel = w > 0.03
        coef = np.polyfit(lx[sel], ly[sel], 2, w=w[sel])
        p = np.poly1d(coef)
        smooth[i] = np.exp(p(u))
        slope[i] = p.deriv(1)(u)
        curve[i] = p.deriv(2)(u)
    return smooth, slope, curve


def msd_to_moduli(tau, msd, radius_um, temperature_K, dim=2, clip=0.03,
                  smooth_width=0.7):
    """Convert MSD(τ) to the complex modulus G*(ω).

    Parameters
    ----------
    tau : array, seconds — lag times
    msd : array, µm² — mean-squared displacement
    radius_um : tracer radius in microns
    temperature_K : absolute temperature
    dim : dimensionality of the MSD data (2 for video tracking)
    clip : moduli below ``clip`` × ``|G*|`` are set to zero (unreliable)
    smooth_width : Gaussian width for the local log-derivative fits;
        increase for noisy data

    Returns
    -------
    omega : array, rad/s
    G_star : complex array, Pa (``.real`` = G′, ``.imag`` = G″)

    Needs ~8 points per decade; warns if the local curvature is high.
    """
    tau = np.asarray(tau, dtype=float)
    msd = np.asarray(msd, dtype=float)
    if tau.shape != msd.shape:
        raise ValueError("tau and msd must have the same shape")
    a_m = radius_um * 1e-6
    msd_m2 = msd * 1e-12
    omega = 1.0 / tau

    prefactor = dim * constants.k * temperature_K / (3 * np.pi * a_m)

    # |G*(ω)| from the MSD's local slope α and curvature β (2nd order)
    msd_s, alpha, beta = _local_log_derivatives(tau, msd_m2, smooth_width)
    g_abs = prefactor / (msd_s * gamma(1 + alpha) * (1 + beta / 2))

    # split into G′ and G″ via the 2nd-order continuation in ω
    g_s, alpha_w, beta_w = _local_log_derivatives(omega, g_abs, smooth_width)
    correction = (np.pi / 2 - 1) * (alpha_w + 1j * (1 - alpha_w)) * beta_w
    g_star = (g_s / (1 + beta_w)
              * (np.exp(0.5j * np.pi * alpha_w) - correction))

    if np.abs(beta).max() > 0.15 or np.abs(beta_w).max() > 0.15:
        warnings.warn("High curvature in MSD data; moduli may be unreliable.")

    g_star = g_star.copy()
    g_star[g_star.real < np.abs(g_star) * clip] = \
        1j * g_star[g_star.real < np.abs(g_star) * clip].imag
    g_star[g_star.imag < np.abs(g_star) * clip] = \
        g_star[g_star.imag < np.abs(g_star) * clip].real + 0j
    return omega, g_star
