"""Particle-tracking microrheology.

Thin, documented wrappers around ``trackpy`` implementing the classic
Crocker–Grier workflow: locate → link → mean-squared displacement →
viscosity / moduli.

Reference: Crocker, J. C. & Grier, D. G., J. Colloid Interface Sci.
**179**, 298–310 (1996). https://doi.org/10.1006/jcis.1996.0217

``trackpy`` is imported lazily so ``import microrheology`` stays light;
these functions raise a clear error if it is missing
(``pip install trackpy``). Note: trackpy ships no wheels, so this path
is desktop/Colab only — the DDM path in :mod:`microrheology.ddm` is the
browser-friendly one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "locate_frame",
    "batch_locate",
    "link_trajectories",
    "mean_squared_displacement",
    "viscosity_from_msd",
]


def _trackpy():
    try:
        import trackpy as tp
    except ImportError as exc:
        raise ImportError("particle tracking needs trackpy "
                          "(pip install trackpy)") from exc
    return tp


def locate_frame(frame, diameter=15, minmass=1000, invert=True):
    """Locate particles in a single frame (Crocker–Grier).

    Returns a DataFrame with ``x``, ``y`` (px) and ``mass`` columns.
    """
    tp = _trackpy()
    return tp.locate(np.asarray(frame), diameter, minmass=minmass,
                     invert=invert)


def batch_locate(video, diameter=15, minmass=1000, invert=True,
                 start_frame=0, end_frame=None):
    """Locate particles in every frame; returns a stacked DataFrame."""
    tp = _trackpy()
    n = len(video)
    if end_frame is None:
        end_frame = n
    parts = []
    for i in range(start_frame, end_frame):
        df = tp.locate(np.asarray(video[i]), diameter, minmass=minmass,
                       invert=invert)
        df["frame"] = i
        parts.append(df)
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    video.result_tracks = out
    return out


def link_trajectories(features, search_range=5, memory=3, min_length=50):
    """Link per-frame features into trajectories; drop short stubs."""
    tp = _trackpy()
    t = tp.link_df(features, search_range, memory=memory)
    return tp.filter_stubs(t, min_length)


def mean_squared_displacement(traj, muperpix, fps, max_lagtime=100,
                              ensemble=True):
    """Mean-squared displacement (µm²) vs lag time (s).

    Uses the ensemble MSD by default (``trackpy.emsd``); set
    ``ensemble=False`` for per-particle MSDs (``trackpy.imsd``).
    Returns a pandas Series indexed by lag time in seconds.
    """
    tp = _trackpy()
    fn = tp.emsd if ensemble else tp.imsd
    msd = fn(traj, muperpix, fps, max_lagtime=max_lagtime)
    msd.index = msd.index  # already in seconds
    return msd


def viscosity_from_msd(tau, msd, temperature_C=25.0, radius_um=0.2, dim=2):
    """Viscosity (Pa·s) from a Newtonian MSD via Stokes–Einstein.

    Fits ``MSD = 2·dim·D·τ`` in log space; ``D`` then gives ``η``.
    For water-like tracers expect ``n ≈ 1`` (diffusive).
    """
    from scipy import constants
    tau = np.asarray(tau, dtype=float)
    msd = np.asarray(msd, dtype=float)
    mask = (tau > 0) & (msd > 0) & np.isfinite(msd)
    slope, intercept = np.polyfit(np.log(tau[mask]), np.log(msd[mask]), 1)
    D_um2_s = np.exp(intercept) / (2 * dim)          # µm²/s
    eta = (constants.k * (273.15 + temperature_C)
           / (6 * np.pi * D_um2_s * radius_um * 1e-18))
    return {"eta_Pa_s": float(eta), "D_um2_s": float(D_um2_s),
            "loglog_slope": float(slope)}
