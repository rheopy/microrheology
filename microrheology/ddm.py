"""Differential dynamic microscopy (DDM).

DDM turns a plain bright-field microscopy movie into a scattering
experiment: for each lag time Δt, the power spectrum of image
differences is computed and azimuthally averaged to give the image
structure function ``D(q, Δt)``. For Brownian particles

.. math:: D(q, \\Delta t) = A(q)\\,[1 - e^{-D q^2 \\Delta t}] + B(q)

so fitting each q-slice gives a relaxation time ``τ(q) = 1/(D q²)``;
a power-law fit ``τ = 1/(D q²)`` over q yields the diffusion
coefficient, and Stokes–Einstein gives viscosity (or radius).

Reference: Cerbino, R. & Trappe, V., Phys. Rev. Lett. **100**, 188102
(2008). https://doi.org/10.1103/PhysRevLett.100.188102

The ``video`` argument is any object with ``len()``, integer frame
indexing and ``frame_shape`` — e.g. :class:`microrheology.video.Video`
or a ``pims`` reader.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import fftpack

__all__ = [
    "azimuthal_average",
    "structure_function",
    "calculate_ddm",
    "ddm_decay",
    "viscosity_from_ddm",
    "radius_from_ddm",
    "browse_ddm",
    "explore_ddm",
]


def azimuthal_average(image, center=None):
    """Azimuthally averaged radial profile of a 2D image.

    Bins pixels by integer radius from ``center`` (defaults to the
    image center) and returns the mean value per bin.
    """
    y, x = np.indices(image.shape)
    if center is None:
        center = np.array([(x.max() - x.min()) / 2.0,
                           (x.max() - x.min()) / 2.0])
    r = np.hypot(x - center[0], y - center[1])
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]
    r_int = r_sorted.astype(int)
    deltar = r_int[1:] - r_int[:-1]
    rind = np.where(deltar)[0]
    nr = rind[1:] - rind[:-1]
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]
    return tbin / nr


def _as_gray(video, frame):
    img = np.asarray(video[frame])
    if img.ndim == 3:
        img = img[:, :, 1]  # green channel
    return img.astype(float)


def _square_crop(img):
    n = min(img.shape[:2])
    return img[:n, :n]


def structure_function(video, delta, n_average=None, start_frame=0,
                       end_frame=None, progress=None):
    """Mean azimuthally-averaged power spectrum of frame differences.

    Averages the ``|FFT(I(t+Δt) − I(t))|²`` radial profile over
    ``n_average`` frame pairs separated by ``delta`` frames.
    """
    n_frames = len(video)
    if end_frame is None:
        end_frame = n_frames
    if n_average is None:
        n_average = n_frames - delta
    pairs = range(start_frame, end_frame - delta,
                  max(1, int((end_frame - delta) / n_average)))
    spectra = []
    for f in pairs:
        im1 = _square_crop(_as_gray(video, f))
        im2 = _square_crop(_as_gray(video, f + delta))
        psd2d = np.abs(fftpack.fftshift(fftpack.fft2(im2 - im1))) ** 2
        spectra.append(azimuthal_average(psd2d))
        if progress is not None:
            progress(f, delta)
    return np.mean(spectra, axis=0)


def calculate_ddm(video, n_average=10, n_deltas=20, interval=None,
                  muperpix=None, start_frame=0, end_frame=None,
                  show_progress=True):
    """Compute the DDM matrix ``D(q, Δt)`` for a video.

    Returns a :class:`pandas.DataFrame` indexed by lag time (s) with
    wavevector ``q`` (µm⁻¹) columns; also stored as
    ``video.result_ddm``. ``interval``/``muperpix`` default to the
    video metadata.
    """
    if interval is None:
        interval = video.info.get("interval", 1.0)
    if muperpix is None:
        muperpix = video.info.get("muperpix", 1.0)
    n_frames = len(video)
    if end_frame is None:
        end_frame = n_frames
    span = end_frame - start_frame

    deltas = np.unique(
        np.logspace(0, np.log10(span - int(span / 2)),
                    num=n_deltas, dtype=int))
    lag_times = deltas * interval

    progress = None
    if show_progress:
        try:
            from IPython.display import display
            import ipywidgets as widgets
            w_dt = widgets.Text(description="Δt (frames)", value="")
            w_f = widgets.Text(description="frame pair", value="")
            display(w_dt, w_f)

            def progress(f, d, _w_dt=w_dt, _w_f=w_f):
                _w_dt.value = str(d)
                _w_f.value = str(f)
        except ImportError:
            pass

    rows = [structure_function(video, int(d), n_average,
                               start_frame, end_frame, progress)
            for d in deltas]
    matrix = np.array(rows)
    frame_len = min(video.frame_shape[:2])
    q = np.arange(matrix.shape[1]) * 2 * np.pi / (muperpix * frame_len)
    result = pd.DataFrame(matrix, index=lag_times, columns=q)
    video.result_ddm = result
    return result


def ddm_decay(x, off, amp, tau):
    """Single-exponential DDM relaxation for Brownian particles."""
    return off + amp * (1.0 - np.exp(-x / tau))


def _fit_tau_q(result, qmin=1, qmax=None):
    """Fit D(q,Δt) at each q → relaxation times τ(q)."""
    try:
        from lmfit import Model
        import lmfit
    except ImportError as exc:
        raise ImportError("DDM fitting needs lmfit (pip install lmfit)") from exc
    model = Model(ddm_decay)
    lags = np.asarray(result.index, dtype=float)
    values = result.to_numpy()
    if qmax is None:
        qmax = values.shape[1]
    taus, qs = [], []
    for i in range(qmin, qmax):
        y = values[1:, i]
        try:
            fit = model.fit(y, x=lags[1:],
                            amp=float(np.max(y) - np.min(y)),
                            off=float(np.min(y)),
                            tau=float(np.mean(lags[1:])))
        except Exception:
            continue
        taus.append(fit.params["tau"].value)
        qs.append(float(result.columns[i]))
    return np.asarray(qs), np.asarray(taus)


def _diffusion_from_tau(q, tau, plot=False):
    """τ(q) = 1/(D q²): power-law fit → diffusion coefficient D (µm²/s)."""
    try:
        import lmfit
    except ImportError as exc:
        raise ImportError("DDM fitting needs lmfit (pip install lmfit)") from exc
    model = lmfit.models.PowerLawModel()
    pars = model.guess(tau, x=q)
    out = model.fit(tau, params=pars, x=q)
    if plot:
        import matplotlib.pyplot as plt
        plt.plot(q, tau, "o")
        plt.plot(q, out.best_fit)
        plt.xlabel("q [µm⁻¹]")
        plt.ylabel("τ [s]")
    return 1.0 / out.params["amplitude"].value


def viscosity_from_ddm(video, radius_um, muperpix=None, qmin=1, qmax=None,
                       temperature_C=25.0, plot=False):
    """Viscosity (Pa·s) from a DDM-analyzed video.

    Needs ``video.result_ddm`` (see :func:`calculate_ddm`).
    ``radius_um`` is the tracer radius in microns. Assumes dilute
    Brownian tracers via Stokes–Einstein.
    """
    from scipy import constants
    result = video.result_ddm
    if result is None:
        raise ValueError("run calculate_ddm(video) first")
    q, tau = _fit_tau_q(result, qmin, qmax)
    D = _diffusion_from_tau(q, tau, plot=plot)          # µm²/s
    eta = (constants.k * (273.15 + temperature_C)
           / (6 * np.pi * D * radius_um * 1e-18))
    return float(eta)


def radius_from_ddm(video, viscosity, muperpix=None, qmin=1, qmax=None,
                    temperature_C=25.0, plot=False):
    """Tracer radius (µm) from a DDM-analyzed video of known viscosity."""
    from scipy import constants
    result = video.result_ddm
    if result is None:
        raise ValueError("run calculate_ddm(video) first")
    q, tau = _fit_tau_q(result, qmin, qmax)
    D = _diffusion_from_tau(q, tau, plot=plot)          # µm²/s
    a = (constants.k * (273.15 + temperature_C)
         / (6 * np.pi * D * viscosity * 1e-18))
    return float(a)


def browse_ddm(video, interval=None, start_frame=0, end_frame=None):
    """Interactive widget: frame pair → difference → FFT → radial profile."""
    from ipywidgets import interactive
    import matplotlib.pyplot as plt
    if interval is None:
        interval = video.info.get("interval", 1.0)
    n_frames = len(video)
    if end_frame is None:
        end_frame = n_frames

    def view_image(framenum, delta):
        im1 = _square_crop(_as_gray(video, framenum))
        im2 = _square_crop(_as_gray(video, framenum + delta))
        imdiff = im2 - im1
        psd2d = np.abs(fftpack.fftshift(fftpack.fft2(imdiff))) ** 2
        psd1d = azimuthal_average(psd2d)
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 3, 1)
        plt.title(f"(a) I(r,t), t={framenum * interval:.2f} s")
        plt.imshow(im1, cmap="gray"); plt.axis("off")
        plt.subplot(2, 3, 2)
        plt.title(f"(b) I(r,t+Δt), t={(framenum + delta) * interval:.2f} s")
        plt.imshow(im2, cmap="gray"); plt.axis("off")
        plt.subplot(2, 3, 3)
        plt.title(f"(c) difference, Δt={delta * interval:.2f} s")
        plt.imshow(imdiff, cmap="gray"); plt.axis("off")
        plt.subplot(2, 3, 4)
        plt.title("(d) FFT power")
        plt.imshow(np.log10(psd2d + 1e-12))
        plt.subplot(2, 3, 5)
        plt.plot(psd1d); plt.yscale("log")
        plt.xlabel("spatial frequency"); plt.ylabel("power")
        plt.subplot(2, 3, 6)
        plt.plot(psd1d); plt.xscale("log"); plt.yscale("log")
        plt.xlabel("spatial frequency"); plt.ylabel("power")
        plt.show()

    return interactive(view_image,
                       framenum=(start_frame, max(start_frame, end_frame - 101)),
                       delta=(1, max(2, end_frame - start_frame - 100)))


def explore_ddm(video, interval=None):
    """Interactive widget: inspect D(q,Δt) slices and per-q fits."""
    from ipywidgets import interactive
    import matplotlib.pyplot as plt
    try:
        from lmfit import Model
    except ImportError as exc:
        raise ImportError("explore_ddm needs lmfit (pip install lmfit)") from exc
    result = video.result_ddm
    if result is None:
        raise ValueError("run calculate_ddm(video) first")
    if interval is None:
        interval = video.info.get("interval", 1.0)
    lags = np.asarray(result.index, dtype=float)
    values = result.to_numpy()
    q = np.asarray(result.columns, dtype=float)
    model = Model(ddm_decay)

    def view_plot(selectq, deltat):
        plt.figure(figsize=(15, 3))
        plt.subplot(1, 3, 1)
        plt.plot(q, values[deltat, :], "o")
        plt.xscale("log"); plt.yscale("log")
        plt.axvline(q[selectq], color="k", linestyle="--")
        plt.xlabel("q [µm⁻¹]"); plt.ylabel("I(q) [a.u.]")
        plt.title(r"I(q) at fixed $\Delta t$")
        plt.subplot(1, 3, 2)
        plt.plot(lags, values[:, selectq], "o")
        plt.xscale("log"); plt.yscale("log")
        plt.xlabel("Δt [s]"); plt.ylabel("g(q,Δt) [a.u.]")
        y = values[:, selectq]
        fit = model.fit(y, x=lags,
                        amp=float(np.max(y[1:]) - np.min(y[1:])),
                        off=float(np.min(y[1:])),
                        tau=float(np.mean(lags[1:])))
        tau = fit.params["tau"].value
        off = fit.params["off"].value
        amp = fit.params["amp"].value
        plt.plot(lags, fit.best_fit, "r-")
        plt.subplot(1, 3, 3)
        plt.plot(lags, -((y - (off + amp)) / amp), "o")
        plt.plot(lags, np.exp(-lags / tau))
        plt.xscale("log")
        plt.xlabel("Δt [s]"); plt.ylabel("f(q,Δt) [a.u.]")
        plt.show()

    return interactive(view_plot,
                       selectq=(1, values.shape[1] - 1),
                       deltat=(1, len(lags) - 1))
