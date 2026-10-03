"""Lightweight video container and readers.

The analysis routines in :mod:`microrheology.ddm` and
:mod:`microrheology.tracking` operate on a :class:`Video`: an object
with ``len()``, integer indexing (frames as numpy arrays), a
``frame_shape`` attribute, and an ``info`` dict for metadata
(fps, microns-per-pixel, temperature, tracer size, ...).

Use :func:`read_video` for mp4 files (OpenCV backend, works in
JupyterLite/Pyodide). For Phantom ``.cin`` files use ``pims.Cine``
directly and wrap it — the analysis code only needs the duck-typed
interface above.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


class Video:
    """Frames + metadata + slots for analysis results.

    Parameters
    ----------
    frames : array-like, shape (n_frames, height, width[, channels])
    fps : frames per second
    muperpix : microns per pixel
    **info : extra metadata (temperature_C, tracer_diameter_um, filename, ...)
    """

    def __init__(self, frames, fps=100.0, muperpix=0.1, **info):
        self._frames = np.asarray(frames)
        if self._frames.ndim not in (3, 4):
            raise ValueError("frames must be (n, h, w) or (n, h, w, c)")
        self.info = {"fps": float(fps), "muperpix": float(muperpix), **info}
        # result slots filled by analysis routines
        self.result_ddm = None
        self.result_tracks = None

    def __len__(self):
        return self._frames.shape[0]

    def __getitem__(self, idx):
        return self._frames[idx]

    @property
    def frame_shape(self):
        """Shape of a single frame (h, w[, c])."""
        return self._frames.shape[1:]

    @property
    def fps(self):
        return self.info["fps"]

    @property
    def interval(self):
        """Seconds between frames."""
        return 1.0 / self.info["fps"]

    def __repr__(self):
        return (f"Video(n_frames={len(self)}, frame_shape={self.frame_shape}, "
                f"fps={self.fps:g}, muperpix={self.info['muperpix']:g})")


def read_video(path, fps=100.0, muperpix=0.1, max_frames=None, to_gray=True):
    """Read a video file into a :class:`Video` (OpenCV backend).

    Works with mp4/avi/mov — and inside JupyterLite/Pyodide, where the
    ``opencv-python`` build is available. For ``.cin`` Phantom files use
    ``pims.Cine`` instead (desktop only).
    """
    try:
        import cv2
    except ImportError as exc:
        raise ImportError("read_video needs opencv-python "
                          "(pip install opencv-python)") from exc

    path = Path(path)
    cap = cv2.VideoCapture(str(path))
    if not cap.isOpened():
        raise OSError(f"cannot open video: {path}")
    file_fps = cap.get(cv2.CAP_PROP_FPS)
    frames = []
    while True:
        ok, frame = cap.read()
        if not ok:
            break
        if to_gray and frame.ndim == 3:
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        frames.append(frame)
        if max_frames is not None and len(frames) >= max_frames:
            break
    cap.release()
    if not frames:
        raise OSError(f"no frames read from {path}")
    return Video(np.asarray(frames), fps=file_fps or fps, muperpix=muperpix,
                 filename=path.name)
