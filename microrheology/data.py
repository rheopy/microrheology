"""Bundled example videos: download-once, cache, reuse.

The two water test videos (400 nm tracers in water, 0.1 µm/px,
100 fps) live in the repo under ``data/``; the full-fidelity
``.cin`` Phantom files stay on the releases page. Videos are never
shipped inside the wheel — :func:`water_video` fetches and caches
them under ``~/.cache/microrheology/``.
"""

from __future__ import annotations

from pathlib import Path
import shutil
import urllib.request

__all__ = ["WATER_LOW", "WATER_HIGH", "water_video", "CIN_URLS"]

REPO_RAW = "https://raw.githubusercontent.com/rheopy/microrheology/main/data"

#: Low tracer concentration — suited to particle tracking.
WATER_LOW = {
    "filename": "400nm_100dil_water_01umpix_100fps_short.mp4",
    "fps": 100.0,
    "muperpix": 0.1,
    "tracer_diameter_um": 0.4,
    "temperature_C": 25.0,
    "description": "400 nm tracers in water, low concentration "
                   "(particle tracking)",
}

#: High tracer concentration — DDM territory (tracking infeasible).
WATER_HIGH = {
    "filename": "400nm_water_01umpix_100fps_short.mp4",
    "fps": 100.0,
    "muperpix": 0.1,
    "tracer_diameter_um": 0.4,
    "temperature_C": 25.0,
    "description": "400 nm tracers in water, high concentration (DDM)",
}

#: Full-fidelity Phantom files (large; desktop only, via pims).
CIN_URLS = {
    "low": "https://github.com/marcocaggioni/microrheology_water_test"
           "/releases/download/1.0/400nm_100dil_water_01umpix_100fps_short.cin",
    "high": "https://github.com/marcocaggioni/microrheology_water_test"
            "/releases/download/1.0/400nm_water_01umpix_100fps_short.cine",
}


def _cache_dir():
    d = Path.home() / ".cache" / "microrheology"
    d.mkdir(parents=True, exist_ok=True)
    return d


def water_video(which="low"):
    """Local path to a water test video, downloading once if needed.

    Parameters
    ----------
    which : "low" (tracking) or "high" (DDM)
    """
    meta = {"low": WATER_LOW, "high": WATER_HIGH}[which]
    dest = _cache_dir() / meta["filename"]
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    # prefer a repo-local copy when running from a checkout
    local = Path(__file__).resolve().parent.parent / "data" / meta["filename"]
    if local.exists():
        shutil.copy(local, dest)
        return dest
    url = f"{REPO_RAW}/{meta['filename']}"
    print(f"downloading {meta['filename']} ...")
    urllib.request.urlretrieve(url, dest)
    return dest
