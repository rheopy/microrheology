"""Smoke tests: imports, GSER on synthetic diffusive MSD, DDM on synthetic video."""

import numpy as np


def test_imports():
    import microrheology as mr
    assert mr.__version__ == "0.1.0"
    for mod in ["video", "ddm", "gser", "tracking", "data"]:
        assert hasattr(mr, mod)


def test_gser_newtonian():
    """Diffusive MSD -> G'' ~ omega, G' ~ 0 (water-like)."""
    from microrheology import gser
    tau = np.logspace(-2, 1, 60)
    D = 2.0  # um^2/s
    msd = 4 * D * tau  # 2D diffusive
    omega, Gs = gser.msd_to_moduli(tau, msd, radius_um=0.2,
                                   temperature_K=298.15, dim=2)
    assert omega.shape == Gs.shape == tau.shape
    # Newtonian: loss dominates, grows ~linearly with omega
    assert np.median(Gs.imag) > np.median(np.abs(Gs.real))
    slope = np.polyfit(np.log(omega), np.log(Gs.imag + 1e-30), 1)[0]
    assert 0.8 < slope < 1.2, slope


def test_azimuthal_average():
    from microrheology.ddm import azimuthal_average
    y, x = np.mgrid[-50:51, -50:51]
    img = np.exp(-(x ** 2 + y ** 2) / 200.0)
    prof = azimuthal_average(img)
    assert prof[0] > prof[-1]  # decreasing radial profile
    assert np.all(np.isfinite(prof))


def test_ddm_viscosity_synthetic():
    """DDM on a synthetic Brownian movie recovers a sane viscosity."""
    from microrheology import ddm
    from microrheology.video import Video
    rng = np.random.default_rng(0)
    n, h, w = 60, 64, 64
    # random walkers as bright dots
    frames = np.zeros((n, h, w))
    pos = rng.uniform(0, min(h, w), size=(40, 2))
    for i in range(n):
        pos += rng.normal(0, 1.2, size=pos.shape)
        pos %= min(h, w)
        for (yy, xx) in pos.astype(int):
            frames[i, yy, xx] = 255
    video = Video(frames, fps=100, muperpix=0.1)
    ddm.calculate_ddm(video, n_average=8, n_deltas=6, show_progress=False)
    assert video.result_ddm is not None
    eta = ddm.viscosity_from_ddm(video, radius_um=0.2, qmin=1, qmax=8)
    assert np.isfinite(eta) and eta > 0


def test_video_container():
    from microrheology.video import Video
    v = Video(np.zeros((10, 32, 32)), fps=100, muperpix=0.1)
    assert len(v) == 10
    assert v.frame_shape == (32, 32)
    assert abs(v.interval - 0.01) < 1e-12
