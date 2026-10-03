# 🔬 microrheology

**Video in, rheology out.** Particle-tracking microrheology and differential
dynamic microscopy (DDM) for extracting diffusion, viscosity, and viscoelastic
moduli from microscopy movies of tracer particles.

```python
import microrheology as mr

# --- DDM: works anywhere, including in the browser ---
video = mr.read_video(mr.data.water_video("high"), fps=100, muperpix=0.1)
mr.ddm.calculate_ddm(video, n_average=20, n_deltas=10)
eta = mr.ddm.viscosity_from_ddm(video, radius_um=0.2)   # Pa·s, water ≈ 1e-3

# --- particle tracking (desktop/Colab; needs trackpy) ---
video = mr.read_video(mr.data.water_video("low"), fps=100, muperpix=0.1)
features = mr.tracking.batch_locate(video, diameter=15, minmass=1000)
traj = mr.tracking.link_trajectories(features)
msd = mr.tracking.mean_squared_displacement(traj, muperpix=0.1, fps=100)
omega, Gstar = mr.gser.msd_to_moduli(msd.index, msd.values,
                                     radius_um=0.2, temperature_K=298.15)
```

## Install

```bash
pip install rheopy-microrheology        # core (DDM without fitting: numpy/scipy/pandas/matplotlib)
pip install "rheopy-microrheology[ddm]"      # + lmfit for DDM fits
pip install "rheopy-microrheology[tracking]" # + trackpy, pims for particle tracking
pip install "rheopy-microrheology[all]"      # everything
```

## Example data

Two water test videos (400 nm tracers, 0.1 µm/px, 100 fps) live in
[`data/`](data/) and are fetched on demand by `microrheology.data`
— never shipped inside the wheel:

- **low concentration** → particle tracking
- **high concentration** → DDM (tracking infeasible)

Full-fidelity Phantom `.cin` files are linked from `microrheology.data.CIN_URLS`.

## Docs

Two walkthroughs follow the notebook workflow:

1. 📍 [Particle-tracking microrheology](https://rheopy.github.io/microrheology/walkthrough-tracking.html)
2. 🌊 [Differential dynamic microscopy](https://rheopy.github.io/microrheology/walkthrough-ddm.html) *(browser-friendly)*

## Citations

- DDM: Cerbino, R. & Trappe, V., Phys. Rev. Lett. **100**, 188102 (2008).
  [doi:10.1103/PhysRevLett.100.188102](https://doi.org/10.1103/PhysRevLett.100.188102)
- GSER: Mason, T. G. & Weitz, D. A., Phys. Rev. Lett. **74**, 1250–1253 (1995).
  [doi:10.1103/PhysRevLett.74.1250](https://doi.org/10.1103/PhysRevLett.74.1250);
  Mason, T. G., Rheol. Acta **39**, 371–378 (2000).
  [doi:10.1007/s003970000094](https://doi.org/10.1007/s003970000094)
- Particle tracking: Crocker, J. C. & Grier, D. G., J. Colloid Interface Sci.
  **179**, 298–310 (1996).
  [doi:10.1006/jcis.1996.0217](https://doi.org/10.1006/jcis.1996.0217)

## License

MIT. The GSER implementation in `microrheology/gser.py` is an independent
implementation of Mason's published method.
