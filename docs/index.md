# 🔬 microrheology

**Video in, rheology out.** Two roads from a microscopy movie of tracer
particles to material properties:

| | 📍 Particle tracking | 🌊 DDM |
|---|---|---|
| needs | sparse tracers, `trackpy` | any concentration |
| runs in browser | ✗ (no wheels) | ✓ |
| walkthrough | [tracking](walkthrough-tracking.md) | [ddm](walkthrough-ddm.md) |

```python
import microrheology as mr

video = mr.read_video(mr.data.water_video("high"), fps=100, muperpix=0.1)
mr.ddm.calculate_ddm(video, n_average=20, n_deltas=10)
mr.ddm.viscosity_from_ddm(video, radius_um=0.2)   # ≈ 1e-3 Pa·s for water
```

```{toctree}
:maxdepth: 2

walkthrough-tracking
walkthrough-ddm
api
```
