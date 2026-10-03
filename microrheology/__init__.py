"""microrheology — video in, rheology out.

Two roads from a microscopy movie to material properties:

- :mod:`microrheology.tracking` — particle-tracking microrheology
  (locate → link → MSD → viscosity / moduli)
- :mod:`microrheology.ddm` — differential dynamic microscopy
  (image structure function → diffusion → viscosity / radius)

Shared pieces:

- :mod:`microrheology.video` — lightweight video container + reading
- :mod:`microrheology.gser` — MSD → G*(ω) via the generalized
  Stokes–Einstein relation
- :mod:`microrheology.data` — bundled example videos
"""

from microrheology.video import Video, read_video
from microrheology import data, ddm, gser, tracking, video

__version__ = "0.1.0"
__all__ = ["Video", "read_video", "data", "ddm", "gser", "tracking", "video",
           "__version__"]
