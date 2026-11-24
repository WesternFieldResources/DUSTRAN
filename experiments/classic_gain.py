"""Smooth vicarious gain for AVIRIS Classic 2019 from a consensus of independent references.

    python experiments/classic_gain.py

References (both at Sentinel's time of day): Sentinel-2A L2A and RadCalNet RVUS. They agree
within their uncertainty, so their mean is the consensus. Gain g(wl) = consensus / DUSTRAN
(DUSTRAN moved to that time of day with RadCalNet's diurnal curve; see three_way.py).
g is fitted as a quadratic in wavelength: smooth, 3 numbers, no per-band freedom.
It is meant for RADIANCE (applied before DUSTRAN); over a bright surface reflectance scales
almost linearly with radiance, so the reflectance-domain ratio is a good first estimate, and the
script re-runs DUSTRAN with the gain to check that, and checks water (not used in the fit).
Writes the gain to classic_2019_gain.npy next to the scene data.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import classic_playa as C

t = np.load(C.SC / "three_way.npy")             # rows: c, s2, rcn, dustran, nasa
c, s2, rcn, dus = t[:, 0], t[:, 1], t[:, 2], t[:, 3]
cons = (s2 + rcn) / 2
g = cons / dus
x = (c - 1500) / 1000
p = np.polyfit(x, g, 2)
fit = np.polyval(p, x)
print("band  nm    S2/RCN  consensus/DUSTRAN  quadratic  residual")
for i in range(len(c)):
    print(f"      {c[i]:6.0f}  {s2[i] / rcn[i]:.3f}   {g[i]:.3f}              {fit[i]:.3f}     {100 * (g[i] - fit[i]):+.1f}%")
print("quadratic coefficients (x = (nm-1500)/1000):", np.round(p, 4))
wl = np.arange(350, 2551, 5.0)
np.save(C.SC / "classic_2019_gain.npy", np.vstack([wl, np.polyval(p, (wl - 1500) / 1000)]))
