"""Derive the conservative default bad-band ranges (dustran.batch.CORE_OPAQUE_NM).

    python experiments/core_opaque.py TAUDIR [--w 0.3] [--threshold 0.01]

Bands are "no data, never data" where almost no sunlight survives the two-way path even in dry
air. Two-way water transmittance (sun and view at zenith) from DUSTRAN's line-by-line optical
depths, water at w x AFGL US standard, smoothed with Gaussian bands of 5 nm (AVIRIS-NG) and
10 nm (AVIRIS Classic) FWHM. A wavelength is opaque if the band-averaged transmittance is below
the threshold at both widths; contiguous runs wider than 20 nm are reported.
"""
import argparse
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument("taudir")
ap.add_argument("--w", type=float, default=0.3); ap.add_argument("--threshold", type=float, default=0.01)
a = ap.parse_args()
z = np.load(Path(a.taudir) / "tau_H2O.npz")
wn, tau = z["wn"], (z["all"] + z["below"]).astype(np.float64) * a.w
wl = 1e7 / wn
grid = np.arange(1200, 2100, 1.0)
opaque = np.ones(len(grid), bool)
for fwhm in (5.0, 10.0):
    s = fwhm / 2.3548
    T = np.exp(-tau)
    out = np.empty(len(grid))
    for i, c in enumerate(grid):
        m = np.abs(wl - c) < 4 * s
        g = np.exp(-0.5 * ((wl[m] - c) / s) ** 2)
        out[i] = np.sum(g * T[m]) / np.sum(g)
    opaque &= out < a.threshold
runs, start = [], None
for x, o in zip(grid, opaque):
    if o and start is None:
        start = x
    if not o and start is not None:
        runs.append((start, x - 1)); start = None
runs = [r for r in runs if r[1] - r[0] >= 20]
print(f"water {a.w} x AFGL, threshold {a.threshold}: opaque ranges (nm):", runs)
