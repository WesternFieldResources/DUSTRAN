"""Wavelength shift on the playa: best shift per gas window, then a straight-line shift.

    python experiments/playa_shift.py TAUDIR

Uses the playa_gas_score model (line-by-line gases, TSIS, no scattering). Per window, the shift
that minimises roughness (baseline-free: rms of the 2nd difference, so the smoothing spline
cannot be bent by the shift). If the shift grows with wavelength, that is a grating/dispersion
error, not a constant offset. Then fits shift(wl) = a + b*(wl - 1500 nm) with water, and rescores.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, str(Path(__file__).resolve().parent))
import playa_gas_score as P

S = P.load_scene(); wl = S["centers"]; M = P.Model(sys.argv[1], S)
W = 2.15
shifts = np.round(np.arange(-1.0, 1.01, 0.05), 2)
print("best shift per window (min roughness), w = 2.15:")
per = {}
for g in P.GAS:
    r = [P.score(wl, M.rho(W, s), S["bbl"])[0][g][1] for s in shifts]
    i = int(np.argmin(r)); per[g] = shifts[i]
    print(f"  {g:16s} {shifts[i]:+.2f} nm   roughness {r[i]:.2f} (at 0: {r[list(shifts).index(0.0)]:.2f})")
mids = np.array([np.mean(v) for v in P.GAS.values()]); vals = np.array(list(per.values()))
b, a = np.polyfit(mids - 1500, vals, 1)
print(f"line through per-window shifts: {a:+.3f} nm at 1500 nm, slope {b*1000:+.3f} nm per 1000 nm")

def cost(x):
    w, a, b = x
    sc = P.score(wl, M.rho(w, a + b * (wl - 1500)), S["bbl"])[0]
    return np.sqrt(np.mean([sc[g][1] ** 2 for g in P.GAS]))
res = minimize(cost, [W, a, b], method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-3, maxiter=150))
w, a, b = res.x
sc = P.score(wl, M.rho(w, a + b * (wl - 1500)), S["bbl"])[0]
print(f"\njoint fit (roughness): w={w:.3f}, shift = {a:+.3f} nm at 1500, slope {b*1000:+.3f} nm/1000 nm"
      f"  -> {a+b*(760-1500):+.2f} nm at 760, {a+b*(2330-1500):+.2f} nm at 2330")
print(f"{'window':16s} {'line shift':>14s}")
print(f"{'clean bands':16s} {sc['clean']:14.2f}")
for g in P.GAS:
    print(f"{g:16s} {sc[g][0]:6.2f} / {sc[g][1]:5.2f}")
