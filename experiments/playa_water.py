"""Water column on the playa, fitted against RadCalNet's measured surface shape.

    python experiments/playa_water.py TAUDIR RCN_INPUT_FILE

The smooth-spline baseline used so far assumes the surface is featureless across each water
window. Here the surface shape comes from RadCalNet's ground reflectance instead (10 nm, linearly
interpolated to band centres). For each water window, the ratio model/ground should be flat
(scattering is smooth); a straight line is fitted to the ratio on clean bands just outside the
window, and water is chosen to minimise the ratio's rms departure from that line inside it.
Reported in cm (= g cm^-2) using the AFGL US-standard column above the playa.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import playa_gas_score as P
import playa_gas_tau as T
from dustran import atmosphere

taudir, rcn = sys.argv[1], sys.argv[2]
rows = [l.split("\t") for l in open(rcn).read().splitlines()]
head = {}
for r in rows:
    if r and r[0].strip() and not r[0].strip().isdigit():
        head.setdefault(r[0].rstrip(":"), r[1:])
col = head["UTC"].index("21:30")
spec = np.array([(float(r[0]), float(r[1 + col])) for r in rows if r and r[0].strip().isdigit()][:191])
spec = spec[spec[:, 1] < 9]
L0 = atmosphere.layers(atmosphere.load_afgl(T.AFGL),
                       surface_km=P.Model.__init__.__globals__.get("GROUND_KM", 1.4357627632017484))
cm1 = L0["col"]["H2O"].sum() / 3.3428e22                    # AFGL column above the playa, g cm^-2
print(f"AFGL US std water above playa: {cm1:.3f} cm; RadCalNet photometer: {head['WV'][col]} cm "
      f"(= w {float(head['WV'][col]) / cm1:.2f})")

S = P.load_scene(); wl = S["centers"]; M = P.Model(taudir, S)
ground = np.interp(wl, spec[:, 0], spec[:, 1])
shift = 0.389 + 0.063e-3 * (wl - 1500)
good = S["bbl"] > 0
WIN = {"H2O 940": ((880, 1000), (850, 880), (1000, 1040)), "H2O 1140": ((1080, 1180), (1040, 1075), (1230, 1260))}
ws = np.round(np.arange(1.0, 3.01, 0.02), 3)
for name, ((lo, hi), e1, e2) in WIN.items():
    inw = good & (wl >= lo) & (wl <= hi)
    edge = good & (((wl >= e1[0]) & (wl <= e1[1])) | ((wl >= e2[0]) & (wl <= e2[1])))
    costs = []
    for w in ws:
        q = M.rho(w, shift) / ground
        b = np.polyfit(wl[edge], q[edge], 1)
        costs.append(100 * np.sqrt(np.mean((q[inw] / np.polyval(b, wl[inw]) - 1) ** 2)))
    i = int(np.argmin(costs))
    c_rcn = costs[int(np.argmin(abs(ws - float(head['WV'][col]) / cm1)))]
    print(f"{name}: best w {ws[i]:.2f} = {ws[i] * cm1:.2f} cm (rms {costs[i]:.2f}%); "
          f"at RadCalNet water the rms is {c_rcn:.2f}%")
