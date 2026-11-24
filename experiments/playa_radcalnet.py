"""Compare playa reflectance with RadCalNet RVUS ground reflectance (2020-07-12, 21:30 UTC).

    python experiments/playa_radcalnet.py TAUDIR RCN_INPUT_FILE [--plot out.png]

RadCalNet's *.input file holds bottom-of-atmosphere (surface) reflectance every 10 nm, 400-2300 nm
(9999 beyond), plus the site's measured water vapour and AOD. DUSTRAN's rho here has no scattering
correction, so only its SHAPE is compared: each spectrum is divided by RadCalNet and the ratio is
judged for smoothness.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import sys as _sys; _sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
from dustran.paths import DATA
sys.path.insert(0, str(Path(__file__).resolve().parent))
import playa_gas_score as P

ap = argparse.ArgumentParser(); ap.add_argument("taudir"); ap.add_argument("rcn"); ap.add_argument("--plot"); ap.add_argument("--cia", action="store_true"); ap.add_argument("--o3", action="store_true")
a = ap.parse_args()
rows = [l.split("\t") for l in open(a.rcn).read().splitlines()]
head = {}
for r in rows:                       # first block only (a second block repeats keys as uncertainties)
    if r and r[0].strip() and not r[0].strip().isdigit():
        head.setdefault(r[0].rstrip(":"), r[1:])
col = head["UTC"].index("21:30")
spec = [(float(r[0]), float(r[1 + col])) for r in rows if r and r[0].strip().isdigit()][:191]
wr, rr = np.array(spec).T; ok = rr < 9
wr, rr = wr[ok], rr[ok]
print(f"RadCalNet 21:30: water {head['WV'][col]} cm, AOD {head['AOD'][col]}, Angstrom {head['Ang'][col]}, P {head['P'][col]} hPa")

S = P.load_scene(); wl = S["centers"]; X = (DATA / "hitran_xsec")
CIA = {X / "O2-O2_2024.cia": ("O2", "O2"), X / "O2-N2_2024.cia": ("O2", "N2")}
O3 = str(DATA / "bremen_o3/serdyuchenkogorshelev5digits.dat")
M = P.Model(a.taudir.split(","), S, cia_files=CIA if a.cia else None,
            o3_file=O3 if a.o3 else None, o3_du=float(head["O3"][col]) if a.o3 else None)   # comma-separated segment dirs
lin = lambda x: 0.389 + 0.063e-3 * (x - 1500)
d = M.rho(2.15, lin(wl))

def to10(y):   # 10 nm boxcar averages at RadCalNet wavelengths, good bands only
    good = (S["bbl"] > 0)
    return np.array([np.nanmean(np.where(good & (abs(wl - c) <= 5), y, np.nan)) for c in wr])
d10 = to10(d)
ok = np.isfinite(d10)
for name, y in (("DUSTRAN (gas only, no scattering)", d10),):
    r = y[ok] / rr[ok]
    rn = r / np.median(r)
    print(f"\n{name}: ratio to RadCalNet, median {np.median(r):.3f}")
    for lo, hi in ((400, 700), (700, 1000), (1000, 1300), (1500, 1750), (2000, 2300), (2250, 2300)):
        m = (wr[ok] >= lo) & (wr[ok] <= hi)
        if m.any():
            print(f"  {lo}-{hi} nm: mean ratio {r[m].mean():.3f}, spread {100*np.std(rn[m]):.2f}%")
print("\nsurface drop 2270 -> 2300 nm:  RadCalNet {:.1f}%  DUSTRAN {:.1f}%".format(
    *[100 * (y[list(wr).index(2300)] / y[list(wr).index(2270)] - 1) for y in (rr, d10)]))

if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    ax[0].plot(wr, rr, "k-", lw=2, label="RadCalNet ground (21:30)")
    ax[0].plot(wr[ok], d10[ok] * np.median(rr[ok] / d10[ok]), label="DUSTRAN gas only (scaled)")
    ax[0].legend(); ax[0].set_ylabel("reflectance")
    ax[1].plot(wr[ok], d10[ok] / rr[ok] / np.median(d10[ok] / rr[ok]), label="DUSTRAN / RadCalNet (normalised)")
    ax[1].axhline(1, c="gray"); ax[1].legend(); ax[1].set_ylabel("ratio"); ax[1].set_xlabel("nm"); ax[1].set_ylim(0.85, 1.15)
    fig.savefig(a.plot, dpi=110, bbox_inches="tight")
