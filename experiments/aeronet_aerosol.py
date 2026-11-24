"""Aerosol for the playa scene from AERONET Railroad_Valley2 (PI J. Czapla-Myers, P. Gupta).

    from aeronet_aerosol import load
    aer, aod500 = load()

Size distribution and refractive index: the Level 1.5 hybrid inversion at 18:50 UTC on
2020-07-12 (closest to the 21:50 flight). AOD: Level 1.5 direct sun at 21:48:49 UTC, 500 nm.
AERONET notes that SSA / imaginary index are unreliable at AOD440 < 0.4 (here ~0.08); the size
distribution is the robust part. Files live outside the repo in dustran_data/aeronet.
"""
import csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dustran.aerosol import Measured

DIR = Path("D:/workspace/dustran_data/aeronet_2020_194")


def _load(f, date="12:07:2020", d=DIR):
    rows = list(csv.reader(open(Path(d) / f)))
    i = [k for k, r in enumerate(rows) if r and (r[0].startswith("AERONET_Site") or r[0].startswith("Date("))][0]
    return rows[i], [r for r in rows[i + 1:] if r and date in ",".join(r[:3])]


def load(inv_time="18:50:28", sun_time="21:48:49", date="12:07:2020", d=DIR, sun_file="aod.txt"):
    """Aerosol (sizes + index from the inversion at inv_time) and direct-sun AOD + water at sun_time."""
    h, rs = _load("SIZ_hyb.txt", date, d); r = [r for r in rs if r[2] == inv_time][0]
    radii = np.array([float(c) for c in h[5:27]]); dv = np.array([float(c) for c in r[5:27]])
    h, rs = _load("RIN_hyb.txt", date, d); dd = dict(zip(h, [r for r in rs if r[2] == inv_time][0]))
    wls = [440, 675, 870, 1020]
    n = [float(dd[f"Refractive_Index-Real_Part[{w}nm]"]) for w in wls]
    k = [float(dd[f"Refractive_Index-Imaginary_Part[{w}nm]"]) for w in wls]
    h, rs = _load(sun_file, date, d); ti = [i for i, c in enumerate(h) if c.startswith("Time")][0]
    dd = dict(zip(h, [r for r in rs if r[ti] == sun_time][0]))
    meas = {w: float(dd[f"AOD_{w}nm"]) for w in (440, 500, 675, 870, 1020, 1640)}
    meas["PW"] = float(dd["Precipitable_Water(cm)"])
    return Measured(radii, dv, wls, n, k), meas


if __name__ == "__main__":
    aer, meas = load()
    e500 = aer.optics(500)["ext"]
    print("wl   AOD model   AOD measured   ssa    g")
    for w in (440, 500, 675, 870, 1020, 1640):
        o = aer.optics(w)
        print(f"{w:5d} {meas[500] * o['ext'] / e500:9.4f} {meas[w]:12.4f}   {o['ssa']:.3f}  {o['g']:.3f}")
