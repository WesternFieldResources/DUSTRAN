"""Does the assumed water profile shape change the fitted water column?

    python experiments/water_profile.py BASE_TAUDIR OUT_TAUDIR --low 1.6 --high 0.6 --zsplit 3.0

Builds a variant AFGL profile with water multiplied by `low` below `zsplit` km and by `high`
above, renormalised to the same total column above the playa, recomputes H2O line by line over
850-1270 nm on BASE's exact grid, and writes OUT_TAUDIR = BASE with that H2O section replaced.
Then run playa_water.py on OUT_TAUDIR: if the fitted column moves a lot, profile shape matters.
"""
import argparse, shutil, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dustran import atmosphere, gas
import playa_gas_tau as T

def variant(low, high, zsplit):
    prof = atmosphere.load_afgl(T.AFGL)
    base = atmosphere.layers(prof, T.GROUND_KM)["col"]["H2O"].sum()
    p = dict(prof); f = np.where(prof["z"] < zsplit, low, high)
    p["h2o"] = prof["h2o"] * f
    new = atmosphere.layers(p, T.GROUND_KM)["col"]["H2O"].sum()
    p["h2o"] = p["h2o"] * base / new
    return p

def run(args):
    prof, w0, w1, step = args
    L = atmosphere.layers(prof, surface_km=T.GROUND_KM, sensor_km=T.SENSOR_KM)
    groups = {"all": np.ones(len(L["p"]), bool), "below": L["below_sensor"]}
    return gas.column_tau(L, "H2O", w0, w1, groups, wstep=step)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("base"); ap.add_argument("out")
    ap.add_argument("--low", type=float, default=1.6); ap.add_argument("--high", type=float, default=0.6)
    ap.add_argument("--zsplit", type=float, default=3.0)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    for g in ("CO2", "CH4", "O2"):
        shutil.copy(Path(a.base) / f"tau_{g}.npz", out)
    z = np.load(Path(a.base) / "tau_H2O.npz"); wn = z["wn"]; step = wn[1] - wn[0]
    i0, i1 = np.searchsorted(wn, [1e7 / 1270, 1e7 / 850])
    edges = np.linspace(i0, i1, 9).astype(int)
    prof = variant(a.low, a.high, a.zsplit)
    jobs = [(prof, wn[e0], wn[e1 - 1] + step / 2, step) for e0, e1 in zip(edges[:-1], edges[1:])]
    allv, below = z["all"].copy(), z["below"].copy()
    with Pool(8) as p:
        for (e0, e1), (w, o) in zip(zip(edges[:-1], edges[1:]), p.map(run, jobs)):
            n = min(len(w), e1 - e0)
            allv[e0:e0 + n], below[e0:e0 + n] = o["all"][:n], o["below"][:n]
    np.savez(out / "tau_H2O.npz", wn=wn, all=allv, below=below)
    print("written", out)
