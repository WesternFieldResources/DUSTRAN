"""Step 1 (go/no-go), part A: line-by-line gas optical depth for the Railroad Valley playa scene.

    python experiments/playa_gas_tau.py OUTDIR [--gases H2O CO2 CH4 O2] [--wl 740 2470]

Scene: AVIRIS-NG ang20200712t214704, ground 1.436 km, sensor 5.313 km (from the NASA L1B navigation).
Atmosphere: AFGL US standard (Anderson et al. 1986).
Writes OUTDIR/tau_<gas>.npz with wn and two columns: 'all' (ground to top of atmosphere, the
sun's path) and 'below' (ground to sensor, the upward path), water at the standard profile
(rescale per fit). One process per gas.
"""
import argparse, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dustran import atmosphere, gas

from dustran.paths import DATA
AFGL = str(DATA / "afgl/afglus.dat")   # AFGL US standard (Anderson et al. 1986)
GROUND_KM, SENSOR_KM = 1.4357627632017484, 5.313446641146695


def run(args):
    g, outdir, wn0, wn1, wstep, ground, sensor = args
    L = atmosphere.layers(atmosphere.load_afgl(AFGL), surface_km=ground, sensor_km=sensor)
    groups = {"all": np.ones(len(L["p"]), bool), "below": L["below_sensor"]}
    t = time.time()
    logf = open(Path(outdir) / f"log_{g}.txt", "w")
    log = lambda s: (logf.write(f"{time.time()-t:7.0f}s {s}\n"), logf.flush())
    wn, out = gas.column_tau(L, g, wn0, wn1, groups, wstep=wstep, log=log)
    np.savez(Path(outdir) / f"tau_{g}.npz", wn=wn, **out)
    log("saved")
    return g, time.time() - t


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--gases", nargs="+", default=["H2O", "CO2", "CH4", "O2"])
    ap.add_argument("--wl", nargs=2, type=float, default=[740, 2470])
    ap.add_argument("--wstep", type=float, default=0.002)
    ap.add_argument("--ground", type=float, default=GROUND_KM, help="ground elevation [km]")
    ap.add_argument("--sensor", type=float, default=SENSOR_KM, help="sensor altitude [km]")
    a = ap.parse_args()
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    wn0, wn1 = 1e7 / a.wl[1], 1e7 / a.wl[0]
    with Pool(len(a.gases)) as p:
        for g, dt in p.imap_unordered(run, [(g, a.outdir, wn0, wn1, a.wstep, a.ground, a.sensor) for g in a.gases]):
            print(g, f"{dt:.0f}s", flush=True)
