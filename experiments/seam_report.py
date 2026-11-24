"""Seam check (free relative calibration) for the validation scenes.

    python experiments/seam_report.py

Runs dustran.instrument.seam_check on the playa median radiance of each scene. Part of every
DUSTRAN workflow: report the seams, flag disagreement over 1 %, never silently correct.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import playa_gas_score as P
from dustran import instrument

t = np.loadtxt(P.TSIS)
sun = lambda c, f: instrument.convolve(1e7 / t[::-1, 0], t[::-1, 1], c, f)


def report(name, c, f, L):
    res = instrument.seam_check(c, f, L, sun(c, f))
    print(f"{name}: {len(res)} seam(s)")
    for s in res:
        flag = "  <-- FLAG" if s["flagged"] else (f"  (low signal: {100 * s['signal_level']:.0f}% of typical, opaque band)" if s["low_signal"] else "")
        print(f"  {s['ends_nm']:.0f} | {s['starts_nm']:.0f} nm: median ratio {s['ratio']:.3f}"
              f" from {len(s['pairs'])} pairs{flag}")
    return res


if __name__ == "__main__":
    import classic_playa as C
    S = C.scene(); report("AVIRIS Classic f190813t01p00r07 (playa)", S["centers"], S["fwhm"], S["L"])
    S = P.load_scene(); report("AVIRIS-NG ang20200712t214704 (playa)", S["centers"], S["fwhm"], S["L"])
