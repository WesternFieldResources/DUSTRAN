"""Leftover gas ripple in reflectance cubes of the same pixels: lower means the gas was removed better.

    python experiments/gas_ripple.py LABEL=CUBE_HDR[@ROW0] [LABEL=CUBE_HDR[@ROW0] ...] --rows R0 R1

Each cube is read over the same radiance rows (ROW0 is the cube's first radiance row, 0 for a full
line). Over bare bright pixels (NDVI < 0.2, 30-95th percentile at 1650 nm), each pixel's ln(reflectance)
in each window gets a quadratic removed; the rms of what is left, as a percentage, is the ripple. Smooth
mineral bands survive the quadratic only partly, so read the numbers relative to each other on the
same pixels, never as absolute residue. The same is done on the median spectrum.
"""
import argparse, json
import numpy as np
from dustran.batch import read_header

WINDOWS = {"O2 760": (745, 785), "H2O 820": (790, 865), "H2O 940": (870, 1000), "H2O 1140": (1080, 1200),
           "CO2 1600": (1540, 1640), "CO2 2010/2060": (1985, 2085), "CH4/H2O 2300-2420": (2290, 2420)}


def load(spec, r0, r1):
    path, _, row0 = spec.partition("@")
    row0 = int(row0 or 0)
    _, h, arr = read_header(path)
    nl, nb, ns = int(h["lines"]), int(h["bands"]), int(h["samples"])
    a = np.memmap(path[:-4], "<f4", "r", shape=(nl, nb, ns))
    x = np.asarray(a[r0 - row0:r1 - row0], np.float64).transpose(0, 2, 1).reshape(-1, nb)
    return arr("wavelength"), x


def ripple(c, X, lo, hi):
    m = (c >= lo) & (c <= hi)
    x = (c[m] - c[m].mean()) / 50
    A = np.c_[np.ones_like(x), x, x ** 2]
    H = np.eye(m.sum()) - A @ np.linalg.pinv(A)
    y = np.log(np.clip(X[:, m], 1e-4, None))
    return 100 * np.sqrt(np.mean((y @ H.T) ** 2, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cubes", nargs="+"); ap.add_argument("--rows", nargs=2, type=int, required=True)
    ap.add_argument("--json")
    a = ap.parse_args()
    data = {}
    for s in a.cubes:
        lab, _, spec = s.partition("=")
        data[lab] = load(spec, *a.rows)
    c, X0 = next(iter(data.values()))
    i = lambda w: int(np.argmin(np.abs(c - w)))
    ok = np.all([(X[:, i(1650)] > 0) & (X[:, i(860)] > 0) for _, X in data.values()], 0)
    nd = (X0[:, i(860)] - X0[:, i(660)]) / (X0[:, i(860)] + X0[:, i(660)] + 1e-9)
    b = X0[:, i(1650)]
    lo, hi = np.percentile(b[ok], [30, 95])
    m = ok & (nd < 0.2) & (b >= lo) & (b <= hi)
    print(f"{m.sum()} bare bright pixels")
    out = {}
    print(f"{'window':20s}" + "".join(f"{k:>16s}" for k in data) + "   (median per-pixel ripple %; median-spectrum ripple %)")
    for wn, (lo_, hi_) in WINDOWS.items():
        row = []
        for k, (cc, X) in data.items():
            p = float(np.median(ripple(cc, X[m], lo_, hi_)))
            q = float(ripple(cc, np.median(X[m], 0)[None], lo_, hi_)[0])
            out.setdefault(wn, {})[k] = dict(pixel=p, median_spectrum=q)
            row.append(f"{p:7.2f};{q:6.2f}  ")
        print(f"{wn:20s}" + "".join(f"{r:>16s}" for r in row))
    if a.json:
        open(a.json, "w").write(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
