"""Pull the median radiance and geometry for a ground point out of an AVIRIS Classic ortho scene,
and run the seam check.

    python experiments/classic_site.py <rdn_dir> <utm_e> <utm_n> [half_m=185]

Reads the ortho radiance (int16, BIP, big-endian, divided by the per-band .gain) and the ortho
obs cube (path, to-sensor az/zen, solar az/zen, ..., UTC) at every ortho pixel within +-half_m
of the point. Writes site.npz next to the data: L (uW cm-2 nm-1 sr-1), centers, fwhm, obs median,
sensor altitude (from the ephemeris), pixel count, and the seam check.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dustran import instrument


def hdr(path):
    t = open(path).read()
    g = lambda k: t.split(k)[1].split("\n")[0].replace("=", "").strip()
    mi = t.split("map info")[1].split("{")[1].split(",")
    return int(g("samples")), int(g("lines")), int(g("bands")), float(mi[3]), float(mi[4]), float(mi[5])


def tsis_band(centers, fwhm):
    T = Path(__file__).resolve().parent
    import playa_gas_score as P
    t = np.loadtxt(P.TSIS)
    return instrument.convolve(1e7 / t[::-1, 0], t[::-1, 1], centers, fwhm)


if __name__ == "__main__":
    d = Path(sys.argv[1]); e, n = float(sys.argv[2]), float(sys.argv[3])
    half = float(sys.argv[4]) if len(sys.argv) > 4 else 185.0
    name = d.name
    img = next(d.glob("*_sc01_ort_img"))
    ns, nl, nb, e0, n0, px = hdr(str(img) + ".hdr")
    x0, y0 = (e - e0) / px, (n0 - n) / px
    r = int(round(half / px))
    ys, xs = np.mgrid[int(y0) - r:int(y0) + r + 1, int(x0) - r:int(x0) + r + 1]
    ok = (ys >= 0) & (ys < nl) & (xs >= 0) & (xs < ns)
    ys, xs = ys[ok], xs[ok]
    cube = np.memmap(img, ">i2", "r", shape=(nl, ns, nb))
    dn = np.stack([np.asarray(cube[y, x]) for y, x in zip(ys, xs)]).astype(float)
    valid = (dn[:, [10, 50, 150]] > 0).all(1)          # fill value is -50 outside the swath
    gain = np.loadtxt(next(d.glob("*.gain")))[:, 0]
    L = np.median(dn[valid], 0) / gain
    obs_ns, obs_nl, obs_nb, *_ = hdr(str(d / f"{name}_obs_ort") + ".hdr")
    obs = np.memmap(d / f"{name}_obs_ort", "<f4", "r", shape=(obs_nl, obs_nb, obs_ns))
    o = np.median(np.stack([np.asarray(obs[y, :, x]) for y, x in zip(ys[valid], xs[valid])]), 0)
    spc = np.loadtxt(next(d.glob("*.spc")))
    eph = np.fromfile(d / f"{name}_eph", "<f8").reshape(-1, 6)
    alt = float(np.median(eph[:, 2])) / 1000.0
    sc = instrument.seam_check(spc[:, 0], spc[:, 1], L, tsis_band(spc[:, 0], spc[:, 1]))
    np.savez(d.parent / "site.npz", L=L, centers=spc[:, 0], fwhm=spc[:, 1], obs=o, sensor_km=alt,
             npx=int(valid.sum()), seam_ratio=[s["ratio"] for s in sc], seam_low=[s["low_signal"] for s in sc])
    print(f"{name}: {valid.sum()} px, UTC {o[9]:.3f}, sza {o[4]:.1f}, vza {o[2]:.1f}, sensor {alt:.2f} km")
    for s in sc:
        print(f"  seam {s['ends_nm']:.0f}|{s['starts_nm']:.0f} nm: ratio {s['ratio']:.4f}"
              f"{'  (opaque band, ignore)' if s['low_signal'] else ''}{'  FLAG' if s['flagged'] else ''}")
