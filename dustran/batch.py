"""Whole-scene atmospheric correction: every pixel of a flight line at once, in numpy.

    python -m dustran.batch RADIANCE OUT_PREFIX --taudirs DIR[,DIR] --ground KM --sensor KM
                            [--obs OBS_ORT_HDR] [--shift A B | --shift fit] [--year 2017]
                            [--o3-du 282] [--aeronet DIR --inv-time T --sun-time T --date DD:MM:YYYY]
                            [--aod 0.05] [--pressure-ratio 1.0] [--lines R0 R1] [--bbl-from HDR]

RADIANCE is either an AVIRIS-NG ortho radiance header (*_rdn_*_img.hdr, float32 BIL, with --obs
pointing at its *_obs_ort.hdr) or an unpacked AVIRIS Classic radiance folder (fYYMMDDtNNpNNrNN_rdn_X:
int16 BIP big-endian scaled by the .gain file, band table from .spc, *_obs_ort beside it).
`classic_geometry(folder)` returns the ground elevation and sensor altitude for the gas run.

Same physics as the validated playa chain (experiments/playa_rayleigh.py, classic_playa.py),
reorganised so that the per-pixel work is array arithmetic:

Once per scene
  * Line-by-line gas optical depths (from experiments/playa_gas_tau.py at the scene's ground
    elevation and sensor altitude), O2 collision-induced absorption, ozone (Bremen cross sections
    scaled to a measured column), TSIS-1 solar spectrum. CO2 scaled to the flight year.
  * Band wavelength shift: given (A + B (nm - 1500)), or fitted from the scene itself ("fit"): one
    shift per spectrometer (Classic has four, split at its seams), chosen to make the gas-corrected
    median spectrum smoothest across the O2 760 nm and CO2 1600 nm bands. Spectrometers
    with no usable gas band take the nearest fitted one's shift.
  * A small gas table built on the fly (no stored lookup tables): ln band<E_sun T_gas> for
    ~40 water scales x a few view zenith angles.
  * Scattering terms R_path, T_dn T_up, S from dustran/sos.py on a coarse wavelength grid, for a
    few view geometries (view zenith x relative azimuth); Rayleigh + aerosol (AERONET-measured
    sizes and index, or a lognormal fitted to an Angstrom exponent).
  * Seam check (spectrometer overlaps): reported and flagged at 1 %, never corrected.

Per pixel (vectorised, in row blocks)
  * Water: the scale that makes the gas-corrected spectrum flattest across the 940 nm band,
    relative to a straight line through clean bands on either side (as experiments/playa_water.py),
    picked from the table and refined with a parabola. 1140 nm is fitted the same way as a
    cross-check; both go into the water map, and their disagreement is the scene's water flag.
  * r = gas-corrected apparent reflectance at that water; x = r - R_path;
    rho = x / (T_dn T_up + S x).

Bad-band list (bbl) in the output, first available of:
  1. the data provider's: the radiance header's bbl, or one copied with --bbl-from;
  2. this scene's: bands whose two-way gas transmittance at the scene's median water is < 1 %
     (no light reaches the sensor there: no data, never data);
  3. the conservative default CORE_OPAQUE_NM (empty: see its note).

Output: ENVI BIL float32 reflectance cube (band table, map info, bbl), a 3-band water map
(940 cm, 1140 cm, 940 fit rms %), and a JSON summary. Ground elevation, solar zenith for the gas
paths and the aerosol are single values per scene.
"""
import argparse, datetime, json, re, sys, time
from pathlib import Path
import numpy as np

from dustran import atmosphere, aerosol, cia, instrument, rayleigh, sos, xsec
from dustran.paths import DATA

AFGL = DATA / "afgl/afglus.dat"
TSIS = DATA / "tsis/tsis_f0_0p1.txt"
CIA = {DATA / "hitran_xsec/O2-O2_2024.cia": ("O2", "O2"), DATA / "hitran_xsec/O2-N2_2024.cia": ("O2", "N2")}
O3 = DATA / "bremen_o3/serdyuchenkogorshelev5digits.dat"
H2O_CM = 1 / 3.3428e22                    # molecules cm^-2 -> g cm^-2 (precipitable cm)
WATER = {"940": ((880, 1000), ((850, 880), (1000, 1040))),
         "1140": ((1080, 1180), ((1040, 1075), (1230, 1260)))}
SCAT_GRID = np.r_[np.arange(400, 1000, 25.0), np.arange(1000, 2501, 150.0)]
# CO2 2010/2060 is not used: on AVIRIS Classic (10 nm bands) its roughness minimum moved the SWIR
# spectrometer by +2.7 nm at Cerbat (2017-12-06), which erased the real 2200 nm Al-OH band and
# made a false 2165 nm one (checked against NASA reflectance at the same pixels).
SHIFT_WINDOWS = {"O2 760": (745, 785), "CO2 1600": (1550, 1630)}
OPAQUE_T = 0.01                           # two-way band transmittance below which a band is no data
# Conservative default: band ranges with no data on any day. Derived with experiments/core_opaque.py
# from DUSTRAN's line-by-line water (ground 1.5 km, sun and view at zenith, 5 and 10 nm bands):
# even the 1380 / 1880 nm cores pass 2-4 % two-way at standard water (light leaks between lines),
# so no band is "never data" at a high desert site. Empty on purpose; scene and provider lists override.
CORE_OPAQUE_NM = ()


# ---------------------------------------------------------------- ENVI helpers
def read_header(path):
    text = Path(path).read_text()
    hdr = {m.group(1).strip().lower(): m.group(2).strip()
           for m in re.finditer(r"^\s*([^=\n]+?)\s*=\s*(\{[^}]*\}|[^\n]*)", text, re.M)}
    arr = lambda k: np.array([float(v) for v in hdr[k].strip("{}").split(",")]) if k in hdr else None
    return text, hdr, arr


def write_header(path, samples, lines, description, map_info, centers=None, fwhm=None, bbl=None,
                 band_names=None, coord=None):
    """A fresh ENVI header for a float32 BIL little-endian file."""
    nb = len(centers) if centers is not None else len(band_names)
    fmt = lambda a: "{ " + " , ".join(f"{float(x):.6g}" for x in a) + " }"
    t = ["ENVI", "description = {\n  " + description + "}", f"samples = {samples}", f"lines = {lines}",
         f"bands = {nb}", "header offset = 0", "file type = ENVI Standard", "data type = 4",
         "interleave = bil", "byte order = 0"]
    if map_info:
        t.append("map info = {" + " ".join(map_info.strip("{}").split()) + "}")
    if coord:
        t.append(f"coordinate system string = {coord}")
    if centers is not None:
        t += ["wavelength units = Nanometers", "wavelength = " + fmt(centers), "fwhm = " + fmt(fwhm)]
    if bbl is not None:
        t.append("bbl = " + fmt(bbl))
    if band_names:
        t.append("band names = {" + ", ".join(band_names) + "}")
    t.append("data ignore value = -9999")
    Path(path).write_text("\n".join(t) + "\n")


# ---------------------------------------------------------------- readers
def earth_sun_au(doy):
    return 1 - 0.01672 * np.cos(np.radians(0.9856 * (doy - 4)))


class NGReader:
    """AVIRIS-NG: float32 BIL radiance [uW cm-2 nm-1 sr-1], float64 BIP obs_ort (11 bands)."""

    def __init__(self, hdr_path, obs_hdr):
        _, h, arr = read_header(hdr_path)
        self.centers, self.fwhm = arr("wavelength"), arr("fwhm")
        self.bbl = arr("bbl")
        self.nl, self.nb, self.ns = int(h["lines"]), int(h["bands"]), int(h["samples"])
        self.map_info, self.coord = h.get("map info"), h.get("coordinate system string")
        self.rdn = np.memmap(str(hdr_path)[:-4], "<f4", "r", shape=(self.nl, self.nb, self.ns))
        _, ho, _ = read_header(obs_hdr)
        nbo = int(ho["bands"])
        self.obs = np.memmap(str(obs_hdr)[:-4], "<f8", "r", shape=(self.nl, self.ns, nbo))
        self.name = Path(hdr_path).name.split("_")[0]
        o = self.sample_obs()
        self.esd = float(np.median(o[:, 10]))

    def sample_obs(self, step=7):
        o = np.asarray(self.obs[::step, ::step, :]).reshape(-1, self.obs.shape[2])
        return o[o[:, 0] > 0]

    def block(self, r0, r1):
        L = np.asarray(self.rdn[r0:r1], np.float64).transpose(0, 2, 1).reshape(-1, self.nb)
        g = np.asarray(self.obs[r0:r1], np.float64).reshape(-1, self.obs.shape[2])
        ok = (L[:, 0] > -9000) & (g[:, 0] > 0) & (np.nanmax(L, 1) > 0)
        return L, g[:, :5], ok


class ClassicReader:
    """AVIRIS Classic: int16 BIP big-endian radiance / .gain = uW cm-2 nm-1 sr-1 (fill -50 outside
    the swath), band table from .spc, float32 BIL obs_ort (10 bands, no Earth-Sun distance)."""

    def __init__(self, folder):
        d = Path(folder)
        self.dir = d
        img = next(d.glob("*_sc01_ort_img"))
        _, h, _ = read_header(str(img) + ".hdr")
        self.nl, self.nb, self.ns = int(h["lines"]), int(h["bands"]), int(h["samples"])
        assert h.get("interleave", "").lower() == "bip" and h.get("byte order") == "1"
        self.map_info, self.coord = h.get("map info"), h.get("coordinate system string")
        self.rdn = np.memmap(img, ">i2", "r", shape=(self.nl, self.ns, self.nb))
        self.gain = np.loadtxt(next(d.glob("*.gain")))[:, 0]
        spc = np.loadtxt(next(d.glob("*.spc")))
        self.centers, self.fwhm = spc[:, 0], spc[:, 1]
        self.bbl = None
        obs = next(d.glob("*_obs_ort"))
        _, ho, _ = read_header(str(obs) + ".hdr")
        assert int(ho["lines"]) == self.nl and int(ho["samples"]) == self.ns
        self.obs = np.memmap(obs, "<f4", "r", shape=(self.nl, int(ho["bands"]), self.ns))
        self.name = d.name.split("rdn")[0].rstrip("_")
        yy, mm, dd = int(self.name[1:3]), int(self.name[3:5]), int(self.name[5:7])
        self.date = datetime.date(2000 + yy if yy < 80 else 1900 + yy, mm, dd)
        self.esd = float(earth_sun_au(self.date.timetuple().tm_yday))

    def sample_obs(self, step=7):
        o = np.asarray(self.obs[::step, :, ::step]).transpose(0, 2, 1).reshape(-1, self.obs.shape[1])
        return o[o[:, 0] > 0]

    def block(self, r0, r1):
        dn = np.asarray(self.rdn[r0:r1]).reshape(-1, self.nb).astype(np.float64)
        g = np.asarray(self.obs[r0:r1], np.float64).transpose(0, 2, 1).reshape(-1, self.obs.shape[1])
        ok = (dn[:, [10, 50, 150]] > 0).all(1) & (g[:, 0] > 0)
        return dn / self.gain, g[:, :5], ok


def classic_geometry(folder):
    """Median ground elevation (IGM) and sensor altitude (ephemeris) [km] of a Classic line."""
    d = Path(folder)
    igm = next(d.glob("*_ort_igm"))
    _, h, _ = read_header(str(igm) + ".hdr")
    nl, ns = int(h["lines"]), int(h["samples"])
    il = h.get("interleave", "bil").lower()
    a = np.memmap(igm, "<f8", "r", shape=(nl, 3, ns) if il == "bil" else (nl, ns, 3))
    z = np.asarray(a[::5, 2, ::5] if il == "bil" else a[::5, ::5, 2]).ravel()
    z = z[(z > -500) & (z < 9000)]
    eph = np.fromfile(next(p for p in d.glob("*_eph") if "lonlat" not in p.name), "<f8").reshape(-1, 6)
    return float(np.median(z)) / 1000.0, float(np.median(eph[:, 2])) / 1000.0


def open_radiance(path, obs=None):
    p = Path(path)
    if p.is_dir():
        return ClassicReader(p)
    return NGReader(p, obs)


# ---------------------------------------------------------------- once per scene
class Scene:
    """Everything that does not depend on the pixel."""

    def __init__(self, centers, fwhm, taudirs, ground_km, sensor_km, sza, esd, vza_nodes, dphi_nodes,
                 shift=(0.389, 0.063e-3), L_med=None, o3_du=None, aer=None, aod=None, aod_ref=550.0,
                 pressure_ratio=1.0, ws=np.linspace(0.1, 2.5, 40), year=2020, log=print):
        self.c, self.f = np.asarray(centers, float), np.asarray(fwhm, float)
        self.ws, self.vza = np.asarray(ws, float), np.asarray(vza_nodes, float)
        self.dphi = np.asarray(dphi_nodes, float)
        self.mu_s = np.cos(np.radians(sza))
        self.esd = esd
        prof = atmosphere.load_afgl(str(AFGL))
        self.L = atmosphere.layers(prof, surface_km=ground_km, sensor_km=sensor_km)
        self.L["air"] = self.L["air"] * pressure_ratio
        self.cm1 = atmosphere.layers(prof, ground_km)["col"]["H2O"].sum() * H2O_CM
        groups = {"all": np.ones(len(self.L["p"]), bool), "below": self.L["below_sensor"]}
        t0 = time.time()
        wn, E, dry_all, dry_below, wet_all, wet_below = self._fine(taudirs, groups,
                                                                   atmosphere.co2_ppm(year) / 330.2, o3_du)
        log(f"fine grid: {len(wn)} points, {time.time() - t0:.0f}s")
        mu_v = np.cos(np.radians(self.vza[len(self.vza) // 2]))
        if isinstance(shift, str) and shift == "fit":
            t0 = time.time()
            self.shift, self.shift_fit = self._fit_shift(wn, E, dry_all, dry_below, wet_all, wet_below,
                                                         mu_v, L_med)
            log(f"shift fit ({time.time() - t0:.0f}s): " + "; ".join(
                f"{k}: {v['shift']:+.2f} nm" + (f" (from {v['from']}; FLAG {v['flag']})" if 'from' in v else "")
                for k, v in self.shift_fit.items()))
        else:
            self.shift = shift[0] + shift[1] * (self.c - 1500)
            self.shift_fit = None
        # gas table: ln band<E T> over (view node, water scale, band)
        t0 = time.time()
        lnET = np.empty((len(self.vza), len(self.ws), len(self.c)))
        lo, hi = 1e7 / wn.max(), 1e7 / wn.min()
        inside = (self.c + self.shift - 2.5 * self.f >= lo) & (self.c + self.shift + 2.5 * self.f <= hi)
        t = np.loadtxt(TSIS)
        sun_only = instrument.convolve(1e7 / t[::-1, 0], t[::-1, 1], self.c, self.f, self.shift)
        self.E_band = sun_only
        for i, v in enumerate(self.vza):
            mv = np.cos(np.radians(v))
            dry = np.exp(-(dry_all / self.mu_s + dry_below / mv)) * E
            wet = wet_all / self.mu_s + wet_below / mv
            for j, w in enumerate(self.ws):
                m = instrument.convolve(wn, dry * np.exp(-w * wet), self.c, self.f, self.shift)
                m = np.where(inside & np.isfinite(m), m, sun_only)
                lnET[i, j] = np.log(np.maximum(m, 1e-30))
        self.lnET = lnET
        log(f"gas table {lnET.shape}: {time.time() - t0:.0f}s")
        t0 = time.time()
        self.scat = self._scattering(aer, aod, aod_ref)
        log(f"scattering {len(SCAT_GRID)} wavelengths x {len(self.vza)} x {len(self.dphi)} geometries: "
            f"{time.time() - t0:.0f}s")

    def transmittance(self, w):
        """Two-way gas transmittance per band at water scale w (middle view node)."""
        j, f = _bracket(self.ws, np.array([w]))
        i = len(self.vza) // 2
        ln = self.lnET[i, j[0]] * (1 - f[0]) + self.lnET[i, j[0] + 1] * f[0]
        return np.exp(ln) / self.E_band

    def _fine(self, taudirs, groups, co2_scale, o3_du):
        segs = []
        for d in taudirs:
            z = {g: np.load(f) for g in ("H2O", "CO2", "CH4", "O2") if (f := Path(d) / f"tau_{g}.npz").exists()}
            wn = z["H2O"]["wn"]
            sc = lambda g: co2_scale if g == "CO2" else 1.0
            segs.append(dict(wn=wn,
                             da=sum(sc(g) * z[g]["all"].astype(np.float64) for g in z if g != "H2O"),
                             db=sum(sc(g) * z[g]["below"].astype(np.float64) for g in z if g != "H2O"),
                             wa=z["H2O"]["all"].astype(np.float64), wb=z["H2O"]["below"].astype(np.float64)))
        segs.sort(key=lambda s: s["wn"][0])
        keep = [segs[0]]
        for s in segs[1:]:
            m = s["wn"] > keep[-1]["wn"][-1]
            keep.append({k: v[m] for k, v in s.items()})
        s = {k: np.concatenate([x[k] for x in keep]) for k in keep[0]}
        wn = s["wn"]
        c = cia.optical_depth(self.L, CIA, wn, groups)
        s["da"] += c["all"]; s["db"] += c["below"]
        scale = (o3_du / xsec.dobson(self.L)) if o3_du else 1.0
        o = xsec.optical_depth(self.L, xsec.read_bremen(str(O3)), "O3", wn, groups, column_scale=scale)
        s["da"] += o["all"]; s["db"] += o["below"]
        t = np.loadtxt(TSIS)
        E = np.interp(1e7 / wn, t[:, 0], t[:, 1])
        return wn, E, s["da"], s["db"], s["wa"], s["wb"]

    def _fit_shift(self, wn, E, dry_all, dry_below, wet_all, wet_below, mu_v, L_med):
        """One shift per spectrometer from the scene's median radiance: minimise the roughness
        (rms 2nd difference / median) of the gas-corrected spectrum across O2 / CO2 bands."""
        assert L_med is not None, "shift fit needs the scene's median radiance"
        y = np.exp(-(dry_all / self.mu_s + dry_below / mu_v)) * E * np.exp(-(wet_all / self.mu_s + wet_below / mu_v))
        bounds = [0] + instrument.seams(self.c) + [len(self.c)]
        shifts = np.full(len(self.c), np.nan)
        info, cand = {}, np.round(np.arange(-4.0, 4.001, 0.1), 2)
        for a, b in zip(bounds[:-1], bounds[1:]):
            name = f"{self.c[a]:.0f}-{self.c[b - 1]:.0f} nm"
            idx = np.arange(a, b)
            wins = [idx[(self.c[idx] >= lo) & (self.c[idx] <= hi)] for lo, hi in SHIFT_WINDOWS.values()]
            wins = [w for w in wins if len(w) >= 4]
            if not wins:
                info[name] = None
                continue
            sel = np.concatenate(wins)
            cost = []
            for s in cand:
                m = instrument.convolve(wn, y, self.c[sel], self.f[sel], s)
                r = L_med[sel] / m
                parts, k = [], 0
                for w in wins:
                    q = r[k:k + len(w)]; k += len(w)
                    parts.append(np.sqrt(np.mean(np.diff(q, 2) ** 2)) / np.median(q))
                cost.append(np.sqrt(np.mean(np.square(parts))))
            cost = np.array(cost); j0 = int(np.argmin(cost))
            if j0 in (0, len(cand) - 1):              # minimum on the grid edge: not a fit
                info[name] = None
                continue
            j = j0
            den = cost[j - 1] - 2 * cost[j] + cost[j + 1]
            d = np.clip(0.5 * (cost[j - 1] - cost[j + 1]) / den, -1, 1) if den > 0 else 0.0
            s = cand[j] + d * 0.1
            shifts[a:b] = s
            info[name] = dict(shift=float(s), roughness_pct=float(100 * cost.min()),
                              windows=[k for k, (lo, hi) in SHIFT_WINDOWS.items()
                                       if any((self.c[idx] >= lo) & (self.c[idx] <= hi))])
        fitted = [(a, b) for a, b in zip(bounds[:-1], bounds[1:]) if np.isfinite(shifts[a])]
        for (a, b), name in zip(zip(bounds[:-1], bounds[1:]), list(info)):
            if info[name] is None:                    # nearest fitted spectrometer
                na, nb_ = min(fitted, key=lambda ab: abs(ab[0] - a))
                shifts[a:b] = shifts[na]
                info[name] = dict(shift=float(shifts[na]), **{"from": f"{self.c[na]:.0f}-{self.c[nb_ - 1]:.0f} nm"},
                                  flag="no usable gas band or no minimum inside +-4 nm")
        return shifts, info

    def _scattering(self, aer, aod, ref):
        L = self.L
        ks = int(np.argmax(L["below_sensor"]))
        out = np.empty((len(self.vza), len(self.dphi), 4, len(SCAT_GRID)))
        for k, wv in enumerate(SCAT_GRID):
            tR = rayleigh.layer_tau(L, wv)
            if aod:
                o = aer.optics(wv)
                bR = np.zeros(len(o["beta"])); bR[:3] = rayleigh.phase_moments(wv)
                tA = aerosol.layer_tau(L, aod, o["ext"] / aer.optics(ref)["ext"])
                tau = tR + tA; sca = tR + o["ssa"] * tA
                beta = (tR[:, None] * bR + (o["ssa"] * tA)[:, None] * o["beta"]) / sca[:, None]
                ssa, n = sca / tau, 16
            else:
                tau, ssa, beta, n = tR, np.ones_like(tR), np.tile(rayleigh.phase_moments(wv), (len(tR), 1)), 8
            for i, v in enumerate(self.vza):
                for j, p in enumerate(self.dphi):
                    r = sos.surface_terms(tau, ssa, beta, self.mu_s, np.cos(np.radians(v)), p, ks, n_streams=n)
                    out[i, j, :, k] = r["R_path"], r["T_dn"] * r["T_up"], r["S"], r["T_dn"]
        # to bands: (vza, dphi, term, band)
        return np.stack([[[np.interp(self.c, SCAT_GRID, out[i, j, t]) for t in range(4)]
                          for j in range(len(self.dphi))] for i in range(len(self.vza))])


# ---------------------------------------------------------------- per pixel
def _bracket(nodes, x):
    """Index i and fraction f with x ~ nodes[i] + f (nodes[i+1]-nodes[i]); clipped to the ends."""
    x = np.clip(x, nodes[0], nodes[-1])
    i = np.clip(np.searchsorted(nodes, x) - 1, 0, len(nodes) - 2)
    return i, (x - nodes[i]) / (nodes[i + 1] - nodes[i])


def _line_operator(wl, inw, edge):
    """Matrix H so that H @ y[edge] is the straight line through y[edge] evaluated at wl[inw]."""
    A = np.c_[np.ones(edge.sum()), wl[edge]]
    return np.c_[np.ones(inw.sum()), wl[inw]] @ np.linalg.pinv(A)


def fit_water(sc, r0, iv, fv, window, good):
    """Per-pixel water scale from one band. r0: (P, B) radiance-based numerator k*L*10 (so that
    r = r0 / band<E T>). Returns scale w (P,) and the rms % at the minimum."""
    (lo, hi), ((a0, a1), (b0, b1)) = window
    wl = sc.c
    inw = good & (wl >= lo) & (wl <= hi)
    edge = good & (((wl >= a0) & (wl <= a1)) | ((wl >= b0) & (wl <= b1)))
    sel = inw | edge
    H = _line_operator(wl[sel], inw[sel], edge[sel])
    ii, ie = inw[sel], edge[sel]
    ln_r0 = np.log(np.maximum(r0[:, sel], 1e-12))
    lnET = sc.lnET[:, :, sel]
    cost = np.empty((len(r0), len(sc.ws)))
    for j in range(len(sc.ws)):
        lnm = lnET[iv, j] * (1 - fv[:, None]) + lnET[iv + 1, j] * fv[:, None]
        r = np.exp(ln_r0 - lnm)
        line = r[:, ie] @ H.T
        cost[:, j] = np.sqrt(np.mean((r[:, ii] / line - 1) ** 2, axis=1))
    j = np.clip(np.argmin(cost, 1), 1, len(sc.ws) - 2)
    P = np.arange(len(r0))
    c0, c1, c2 = cost[P, j - 1], cost[P, j], cost[P, j + 1]
    den = c0 - 2 * c1 + c2
    d = np.where(den > 0, 0.5 * (c0 - c2) / np.where(den > 0, den, 1), 0.0)
    d = np.clip(d, -1, 1)
    step = sc.ws[1] - sc.ws[0]
    return sc.ws[j] + d * step, 100 * np.min(cost, 1)


def correct_block(sc, Lrad, vza, dphi, sza_px, good):
    """Lrad (P, B) radiance [uW cm-2 sr-1 nm-1]; vza, dphi, sza_px (P,). Returns rho (P,B),
    water scale at 940 and 1140 (P,), 940 fit rms %."""
    k = np.pi * sc.esd ** 2 / np.cos(np.radians(sza_px))
    r0 = k[:, None] * Lrad * 10.0                      # uW cm-2 -> mW m-2 (TSIS units)
    iv, fv = _bracket(sc.vza, vza)
    w940, rms940 = fit_water(sc, r0, iv, fv, WATER["940"], good)
    w1140, _ = fit_water(sc, r0, iv, fv, WATER["1140"], good)
    jw, fw = _bracket(sc.ws, w940)
    T = sc.lnET
    lnm = ((T[iv, jw] * (1 - fw[:, None]) + T[iv, jw + 1] * fw[:, None]) * (1 - fv[:, None])
           + (T[iv + 1, jw] * (1 - fw[:, None]) + T[iv + 1, jw + 1] * fw[:, None]) * fv[:, None])
    r = r0 * np.exp(-lnm)
    ip, fp = _bracket(sc.dphi, dphi)
    S4 = sc.scat
    wv, wp = fv[:, None, None], fp[:, None, None]
    terms = ((S4[iv, ip] * (1 - wp) + S4[iv, ip + 1] * wp) * (1 - wv)
             + (S4[iv + 1, ip] * (1 - wp) + S4[iv + 1, ip + 1] * wp) * wv)     # (P, 4, B)
    x = r - terms[:, 0]
    rho = x / (terms[:, 1] + terms[:, 2] * x)
    return rho, w940, w1140, rms940


def relative_azimuth(vaz, saz):
    """|view azimuth - (sun azimuth + 180)| folded into [0, 180] degrees."""
    return np.abs(((vaz - (saz + 180.0)) + 180.0) % 360.0 - 180.0)


def core_opaque(centers):
    c = np.asarray(centers, float)
    out = np.zeros(len(c), bool)
    for lo, hi in CORE_OPAQUE_NM:
        out |= (c >= lo) & (c <= hi)
    return out


def median_radiance(rd, n=4000, seed=0):
    """Median radiance of up to n valid pixels spread over the scene (for seams and the shift fit)."""
    rng = np.random.default_rng(seed)
    rows = np.sort(rng.choice(rd.nl, min(rd.nl, 200), replace=False))
    picks = []
    for r in rows:
        L, _, ok = rd.block(r, r + 1)
        if ok.any():
            picks.append(L[ok][rng.choice(ok.sum(), min(ok.sum(), max(1, n // len(rows))), replace=False)])
    return np.median(np.concatenate(picks), 0)


# ---------------------------------------------------------------- whole scene
def run(radiance, out_prefix, scene_kw, obs=None, rows=64, log=print, lines=None, bbl_from=None):
    rd = open_radiance(radiance, obs)
    nl, nb, ns = rd.nl, rd.nb, rd.ns
    r_lo, r_hi = (0, nl) if not lines else (max(0, lines[0]), min(nl, lines[1]))
    o = rd.sample_obs()
    sza = float(np.median(o[:, 4]))
    vmax = float(np.max(o[:, 2]))
    vza_nodes = np.linspace(0, max(vmax, 1.0) + 0.5, 3)
    dphi_nodes = np.array([0.0, 90.0, 180.0])
    log(f"{rd.name}: {nl} x {ns} px, {nb} bands, rows {r_lo}-{r_hi}; sza {sza:.2f}, "
        f"view zenith up to {vmax:.1f}, Earth-Sun {rd.esd:.5f} AU")
    L_med = median_radiance(rd)
    k_med = np.pi * rd.esd ** 2 / np.cos(np.radians(sza))
    sc = Scene(rd.centers, rd.fwhm, sza=sza, esd=rd.esd, vza_nodes=vza_nodes, dphi_nodes=dphi_nodes,
               L_med=k_med * L_med * 10.0, log=log, **scene_kw)
    good = (rd.bbl > 0) if rd.bbl is not None else np.ones(nb, bool)
    good &= (rd.centers > 395)
    seams = instrument.seam_check(rd.centers, rd.fwhm, L_med, sc.E_band)
    log("seams: " + ("none (single spectrometer)" if not seams else
                     "; ".join(f"{s['ends_nm']:.0f}|{s['starts_nm']:.0f} ratio {s['ratio']:.4f}"
                               f"{' FLAG' if s['flagged'] else ''}{' low-signal' if s['low_signal'] else ''}"
                               for s in seams)))

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    out_img = Path(str(out_prefix) + "_rfl_img")
    wat_img = Path(str(out_prefix) + "_h2o_img")
    n_out = r_hi - r_lo
    rfl = np.memmap(out_img, "<f4", "w+", shape=(n_out, nb, ns))
    wat = np.memmap(wat_img, "<f4", "w+", shape=(n_out, 3, ns))
    t0 = time.time()
    stats = []
    for r0 in range(r_lo, r_hi, rows):
        r1 = min(r_hi, r0 + rows)
        L, g, ok = rd.block(r0, r1)
        R = np.full((len(L), nb), -9999.0, np.float32)
        W = np.full((len(L), 3), -9999.0, np.float32)
        if ok.any():
            rho, w940, w1140, rms = correct_block(sc, L[ok], g[ok, 2], relative_azimuth(g[ok, 1], g[ok, 3]),
                                                  g[ok, 4], good)
            R[ok] = rho
            W[ok] = np.c_[w940 * sc.cm1, w1140 * sc.cm1, rms]
            stats.append(np.c_[w940, w1140])
        rfl[r0 - r_lo:r1 - r_lo] = R.reshape(r1 - r0, ns, nb).transpose(0, 2, 1)
        wat[r0 - r_lo:r1 - r_lo] = W.reshape(r1 - r0, ns, 3).transpose(0, 2, 1)
        if ((r0 - r_lo) // rows) % 16 == 0:
            el = time.time() - t0
            log(f"rows {r1 - r_lo}/{n_out}  {el / 60:.1f} min  eta {el / (r1 - r_lo) * (r_hi - r1) / 60:.1f} min")
    rfl.flush(); wat.flush()
    s = np.concatenate(stats)
    rel = (s[:, 1] - s[:, 0]) / s[:, 0]

    # bad-band list: provider > this scene > conservative default
    if rd.bbl is not None:
        bbl, bbl_src = rd.bbl, "provider (radiance header)"
    elif bbl_from:
        _, h2, a2 = read_header(bbl_from)
        assert len(a2("wavelength")) == nb and np.allclose(a2("wavelength"), rd.centers, atol=0.5), \
            "bbl source has a different band table"
        bbl, bbl_src = a2("bbl"), f"provider ({Path(bbl_from).name})"
    else:
        T = sc.transmittance(float(np.median(s[:, 0])))
        opaque = T < OPAQUE_T
        if opaque.any():
            bbl, bbl_src = (~opaque).astype(float), f"scene (two-way gas transmittance < {OPAQUE_T:.0%} at median water)"
        else:
            bbl, bbl_src = (~core_opaque(rd.centers)).astype(float), "conservative default"
    log(f"bbl: {int((bbl == 0).sum())} bands off, from {bbl_src}")

    # map info of a row subset: shift the reference pixel
    mi = rd.map_info
    if mi and r_lo:
        p = [x.strip() for x in mi.strip("{}").split(",")]
        p[2] = str(float(p[2]) - r_lo)
        mi = "{" + ", ".join(p) + "}"
    write_header(str(out_img) + ".hdr", ns, n_out, "DUSTRAN surface reflectance (Lambertian, plane-parallel)",
                 mi, rd.centers, rd.fwhm, bbl, coord=rd.coord)
    write_header(str(wat_img) + ".hdr", ns, n_out, "DUSTRAN water column", mi,
                 band_names=["water 940 (cm)", "water 1140 (cm)", "940 fit rms (%)"], coord=rd.coord)
    summary = dict(
        radiance=str(radiance), reflectance=str(out_img) + ".hdr", water=str(wat_img) + ".hdr",
        rows=[r_lo, r_hi], sza=sza, esd=rd.esd, vza_nodes=vza_nodes.tolist(), afgl_water_cm=sc.cm1,
        shift=(sc.shift_fit if sc.shift_fit else
               dict(zip(["A", "B"], map(float, scene_kw.get("shift", (0.389, 0.063e-3)))))),
        bbl_source=bbl_src, bbl_off_nm=[float(x) for x in rd.centers[bbl == 0]],
        water_940_cm_pct=dict(zip(["p5", "p50", "p95"], (np.percentile(s[:, 0], [5, 50, 95]) * sc.cm1).tolist())),
        water_1140_cm_pct=dict(zip(["p5", "p50", "p95"], (np.percentile(s[:, 1], [5, 50, 95]) * sc.cm1).tolist())),
        water_1140_vs_940_median_rel=float(np.median(rel)),
        water_flag_frac_gt10pct=float(np.mean(np.abs(rel) > 0.10)),
        water_flag=bool(abs(np.median(rel)) > 0.10),
        seams=[{k: v for k, v in x.items() if k != "pairs"} for x in seams],
        seam_flagged=any(x["flagged"] for x in seams),
        minutes=(time.time() - t0) / 60,
        **{k: (v if isinstance(v, (int, float, str, list)) else str(v)) for k, v in scene_kw.items()
           if k not in ("aer", "shift")})
    Path(str(out_prefix) + "_summary.json").write_text(json.dumps(summary, indent=1))
    log(json.dumps({k: summary[k] for k in ("water_940_cm_pct", "water_1140_cm_pct",
                                             "water_1140_vs_940_median_rel", "water_flag", "minutes")}))
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("radiance", help="AVIRIS-NG *_img.hdr, or an AVIRIS Classic rdn folder")
    ap.add_argument("out_prefix")
    ap.add_argument("--obs", help="AVIRIS-NG *_obs_ort.hdr")
    ap.add_argument("--taudirs", required=True, help="comma-separated line-by-line tau folders")
    ap.add_argument("--ground", type=float, required=True, help="ground elevation [km] (gas tau run)")
    ap.add_argument("--sensor", type=float, required=True, help="sensor altitude [km] (gas tau run)")
    ap.add_argument("--year", type=int, default=2020, help="flight year (CO2)")
    ap.add_argument("--o3-du", type=float, default=None, help="measured ozone column [DU]")
    ap.add_argument("--pressure-ratio", type=float, default=1.0, help="scale the air column (Rayleigh, CIA)")
    ap.add_argument("--shift", nargs="+", default=["0.389", "0.063e-3"],
                    help="band shift A B [nm, + B (nm - 1500)], or 'fit' to fit per spectrometer")
    ap.add_argument("--aeronet", help="AERONET folder (SIZ_hyb.txt, RIN_hyb.txt, aod file)")
    ap.add_argument("--aeronet-aod-file", default="aod.txt")
    ap.add_argument("--inv-time"); ap.add_argument("--sun-time"); ap.add_argument("--date")
    ap.add_argument("--aod", type=float, default=None, help="AOD at 550 nm (lognormal dust, Angstrom below)")
    ap.add_argument("--angstrom", type=float, default=1.2)
    ap.add_argument("--rows", type=int, default=64)
    ap.add_argument("--bbl-from", help="header whose bad-band list (bbl) to copy, e.g. the provider's reflectance product")
    ap.add_argument("--lines", nargs=2, type=int, default=None, help="process only rows R0..R1 (testing)")
    a = ap.parse_args(argv)
    shift = "fit" if a.shift[0] == "fit" else tuple(map(float, a.shift))
    kw = dict(taudirs=a.taudirs.split(","), ground_km=a.ground, sensor_km=a.sensor, o3_du=a.o3_du,
              pressure_ratio=a.pressure_ratio, shift=shift, year=a.year)
    if a.aeronet:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
        import aeronet_aerosol
        aer, meas = aeronet_aerosol.load(a.inv_time, a.sun_time, a.date, a.aeronet, a.aeronet_aod_file)
        kw.update(aer=aer, aod=meas[500], aod_ref=500.0)
        print(f"AERONET: AOD500 {meas[500]:.3f}, photometer water {meas['PW']:.3f} cm")
    elif a.aod:
        kw.update(aer=aerosol.Dust.fit_angstrom(a.angstrom), aod=a.aod, aod_ref=550.0)
    return run(a.radiance, a.out_prefix, kw, obs=a.obs, rows=a.rows, log=lambda s: print(s, flush=True),
               lines=a.lines, bbl_from=a.bbl_from)


if __name__ == "__main__":
    main()
