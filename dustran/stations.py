"""Station check: compare a run's atmosphere with independent sun-photometer measurements.

    python -m dustran.stations SCENE_FOLDER SUMMARY_JSON [--km 200] [--minutes 60]

Reported, never fed back: no DUSTRAN setting is changed by these numbers (docs/station_checks.md).

For an AVIRIS Classic line: scene centre and overpass time from its own IGM and obs files; every
AERONET site within --km; Level 1.5 direct-sun points within +-minutes of the overpass. Compared:
  - AOD at 500 nm vs the run's assumed dust AOD (converted from 550 nm with its Angstrom exponent)
  - Angstrom 440-870 vs the assumed Angstrom
  - column water vs the run's 940 nm median, raw and scaled for the height difference with a
    2 km water scale height (stated, not fitted).
The result is added to the summary JSON as "stations" and printed.
"""
import argparse, csv, datetime, io, json, math, urllib.request
from pathlib import Path
import numpy as np

from .batch import read_header

SITES = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt"
DATA = ("https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site={s}&year={y}&month={m}&day={d}"
        "&year2={y}&month2={m}&day2={d}&AOD15=1&AVG=10&if_no_html=1")
WATER_SCALE_KM = 2.0


def _get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def scene_centre(folder):
    """Median lon, lat (IGM) and median UTC hour (obs band 'UTC time') of a Classic line."""
    d = Path(folder)
    igm = next(d.glob("*_ort_igm"))
    _, h, _ = read_header(str(igm) + ".hdr")
    nl, ns = int(h["lines"]), int(h["samples"])
    a = np.memmap(igm, "<f8", "r", shape=(nl, 3, ns))
    lon, lat = np.asarray(a[::7, 0, ::7]).ravel(), np.asarray(a[::7, 1, ::7]).ravel()
    ok = (np.abs(lat) <= 90) & (np.abs(lon) <= 180) & (lat != 0)
    obs = next(d.glob("*_obs_ort"))
    _, ho, _ = read_header(str(obs) + ".hdr")
    o = np.memmap(obs, "<f4", "r", shape=(int(ho["lines"]), int(ho["bands"]), int(ho["samples"])))
    t = np.asarray(o[::7, 9, ::7]).ravel()
    t = t[(t > 0) & (t < 24)]
    return float(np.median(lon[ok])), float(np.median(lat[ok])), float(np.median(t))


def km(la1, lo1, la2, lo2):
    p = math.radians
    a = math.sin(p(la2 - la1) / 2) ** 2 + math.cos(p(la1)) * math.cos(p(la2)) * math.sin(p(lo2 - lo1) / 2) ** 2
    return 6371 * 2 * math.asin(math.sqrt(a))


def check(folder, summary, max_km=200, minutes=60, angstrom=None, log=print):
    s = json.loads(Path(summary).read_text())
    if s.get("aod") is not None:
        s.setdefault("aod_550", s["aod"] * (550 / s.get("aod_ref", 550)) ** (-(angstrom or 0)))
    if angstrom is not None:
        s["angstrom"] = angstrom
    name = Path(folder).name.split("rdn")[0].rstrip("_")
    yy, mm, dd = int(name[1:3]), int(name[3:5]), int(name[5:7])
    date = datetime.date(2000 + yy, mm, dd)
    lon, lat, utc = scene_centre(folder)
    ground = s.get("ground_km")
    aod550, ang = s.get("aod_550"), s.get("angstrom")
    w940 = s["water_940_cm_pct"]["p50"]

    rows = list(csv.reader(io.StringIO(_get(SITES))))[2:]
    near = sorted((km(lat, lon, float(r[2]), float(r[1])), r[0], float(r[3]) / 1000) for r in rows if len(r) > 3)
    near = [n for n in near if n[0] <= max_km]
    out = dict(scene=dict(lon=lon, lat=lat, utc_hour=utc, date=str(date)), source="AERONET v3 Level 1.5 direct sun",
               window_min=minutes, water_scale_km=WATER_SCALE_KM, sites=[])
    for dist, site, z in near:
        txt = _get(DATA.format(s=site, y=date.year, m=date.month, d=date.day))
        lines = txt.splitlines()
        i = next((k for k, l in enumerate(lines) if l.startswith("AERONET_Site")), None)
        if i is None:
            out["sites"].append(dict(site=site, km=round(dist), points=0)); continue
        pts = []
        for r in csv.DictReader(lines[i:]):
            hh, mi, ss = map(int, r["Time(hh:mm:ss)"].split(":"))
            dt = (hh + mi / 60 + ss / 3600 - utc) * 60
            if abs(dt) <= minutes:
                f = lambda k: float(r[k]) if r.get(k) not in (None, "", "-999.", "-999.000000") and float(r[k]) > -900 else np.nan
                pts.append((dt, f("AOD_500nm"), f("440-870_Angstrom_Exponent"), f("Precipitable_Water(cm)")))
        if not pts:
            out["sites"].append(dict(site=site, km=round(dist), points=0)); continue
        p = np.array(pts)
        aod, an, pw = np.nanmedian(p[:, 1]), np.nanmedian(p[:, 2]), np.nanmedian(p[:, 3])
        e = dict(site=site, km=round(dist), site_elev_km=round(z, 3), points=len(p),
                 nearest_min=round(float(p[np.argmin(np.abs(p[:, 0])), 0]), 1),
                 aod500=round(float(aod), 4), angstrom=round(float(an), 3), water_cm=round(float(pw), 3))
        if aod550 is not None and ang is not None:
            e["dustran_aod500_assumed"] = round(aod550 * (500 / 550) ** (-ang), 4)
            e["dustran_angstrom_assumed"] = ang
        if ground is not None:
            e["water_cm_scaled_to_scene"] = round(float(pw * math.exp(-(ground - z) / WATER_SCALE_KM)), 3)
        e["dustran_water940_cm"] = round(w940, 3)
        out["sites"].append(e)
    s["stations"] = out
    Path(summary).write_text(json.dumps(s, indent=1))
    log(f"stations: scene {lat:.3f},{lon:.3f} {date} {utc:.2f} UTC; DUSTRAN water 940 median {w940:.3f} cm")
    for e in out["sites"]:
        if not e["points"]:
            log(f"  {e['site']:<20} {e['km']:>4} km  no points within +-{minutes} min"); continue
        log(f"  {e['site']:<20} {e['km']:>4} km  {e['points']} pts (nearest {e['nearest_min']:+.0f} min)  "
            f"AOD500 {e['aod500']:.3f} (assumed {e.get('dustran_aod500_assumed', float('nan')):.3f})  "
            f"Angstrom {e['angstrom']:.2f} (assumed {e.get('dustran_angstrom_assumed', float('nan'))})  "
            f"water {e['water_cm']:.2f} cm, scaled {e.get('water_cm_scaled_to_scene', float('nan')):.2f} "
            f"(DUSTRAN {e['dustran_water940_cm']:.2f})")
    if not any(e["points"] for e in out["sites"]):
        log(f"  no AERONET site within {max_km} km had points within +-{minutes} min")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("folder"); ap.add_argument("summary")
    ap.add_argument("--km", type=float, default=200); ap.add_argument("--minutes", type=float, default=60)
    ap.add_argument("--angstrom", type=float, default=None, help="the run's assumed Angstrom exponent")
    a = ap.parse_args()
    check(a.folder, a.summary, a.km, a.minutes, a.angstrom)
