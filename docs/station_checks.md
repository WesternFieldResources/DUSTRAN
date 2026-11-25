# Station checks: independent atmosphere measurements for each scene

Every disciplined run compares DUSTRAN's own retrieved or assumed atmosphere with a measurement
that does not come from the spectrometer. These numbers are **reported, never fed back**: no
setting is tuned to match a station. A disagreement goes into `KNOWN_ISSUES.md` and helps decide
whether an issue is atmosphere (class 1), processing (class 2), or instrument (class 3).

Each source is used only for the specific question it can answer.

## Sources: the question each answers, and the answer it gives

| Source | Question | Answer | Limits |
|---|---|---|---|
| **AERONET sun photometer** (NASA GSFC, v3, Level 1.5 direct sun) | Is our assumed aerosol (dust AOD 0.05 at 550 nm, Angstrom 1.2) right for this day? Is our 940 nm water in the right range? | AOD at 500 nm, the 440–870 nm Angstrom exponent, and column precipitable water (cm). We use the points within ±60 min of the overpass. | A column at one point, often 100+ km from the scene and at a different elevation. Water scales strongly with ground height, so we compare after scaling for the height difference, and we report both raw and scaled. |
| **Radiosonde** (NWS, via the University of Wyoming archive) | Is the *shape* of our water profile (standard atmosphere, scaled) wrong? Do 940 and 1140 nm disagree because water sits at a different height? | Water vapour mixing ratio and temperature against height, plus total column water above the scene's ground height. | Launches at 00 and 12 UTC, so 5–7 h from a typical overpass. Launch sites are 50–200 km away. Good for profile shape, weaker for the exact amount. |
| **GPS precipitable water** (ground GNSS stations) | Is our per-pixel water right at the overpass time? | Column water every 5–30 min at a fixed antenna. | Needs surface pressure and temperature to convert delay into water. Not yet checked for these scenes (open). |

## What a run reports

One `stations` block in the run's `_summary.json` and in the pipeline log. The layout below is
an example of the format only. The station values are the real Goldstone and Las Vegas readings
for 3 May 2013; the DUSTRAN column is placeholder until the run exists:

```
station        source     dist_km  dt_min  quantity        station  dustran  diff
Goldstone      AERONET    120      -1      AOD 500 nm      0.051    0.05*    +0.001   (*assumed)
Goldstone      AERONET    120      -1      water cm        0.27     0.31     +0.04    (height-scaled)
Las Vegas      radiosonde 80       -6h     water above gnd 0.29     0.31     +0.02
```

## Availability for the planned scenes (checked 2026-11-24)

| Scene | Overpass (UTC) | AERONET within ±60 min | Radiosonde |
|---|---|---|---|
| Mountain Pass, 3 May 2013 (f130503t01p00r10) | 18:25 | Goldstone, 120 km, 38 points; nearest −1 min: AOD500 0.051, PW 0.27 cm, Angstrom 0.95 | Las Vegas 72388 |
| Mountain Pass, 7 Nov 2014 (f141107t01p00r10) | 19:06 | Goldstone, 121 km, 10 points; nearest −5 min: AOD500 0.015, PW 0.51 cm, Angstrom 0.55 | Las Vegas 72388 |
| Mountain Pass and Ivanpah, 6 Sep 2018 (f180906t01p00r07, r12, r13) | 16:58–18:01 | Goldstone, 120–141 km, 33 points; nearest −3 min: AOD500 0.029, PW 1.27 cm, Angstrom 1.09 | Las Vegas 72388: **confirmed**, 12 UTC profile present |
| Salton Sea, 21 Oct 2015 (f151021t01p00r05) | 16:27 | **None within ±60 min.** Yuma (151 km) only 4 h later: AOD500 0.045, PW 2.50 cm. The Salton_Sea site had no data that day. | San Diego 72293 (coastal, so weak for this site) |
| Cerbat North, 6 Dec 2017 (f171206t01p00r07) | 18:10 | **None.** NEON_Ivanpah (71 km) had no data; Goldstone had no Level 1.5 data that day. | Las Vegas 72388 or Flagstaff 72376 |
| Railroad Valley, 13 Aug 2019 (f190813t01p00r07) | 22:12 | Railroad Valley sites, 14 km (already on disk) | Las Vegas 72388 |

The NEON_Ivanpah AERONET site, 13–20 km from Mountain Pass and Ivanpah, had no data on any
of these dates.

### Status

- AERONET: checked for every scene. The day files are small text, kept with the run.
- Radiosondes: the Las Vegas 2018 profile is confirmed. The others are listed, not yet
  fetched.
- GPS water: not yet checked.
- The `stations` block: built (`dustran/stations.py`). It runs automatically after every
  correction in the one-command Classic wrapper. Tested on Cerbat North, where it reproduces the
  hand check: no site had points within ±60 min.
