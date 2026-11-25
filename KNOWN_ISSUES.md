# Known issues

Open problems found on real scenes, sorted into three classes. "Faithful" means faithful to the
measured radiance: instrument errors are reported, never patched inside the atmosphere step.

How to tell the classes apart: run the same site on several dates.
- **Atmosphere** issues follow the conditions (water, sun angle, season, haze).
- **Processing** issues are the same on every date; they come from our choices.
- **Instrument** issues follow the flight year or the sensor's state, not the weather.

An issue whose class is not yet known is marked *(class?)*.

## 1. Atmosphere physics (water, gases, aerosol, scattering)

- **940 vs 1140 nm water disagree.** On Cerbat North (f171206t01p00r07), 940 nm gives 0.45 cm
  and 1140 nm gives 0.56 cm, which is flagged. The cause could be line data, the water
  profile shape, or a band shift near 1140 nm. *(class?)*
- **No per-scene aerosol.** We assume dust with AOD 0.05 and Angstrom 1.2 unless
  sun-photometer data are given.
- **Blue level 5–7% below the other product at 400–500 nm** (Cuprite 2, AVIRIS-NG). The cause
  could be aerosol, Rayleigh, ozone, or radiance calibration. *(class?)*

## 2. Processing finesse and global coefficients

- **The band-shift fit can lock onto the wrong feature.** On Cerbat North, the CO2 2060 nm
  window pulled the SWIR spectrometer to +2.67 nm. That created a false 2165 nm dip and erased
  the real 2200 nm Al-OH band.
  - Fixed by fitting only the O2 760 and CO2 1600 windows.
  - Still open: spectrometers A and D have no window and take a neighbour's shift.
- **Global settings that are fixed, not fitted:**
  - the scene opacity threshold for bad bands (two-way T < 0.01);
  - shift search range ±4 nm in 0.1 nm steps;
  - the scattering geometry grid (3 view angles × 3 azimuths);
  - the water-table range (0.1–2.5 × standard).
- **Classifier cutoff 0.30** (in bandcheck, not DUSTRAN) decides how much ground is "no
  winner". NASA's reflectance and ours are scored at the same cutoff.

## 3. Instrument (AVIRIS Classic and AVIRIS-NG; varies by year and condition)

- **Per-spectrometer wavelength shifts.** Measured on Cerbat North: B −0.23 nm, C −0.03 nm.
  These are expected to change from year to year.
- **Seam ratios.** Reported on every run and never corrected. The 1872 | 1866 nm seam is
  low-signal, so its ratio is noisy and doesn't show a calibration error.
- **Radiance calibration and atmosphere can't be separated yet.** Any per-band radiance error
  passes straight through to reflectance. If we ever want a fix, it belongs in a separate
  radiance step, labelled and switchable.
- **Shallower 900 nm iron band than the other product** (Cuprite 2, AVIRIS-NG). The ratio
  rises to 1.02 at 840–960 nm. The water fit has been ruled out. NASA's empirical band gains
  are one candidate. *(class?)*
- **Short cross-track dashes in the Ivanpah 2018 classifier map** (f180906t01p00r07). They are a
  few pixels tall, cross the mapped units, and are densest in the north half and along the
  bottom edge. Judged instrument from the map: too uniform to be processing. Not yet confirmed
  in the radiance rows.
- **Year-to-year differences on the same ground, present in NASA's own reflectance too**
  (Mountain Pass: Nov 2014, May 2014, Sep 2018). At similar sun angles, May 2014 to Sep 2018
  gives a level ratio of 0.84 (750–1600 nm) and 0.67 at 1997 nm. Gas-band ripple is 1.3–4.4%
  in 2014 and 7.9–8.2% in 2018. Within 2014 the level ratio is 1.06, so the drop follows the
  year, not the sun. The 2080–2290 nm shape differs by 2–4% rms in every pair, including
  the same-year one. NASA's band tables are identical in every year from 2014 to 2021.
  Details: `C:\WFR\aviris\areas\nasa_mp_2014_vs_2018\crossover_results.md`.
- **Along-track stripes in the classifier maps** (Cuprite 2, AVIRIS-NG). They appear with both
  corrections. Methane has been ruled out; they are most likely detector columns.

## Not yet assigned

- **Strong alteration signals on Cerbat North:** high sulfidation 4.5% lit, low sulfidation
  12.1% lit. Not yet checked against the diagnostic bands in the spectra, or against whether
  the pattern follows the geology or the instrument.

## Planned multi-date test

The same site flown in different years and conditions, each with the provider's reflectance:

| Line | Site | Date | Sun elevation |
|---|---|---|---|
| ~~f130503t01p00r10~~ | Mountain Pass | May 2013 | 64° (skipped: the provider's package has no `.gain` or `.spc` file, so integers can't be converted to radiance without guessing) |
| f141107t01p00r10 | Mountain Pass | Nov 2014 | 38° |
| f180906t01p00r12 | Mountain Pass | Sep 2018 | 51° |
| f180906t01p00r13 | Mountain Pass 2 | Sep 2018 | — |
| f180906t01p00r07 | Ivanpah playa | Sep 2018 | — |
| f151021t01p00r05 | Salton Sea | Oct 2015 | 28°, below sea level, humid |

### Results (2026-11-24)

| Run | Sun zen | Shift B / C (nm) | Seams VNIR, SWIR1 | 1872 seam | Water 940 / 1140 (cm) | 1140 vs 940 | Station water, scaled to scene | Station AOD500 (assumed 0.056) |
|---|---|---|---|---|---|---|---|---|
| Cerbat N, Dec 2017 | 61° | −0.23 / −0.03 | 0.992, 1.008 | 1.016 low | 0.45 / 0.56 | +25% | none | none |
| Mountain Pass, Nov 2014 | 52° | −0.61 / +0.20 | 0.990, 1.001 | 0.76 low | 0.64 / 0.80 | +25% | 0.49 (Goldstone, 121 km) | 0.014 |
| Mountain Pass, Sep 2018 | 39° | −0.94 / −0.41 | 0.993, 1.002 | 0.51 low | 1.90 / 2.01 | +6% | 1.20 (Goldstone, 121 km) | 0.030 |
| Mountain Pass 2, Sep 2018 | 37° | −0.92 / −0.46 | 0.992, 1.002 | 0.51 low | 1.91 / 2.09 | +9% | 1.25 (Goldstone, 108 km) | 0.032 |
| Ivanpah, Sep 2018 | 47° | −0.96 / −0.43 | 0.992, 1.005 | 0.51 low | 2.08 / 2.37 | +13% | 1.41 (Goldstone, 137 km) | 0.029 |
| Salton Sea, Oct 2015 | 62° | −0.71 / +0.02 | 0.992, 1.005 | n/a | 2.21 / 3.02 | +39% | 1.64 (Table Mtn, 195 km, 2.2 km up) | 0.025 |

Patterns, first sort:
- **Class 3 (instrument):** the shifts track the flight. The three 2018 lines agree within
  0.04 nm, while the 2014, 2015 and 2017 lines differ by up to 0.7 nm. The VNIR seam is
  ~0.992 on every line.
- **Classes 1 and 2 (atmosphere and processing):** 940 nm water is 30–60% above the
  height-scaled station value on every line that has one. 1140 nm runs above 940 nm on every
  line. The gap is largest at high sun zenith (52–62°) and smallest near 37–39°, which points at
  path geometry or the water profile shape more than the amount of water.
- **Our fixed AOD of 0.056 at 500 nm is 2–4 times every station measurement.**
- **The 1872 nm seam ratio is meaningless at this signal level.**

How to read the results:
- An issue that tracks the year belongs to class 3 (instrument).
- One that tracks the water or the sun angle belongs to class 1 (atmosphere).
- One that never moves belongs to class 2 (processing).
