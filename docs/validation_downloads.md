# Validation downloads: AVIRIS Classic + same-day Sentinel-2 (+ RadCalNet where it exists)

Found 2026-11-23 from the NASA flight catalogues (footprint contains the target) and the
Copernicus catalogue. Ivanpah Playa is a NASA calibration site that Classic flew almost every
year. Save everything to `D:\Downloads`.

## Priority 1: Classic + same-day Sentinel-2 (Ivanpah Playa, 35.570 N, 115.397 W)

| Date | AVIRIS Classic (one line over the playa) | Sentinel-2 L2A same day |
|---|---|---|
| 2016-06-15 | `https://popo.jpl.nasa.gov/avcl/y16_data/f160615t01p00r08.tar.gz` (0.9 GB) | `S2A_MSIL2A_20160615T183312_N0500_R127_T11SPV_20231008T090400` |
| 2017-12-07 | `https://popo.jpl.nasa.gov/avcl/y17_data/f171207t01p00r07.tar.gz` (0.8 GB) | `S2A_MSIL2A_20171207T183301_N0500_R127_T11SPV_20230718T164104` |
| 2021-02-24 | `f210224t01p00r07` (FTP link in the Classic catalogue, 0.5 GB) | `S2B_MSIL2A_20210224T182329_N0500_R127_T11SPV_20230531T060733` |

## Priority 2: Classic + Sentinel-2 within 1–2 days

| Date | AVIRIS Classic | Sentinel-2 L2A |
|---|---|---|
| 2018-09-06 | `https://popo.jpl.nasa.gov/avcl/y18_data/f180906t01p00r07.tar.gz` (1.0 GB) | `S2B_MSIL2A_20180908T181909_N0500_R127_T11SPV_20230824T183102` (+2 days) |
| 2020-10-05 | `https://popo.jpl.nasa.gov/avcl/y20_data/f201005t01p00r09.tar.gz` (0.9 GB) | `S2B_MSIL2A_20201007T182239_N0500_R127_T11SPV_20230326T194159` (+2 days) |

## Priority 3: Railroad Valley (has RadCalNet)

| Date | AVIRIS | Sentinel-2 | RadCalNet |
|---|---|---|---|
| 2019-08-13 | Classic `f190813t01p00r08.tar.gz` (second line, 2.8 GB, y19_data) | already have 2019-08-12 | RVUS 2019 day 224 (to match Sentinel's day) |
| 2020-07-12 | NG `ang20200712t215355.tar.gz` (second line, y20_data) | already have | already have |

## Notes
- Several other Classic lines exist for each Ivanpah date (r06–r19). The one listed is enough;
  more lines on the same date are a cheap consistency check later.
- Classic Ivanpah flights are 16–19 UTC and Sentinel-2 passes at ~18:20–18:33 UTC, so time of
  day is much closer than at Railroad Valley.
- AERONET has a site `NEON_Ivanpah` next to the playa (dates not checked yet). I can fetch that
  myself.
- No Classic + Sentinel-2 match within 2 days for 2015-10-26 or 2017-06-28 at Ivanpah, or for
  the Rogers Dry Lake and Black Rock flights.
