# The DUSTRAN stack

DUSTRAN turns imaging-spectrometer radiance (AVIRIS-NG, AVIRIS Classic) into surface reflectance
over desert scenes. It runs in plain Python/numpy with no containers and no stored lookup tables,
and every step can be inspected.

Status key: ✅ working · 🟡 partial · ⬜ planned

## What each layer does

| Layer | Data in | Engine | Module | Status |
|---|---|---|---|---|
| Atmosphere | AFGL standard atmospheres (Anderson et al. 1986); CO2 scaled to the flight year | numpy | `atmosphere.py` | 🟡 |
| Gas absorption | HITRAN line lists | RADIS Voigt cross sections per layer p, T | `gas.py` | ✅ |
| O2 collision-induced absorption | HITRAN CIA (O2–O2, O2–N2) | numpy | `cia.py` | ✅ |
| Ozone | Serdyuchenko et al. 2014 cross sections | numpy | `xsec.py` | ✅ |
| Water continuum | published continuum physics + measurements | numpy | `dimer.py` (self, dimers) | 🟡 see RESULTS.md |
| Solar | TSIS-1 HSRS (Coddington et al. 2021) | table, convolved with gases on the fine grid | `experiments/playa_gas_score.py` | ✅ |
| Rayleigh | Bodhaine et al. 1999 | closed form | `rayleigh.py` | ✅ |
| Dust | lognormal or measured (AERONET) sizes, provisional refractive index | our Mie code | `mie.py`, `aerosol.py` | 🟡 |
| Scattering | layers + Rayleigh + dust | our successive-orders solver | `sos.py` | ✅ scalar |
| Instrument | band centres/FWHM per flight | Gaussian bands + fitted shift; seam check | `instrument.py` | ✅ |
| Inversion | radiance + all of the above | closed form + small fit | to write (per-pixel / batched) | ⬜ |

### How they fit together

```
per scene (once):
  AFGL profile ──atmosphere.py──► layers (p, T, columns)
  layers + HITRAN ──gas.py──► τ_gas(ν) per gas, 0.002 cm⁻¹ grid      (k-distribution per band: planned)
  layers + HITRAN CIA + ozone ──cia.py, xsec.py──► τ(ν)
  layers + Rayleigh + dust ──sos.py──► R_path, T↓, T↑, S on a coarse grid
  TSIS ──► E_sun(ν) on the fine grid
  instrument.py ──► band responses (with fitted shift); seam check reported

per pixel:
  water w (fit at 940/1140) ──► band-averaged E_sun · T_gas
  ρ = (r − R_path) / (T↓T↑ + S (r − R_path)),   r = gas-corrected apparent reflectance
```

## The physics in one line

At-sensor radiance for surface reflectance ρ (Chandrasekhar 1960; Tanré et al. 1979):

```
L = L_path + (E_sun · cos θs · T↓ · T↑ / π) · ρ / (1 − S · ρ)
```

- `L_path`: light scattered into the sensor without reaching the ground
- `T↓`, `T↑`: transmittance down to the ground and back up
- `S`: spherical albedo, the fraction of ground-reflected light the sky sends back down

**The key design split** (Gao et al. 1993): gas absorption is sharp and well measured, so it is
computed line by line at 0.002 cm⁻¹. Scattering is smooth and its inputs (the dust) are
uncertain, so it is computed coarsely. The two are multiplied together.

## Layer notes

### Atmosphere: `dustran/atmosphere.py`
- Level profile → layers with pressure, temperature and per-gas columns, truncated at the
  surface and split at the sensor altitude. Each column is integrated exactly for a density
  that falls exponentially between levels.
- The AFGL tables carry 1976 CO2 (330 ppm); `co2_ppm(year)` scales it to the flight year.
- **To do:** per-scene profiles (radiosonde or reanalysis). The water profile shape is a prime
  suspect for the remaining water bias (RESULTS.md).

### Gas absorption: `dustran/gas.py`
- HITRAN line list → Voigt cross section at each layer's p and T → × layer column, per gas, so
  `τ_H2O(w) = w · τ_H2O(1)`. Lines are cut 25 cm⁻¹ from their centres.
- Grid 0.002 cm⁻¹, computed in 200 cm⁻¹ chunks; once per site and altitude (15–30 min).
- **To do:** per-band k-distribution, so per-pixel water fitting is fast.

### Water continuum
- Self-continuum (water–water, dominated by dimers): `dimer.py`, from Simonova et al. It
  scales with water squared, so it is negligible over dry deserts and matters on humid scenes.
- Foreign continuum (water–N2/O2): needed. On both playa test days a continuum of the expected
  size moves fitted water 10–15 % toward the photometer. Line wings beyond 25 cm⁻¹ supply only a
  seventh of it. A collision band under each line (`foreign.py`), measured on one playa day,
  does not transfer to the other (RESULTS.md §4), so it is not adopted. Next: per-scene water
  profiles, then retest.

### Scattering: `dustran/sos.py`
- Successive orders of scattering, Fourier azimuth modes, Gauss streams, delta-M; outputs
  R_path, T↓, T↑, S. Scalar (no polarisation). Seven physics tests in `tests/test_sos.py`.

### Instrument: `dustran/instrument.py`
- Gaussian bands with a wavelength shift. AVIRIS-NG over Railroad Valley: +0.30 nm in the
  near-infrared rising to +0.50 nm at 2330 nm.
- **Seam check** (`seam_check`): where spectrometers overlap, their bands see the same light, so
  the ratio gives the relative calibration from the scene itself.

## Design rules

- **Seam check in every workflow.** Report the seam ratios for every multi-spectrometer scene
  and flag anything over 1 %. Never correct silently: a correction must be consistent or
  calculable, and agreed first.
- No containers. No stored lookup tables. Small per-scene grids built on the fly are fine.
- Every module written from papers and first principles; every number traceable to a paper or
  a public dataset (`BIBLIOGRAPHY.md`, `PROVENANCE.md`).
- Judge every piece by measurements on real scenes. Real observations only, no synthetic test
  spectra.

## Build order

1. ✅ Line-by-line gases, CIA, ozone, wavelength shift (AVIRIS-NG, Railroad Valley 2020).
2. ✅ Rayleigh, dust, scattering solver; validated against RadCalNet.
3. ✅ AVIRIS Classic: seam check, multi-flight brightness test (Railroad Valley, Ivanpah).
4. 🟡 Water: foreign continuum, per-scene profiles.
5. ⬜ k-distribution and batched per-pixel inversion; full-scene runs.
6. ⬜ Classic calibration policy (needs more flights; see RESULTS.md).
