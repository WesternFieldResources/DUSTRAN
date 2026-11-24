"""Water vapour self-continuum in the 0.94 and 1.14 um bands from a water-dimer model.

Model (Simonova, Ptashnik et al., JQSRT, "Water vapour self-continuum in near-visible IR absorption
bands: measurements and semiempirical model of water dimer absorption", eqs. 3 and 5, Table 2,
Table 3). Written for DUSTRAN from the published equations and parameters:

    Cs(nu, T) = Kb(T) * sum_i Sb_i fb(nu - nu_i, 20 cm-1) + Kq(T) * sum_j Sq_j fq(nu - nu_j, 10 cm-1)

    Cs        cross-section per water molecule per atm of water [cm2 molec-1 atm-1]
    bound (b) dimers: the 14 overtone transitions of Table 2, HWHM 20 cm-1
    quasibound (q) dimers: HITRAN water lines at doubled intensity, HWHM 10 cm-1
    total equilibrium constant (eq. 5): Kb + Kq = 3.717e8 * T^-3.886  [atm-1]
    split between b and q: the paper's Kb(aver) at 398 K (0.0031) and 431 K (0.0019) gives an
        effective binding temperature of 2540 K; Kb(T) = 0.0031 exp(2540 (1/T - 1/398)).
        That extrapolation (and so the b/q split at 260-300 K) is OUR assumption, flagged.
Profiles: Lorentzian (Doppler width is negligible at 10-20 cm-1 HWHM), cut at +-250 cm-1.
Absorption coefficient: k [cm-1] = Cs * n_w [molec cm-3] * P_w [atm]; tau = Cs * P_w * column_w.
This is a self-continuum only; the foreign continuum is not in this model.
"""
import numpy as np

C2 = 1.4387769               # cm K
BOUND = np.array([            # Table 2: position [cm-1], intensity [cm/molec]
    (8530.5, 2.35e-21), (8754.9, 6e-22), (8804.8, 7e-21), (8806.9, 5.5e-20), (8930.1, 1.9e-21),
    (9006.9, 4.85e-24), (10057.5, 2.65e-22), (10161.1, 1.95e-21), (10601.0, 1.85e-21),
    (10611.0, 6.5e-21), (10615.3, 1.8e-20), (10673.7, 9.5e-21), (10869.7, 8e-22),
    (10889.1, 3.5e-21), (11042.0, 2.25e-21)])
HW_B, HW_Q, CUT = 20.0, 10.0, 250.0


def k_total(T):
    return 3.717e8 * np.asarray(T, float) ** -3.886


def k_bound(T):
    return 0.0031 * np.exp(2540.0 * (1.0 / np.asarray(T, float) - 1.0 / 398.0))


def _lorentz(d, hw):
    return hw / np.pi / (d * d + hw * hw)


def water_lines(wn_min=7800.0, wn_max=11600.0):
    """HITRAN H2O lines (position, S at 296 K, lower-state energy) via RADIS' local databank."""
    from radis import SpectrumFactory
    sf = SpectrumFactory(wavenum_min=wn_min, wavenum_max=wn_max, molecule="H2O", isotope="1",
                         wstep=0.01, truncation=25, verbose=0)
    sf.fetch_databank("hitran")
    df = sf.df0
    return np.asarray(df["wav"]), np.asarray(df["int"]), np.asarray(df["El"])


def line_strength(S296, El, nu, T):
    """Monomer line intensity at T (partition function approximated as T^1.5 for H2O)."""
    return (S296 * (296.0 / T) ** 1.5 * np.exp(-C2 * El * (1 / T - 1 / 296.0))
            * (1 - np.exp(-C2 * nu / T)) / (1 - np.exp(-C2 * nu / 296.0)))


def cross_section(wn, T, lines):
    """Cs(wn, T) [cm2 molec-1 atm-1] on grid wn (coarse, ~1 cm-1 is enough)."""
    nu, S296, El = lines
    m = S296 > 1e-27                              # weak lines add nothing at 10 cm-1 widths
    nu, S = nu[m], 2.0 * line_strength(S296[m], El[m], nu[m], T)
    q = np.zeros(len(wn))
    for k in range(len(wn)):
        d = wn[k] - nu
        sel = np.abs(d) < CUT
        q[k] = np.sum(S[sel] * _lorentz(d[sel], HW_Q))
    b = np.zeros(len(wn))
    for nu_i, s in BOUND:
        d = wn - nu_i
        b += np.where(np.abs(d) < CUT, s * _lorentz(d, HW_B), 0.0)
    kb = k_bound(T); kq = max(k_total(T) - kb, 0.0)
    return kb * b + kq * q


def optical_depth(L, wn, groups, lines=None, coarse_step=1.0):
    """Self-continuum optical depth (water at the layer table's amounts) summed over layer groups."""
    lines = lines if lines is not None else water_lines()
    grid = np.arange(wn.min() - 1, wn.max() + coarse_step + 1, coarse_step)
    out = {g: np.zeros(len(grid)) for g in groups}
    kB = 1.380649e-16
    for i, T in enumerate(L["T"]):
        nw = L["col"]["H2O"][i] / L["dz"][i]
        pw_atm = nw * kB * T / 1.01325e6
        tau = cross_section(grid, T, lines) * pw_atm * L["col"]["H2O"][i]
        for g, m in groups.items():
            if m[i]:
                out[g] += tau
    return {g: np.interp(wn, grid, v) for g, v in out.items()}
