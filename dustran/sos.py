"""Successive orders of scattering (SOS): a plane-parallel, scalar radiative transfer solver.

Written for DUSTRAN from first principles. No external RT code.

The atmosphere is a stack of homogeneous sublayers (thin enough that the source function is
nearly constant across each). Radiance is expanded in azimuth Fourier modes,
    I(tau, mu, phi) = sum_m (2 - delta_m0) I_m(tau, mu) cos(m * dphi),
and directions are Gauss-Legendre nodes in each hemisphere, plus any extra "output" directions
(zero weight). Order 1 is light scattered once out of the direct solar beam; order n+1 is light
scattered once out of order n. Each order is propagated exactly through each sublayer for a
constant source. Orders are summed until they stop contributing.

Phase function: P(cos Theta) = sum_l beta_l P_l(cos Theta), normalised so that
(1/4 pi) * integral P dOmega = 1 (beta_0 = 1). Per sublayer: optical depth, single-scattering
albedo, and Legendre moments (Rayleigh + aerosol mixed by scattering optical depth).

Units: incident solar flux on a surface normal to the beam = 1. Reflectance-like quantities are
pi * I / mu0.

Scalar only (polarisation ignored): for Rayleigh this misstates path radiance by a few percent.
Strongly forward-peaked phase functions (aerosol) are handled with delta-M scaling in
surface_terms (Wiscombe 1977); no single-scattering (TMS) correction yet.
"""
import numpy as np
from math import factorial
from scipy.special import lpmv


def _norm_assoc(lmax, m, x):
    """(lmax, len(x)) table of sqrt((l-m)!/(l+m)!) P_l^m(x), no Condon-Shortley phase (0 for l < m)."""
    x = np.asarray(x, float)
    out = np.zeros((lmax, len(x)))
    for l in range(m, lmax):
        out[l] = lpmv(m, l, x) * (-1) ** m * np.sqrt(factorial(l - m) / factorial(l + m))
    return out


def _pm_matrix(beta, m, x, y):
    """P^m(x_i, y_j) = sum_{l>=m} beta_l (l-m)!/(l+m)! P_l^m(x_i) P_l^m(y_j)."""
    L = len(beta)
    A, B = _norm_assoc(L, m, x), _norm_assoc(L, m, y)
    return (A * np.asarray(beta)[:, None]).T @ B


def delta_m(dtau, ssa, beta, n_terms):
    """Delta-M scaling (Wiscombe 1977): move the forward peak into the direct beam.
    beta: (n, L) moments with phase = sum beta_l P_l; keeps n_terms moments.
    Returns scaled dtau, ssa, beta (n, n_terms)."""
    beta = np.atleast_2d(beta)
    L = beta.shape[1]
    l = np.arange(L)
    chi = beta / (2 * l + 1)
    if L <= n_terms:
        return dtau, ssa, beta
    f = chi[:, n_terms]
    chi_s = (chi[:, :n_terms] - f[:, None]) / (1 - f[:, None])
    beta_s = chi_s * (2 * np.arange(n_terms) + 1)
    tau_s = (1 - ssa * f) * dtau
    ssa_s = (1 - f) * ssa / (1 - ssa * f)
    return tau_s, ssa_s, beta_s


class Atmosphere:
    """Sublayers top to bottom: dtau (n,), ssa (n,), beta (n, L) Legendre moments."""

    def __init__(self, dtau, ssa, beta):
        self.dtau = np.asarray(dtau, float)
        self.ssa = np.asarray(ssa, float)
        self.beta = np.atleast_2d(beta)
        self.tau = np.r_[0.0, np.cumsum(self.dtau)]        # level optical depths, top = 0

    @staticmethod
    def from_layers(tau_layers, ssa_layers, beta_layers, max_dtau=0.01, levels=()):
        """Split layers into sublayers of at most `max_dtau`. `levels` are layer-boundary
        indices whose sublayer index is recorded in .level_index (e.g. the sensor level)."""
        d, w, b, idx, k = [], [], [], {}, 0
        for i, (t, s, bb) in enumerate(zip(tau_layers, ssa_layers, beta_layers)):
            if i in levels:
                idx[i] = k
            n = max(1, int(np.ceil(t / max_dtau)))
            d += [t / n] * n; w += [s] * n; b += [bb] * n; k += n
        if len(tau_layers) in levels:
            idx[len(tau_layers)] = k
        a = Atmosphere(d, w, np.array(b))
        a.level_index = idx
        return a


def gauss_nodes(n):
    x, w = np.polynomial.legendre.leggauss(2 * n)
    m = x > 0
    return x[m], w[m]                    # mu in (0,1), weights summing to 1


class SOS:
    def __init__(self, atm, n_streams=8, extra_mu=(), max_order=200, tol=1e-7):
        self.atm = atm
        mu, w = gauss_nodes(n_streams)
        self.mu = np.r_[mu, np.asarray(extra_mu, float)]
        self.w = np.r_[w, np.zeros(len(extra_mu))]
        self.n_gauss = len(mu)
        self.max_order, self.tol = max_order, tol
        self.n_modes = atm.beta.shape[1]

    # ---- one azimuth mode -------------------------------------------------------------
    def _propagate(self, J_up, J_dn, I_top_dn, I_bot_up):
        """Exact transfer through constant-source sublayers.
        J_*: (n_sub, n_mu) sources. Returns I_up, I_dn at levels (n_sub+1, n_mu)."""
        a = self.atm
        n = len(a.dtau)
        t = np.exp(-a.dtau[:, None] / self.mu[None, :])           # sublayer transmission
        I_dn = np.zeros((n + 1, len(self.mu))); I_up = np.zeros_like(I_dn)
        I_dn[0] = I_top_dn
        for k in range(n):
            I_dn[k + 1] = I_dn[k] * t[k] + J_dn[k] * (1 - t[k])
        I_up[n] = I_bot_up
        for k in range(n - 1, -1, -1):
            I_up[k] = I_up[k + 1] * t[k] + J_up[k] * (1 - t[k])
        return I_up, I_dn

    def _scatter(self, m, I_up, I_dn, Pm):
        """Source in each sublayer from scattering the (level-averaged) diffuse field."""
        a = self.atm
        Iu = 0.5 * (I_up[:-1] + I_up[1:]); Id = 0.5 * (I_dn[:-1] + I_dn[1:])
        wu, wd = Iu * self.w, Id * self.w                          # weight by quadrature
        # P^m(out, in) for out/in up (+mu) and down (-mu): Pm[layer_group][...]
        J_up = np.zeros_like(Iu); J_dn = np.zeros_like(Id)
        for g, (rows, PPuu, PPud, PPdu, PPdd) in Pm.items():
            J_up[rows] = 0.5 * a.ssa[rows, None] * (wu[rows] @ PPuu.T + wd[rows] @ PPud.T)
            J_dn[rows] = 0.5 * a.ssa[rows, None] * (wu[rows] @ PPdu.T + wd[rows] @ PPdd.T)
        return J_up, J_dn

    def _phase_tables(self, m, mu0=None):
        """P^m tables per distinct phase function (sublayers grouped by identical beta)."""
        a = self.atm
        keys, inv = np.unique(np.round(a.beta, 12), axis=0, return_inverse=True)
        inv = np.ravel(inv)
        L = a.beta.shape[1]
        Au = _norm_assoc(L, m, self.mu); Ad = _norm_assoc(L, m, -self.mu)
        A0 = _norm_assoc(L, m, [-mu0]) if mu0 is not None else None
        Pm, P0 = {}, {}
        for g, beta in enumerate(keys):
            rows = np.where(inv == g)[0]
            bu, bd = (Au * beta[:, None]).T, (Ad * beta[:, None]).T
            Pm[g] = (rows, bu @ Au, bu @ Ad, bd @ Au, bd @ Ad)
            if mu0 is not None:
                P0[g] = (rows, (bu @ A0)[:, 0], (bd @ A0)[:, 0])
        return Pm, P0

    def _iterate(self, m, J_up, J_dn, Pm, I_top_dn, I_bot_up):
        I_up, I_dn = self._propagate(J_up, J_dn, I_top_dn, I_bot_up)
        tot_up, tot_dn = I_up.copy(), I_dn.copy()
        for order in range(2, self.max_order + 1):
            J_up, J_dn = self._scatter(m, I_up, I_dn, Pm)
            I_up, I_dn = self._propagate(J_up, J_dn, 0.0, 0.0)
            tot_up += I_up; tot_dn += I_dn
            if max(np.abs(I_up).max(), np.abs(I_dn).max()) < self.tol * max(np.abs(tot_up).max(), 1e-30):
                break
        return tot_up, tot_dn

    # ---- problems ---------------------------------------------------------------------
    def solar(self, mu0, n_modes=None):
        """Diffuse radiance modes for a solar beam (flux 1) at cosine mu0, black surface.
        Returns list over m of (I_up, I_dn) at levels, shape (n_levels, n_mu)."""
        a = self.atm
        n_modes = n_modes or self.n_modes
        # sublayer-mean direct beam
        e0 = np.exp(-a.tau[:-1] / mu0) * mu0 / a.dtau * (1 - np.exp(-a.dtau / mu0))
        out = []
        for m in range(n_modes):
            Pm, P0 = self._phase_tables(m, mu0)
            J_up = np.zeros((len(a.dtau), len(self.mu))); J_dn = np.zeros_like(J_up)
            for g, (rows, pu, pd) in P0.items():
                J_up[rows] = a.ssa[rows, None] / (4 * np.pi) * e0[rows, None] * pu[None, :]
                J_dn[rows] = a.ssa[rows, None] / (4 * np.pi) * e0[rows, None] * pd[None, :]
            out.append(self._iterate(m, J_up, J_dn, Pm, 0.0, 0.0))
        return out

    def isotropic_bottom(self):
        """Diffuse field for unit isotropic radiance entering at the bottom (m = 0 only).
        Returns (I_up, I_dn) of the scattered + transmitted field (direct part included in I_up)."""
        Pm, _ = self._phase_tables(0)
        z = np.zeros((len(self.atm.dtau), len(self.mu)))
        return self._iterate(0, z, z, Pm, 0.0, 1.0)

    def flux_dn(self, I_dn_level):
        """Downward diffuse flux from the m = 0 radiance at one level: 2 pi sum w mu I."""
        return 2 * np.pi * np.sum(self.w * self.mu * I_dn_level)

    def flux_up(self, I_up_level):
        return 2 * np.pi * np.sum(self.w * self.mu * I_up_level)


def radiance(modes, level, j, dphi_deg, direction="up"):
    """Sum Fourier modes at sublayer level `level`, direction index j, relative azimuth."""
    k = 0 if direction == "up" else 1
    return sum((1 if m == 0 else 2) * md[k][level, j] * np.cos(np.radians(m * dphi_deg))
               for m, md in enumerate(modes))


def surface_terms(tau_layers, ssa_layers, beta_layers, mu0, muv, dphi_deg, sensor_layer,
                  n_streams=8, max_dtau=0.01):
    """The four numbers the radiance equation needs, at one wavelength.

    Layers top to bottom; `sensor_layer` = index of the first layer below the sensor.
      R_path  path reflectance at the sensor (pi I_up / mu0), black surface
      T_dn    total (direct + diffuse) transmittance sun -> ground, per unit mu0
      T_up    total transmittance ground -> sensor for a Lambertian surface (by reciprocity:
              the slab below the sensor lit from above at muv)
      S       spherical albedo of the whole atmosphere (lit isotropically from below)
    so that   rho_toa_like = R_path + T_dn T_up rho / (1 - S rho).
    """
    tau_layers, ssa_layers, beta_layers = delta_m(np.asarray(tau_layers, float), np.asarray(ssa_layers, float),
                                                  np.asarray(beta_layers, float), 2 * n_streams)
    atm = Atmosphere.from_layers(tau_layers, ssa_layers, beta_layers, max_dtau, levels=(sensor_layer,))
    ks = atm.level_index[sensor_layer]
    s = SOS(atm, n_streams, extra_mu=[muv])
    jv = s.n_gauss                                   # index of the view direction
    modes = s.solar(mu0)
    R_path = np.pi * radiance(modes, ks, jv, dphi_deg, "up") / mu0
    T_dn = np.exp(-atm.tau[-1] / mu0) + s.flux_dn(modes[0][1][-1]) / mu0
    up, dn = s.isotropic_bottom()
    S = s.flux_dn(dn[-1]) / np.pi
    below = Atmosphere.from_layers(tau_layers[sensor_layer:], ssa_layers[sensor_layer:],
                                   beta_layers[sensor_layer:], max_dtau)
    sb = SOS(below, n_streams)
    mb = sb.solar(muv, n_modes=1)
    T_up = np.exp(-below.tau[-1] / muv) + sb.flux_dn(mb[0][1][-1]) / muv
    return dict(R_path=R_path, T_dn=T_dn, T_up=T_up, S=S)
