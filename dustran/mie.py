"""Mie scattering by homogeneous spheres, and size-distribution averages.

Written for DUSTRAN from the textbook series (Bohren & Huffman 1983, ch. 4): Riccati-Bessel
functions by upward recurrence, the logarithmic derivative D_n(mx) by downward recurrence
(stable for absorbing particles), and the amplitude functions S1, S2 from the angular functions
pi_n, tau_n.

Dust is not spherical; spheres get extinction and absorption about right but misplace side- and
back-scatter by tens of percent. A non-spherical (spheroid) correction is a known upgrade.
"""
import numpy as np


def mie(x, m, mu):
    """Single sphere: size parameter x = 2 pi r / lambda, complex index m = n + ik, cosines mu.
    Returns Qext, Qsca, and |S1|^2 + |S2|^2 at the given scattering-angle cosines."""
    nstop = int(x + 4 * x ** (1 / 3) + 2)
    nmx = int(max(nstop, abs(m * x)) + 16)
    y = m * x
    D = np.zeros(nmx + 1, complex)
    for n in range(nmx, 0, -1):
        D[n - 1] = n / y - 1 / (D[n] + n / y)
    n = np.arange(1, nstop + 1)
    psi = np.zeros(nstop + 1); chi = np.zeros(nstop + 1)
    psi[0], chi[0] = np.sin(x), np.cos(x)
    psim1, chim1 = np.cos(x), -np.sin(x)
    p0, c0, pm, cm = psi[0], chi[0], psim1, chim1
    for k in range(1, nstop + 1):
        p1 = (2 * k - 1) / x * p0 - pm
        c1 = (2 * k - 1) / x * c0 - cm
        psi[k], chi[k] = p1, c1
        pm, cm, p0, c0 = p0, c0, p1, c1
    xi = psi - 1j * chi
    Dn = D[1:nstop + 1]
    a = ((Dn / m + n / x) * psi[1:] - psi[:-1]) / ((Dn / m + n / x) * xi[1:] - xi[:-1])
    b = ((m * Dn + n / x) * psi[1:] - psi[:-1]) / ((m * Dn + n / x) * xi[1:] - xi[:-1])
    f = 2 * n + 1
    qext = 2 / x ** 2 * np.sum(f * (a.real + b.real))
    qsca = 2 / x ** 2 * np.sum(f * (np.abs(a) ** 2 + np.abs(b) ** 2))
    # angular functions
    mu = np.asarray(mu, float)
    pi_prev = np.zeros_like(mu); pi_n = np.ones_like(mu)
    S1 = np.zeros(len(mu), complex); S2 = np.zeros(len(mu), complex)
    for k in range(1, nstop + 1):
        tau_n = k * mu * pi_n - (k + 1) * pi_prev
        c = (2 * k + 1) / (k * (k + 1))
        S1 += c * (a[k - 1] * pi_n + b[k - 1] * tau_n)
        S2 += c * (a[k - 1] * tau_n + b[k - 1] * pi_n)
        pi_prev, pi_n = pi_n, ((2 * k + 1) * mu * pi_n - (k + 1) * pi_prev) / k
    return qext, qsca, np.abs(S1) ** 2 + np.abs(S2) ** 2


def lognormal_optics(wl_um, m, r_g, sigma_g, n_moments=64, n_radii=120, n_angles=721):
    """Optical properties of a lognormal number distribution (median radius r_g [um],
    geometric std sigma_g) at one wavelength. See distribution_optics for the output."""
    lnr = np.linspace(np.log(r_g) - 4 * np.log(sigma_g), np.log(r_g) + 4 * np.log(sigma_g), n_radii)
    w = np.exp(-0.5 * ((lnr - np.log(r_g)) / np.log(sigma_g)) ** 2)
    return distribution_optics(wl_um, m, np.exp(lnr), w / w.sum(), n_moments, n_angles)


def distribution_optics(wl_um, m, r, w, n_moments=64, n_angles=721):
    """Optical properties of particles with radii r [um] and number weights w (summing to 1).

    Returns dict: ext [um^2 per particle], ssa, g, beta (Legendre moments, beta_0 = 1,
    phase = sum beta_l P_l(cos)), and the phase function on its angle grid.
    """
    mu_g, mu_w = np.polynomial.legendre.leggauss(n_angles)
    k = 2 * np.pi / wl_um
    Cext = Csca = 0.0
    P = np.zeros(n_angles)
    for ri, wi in zip(r, w):
        qe, qs, s = mie(k * ri, m, mu_g)
        Cext += wi * qe * np.pi * ri ** 2
        Csca += wi * qs * np.pi * ri ** 2
        P += wi * s / k ** 2                  # differential scattering cross section x 2
    P = P / (0.5 * np.sum(mu_w * P))         # normalise: (1/2) int P dmu = 1
    beta = np.array([(2 * l + 1) / 2 * np.sum(mu_w * P * np.polynomial.legendre.Legendre.basis(l)(mu_g))
                     for l in range(n_moments)])
    return dict(ext=Cext, ssa=Csca / Cext, g=beta[1] / 3, beta=beta, mu=mu_g, phase=P)
