"""Physics checks for the SOS solver (no reference code: conservation laws and limits)."""
import numpy as np
from dustran import sos, rayleigh

RAY = rayleigh.phase_moments(550)


def fluxes(tau, ssa, beta, mu0, n=8, dt=0.01):
    a = sos.Atmosphere.from_layers([tau], [ssa], [beta], dt)
    s = sos.SOS(a, n)
    md = s.solar(mu0, n_modes=1)
    R = s.flux_up(md[0][0][0]) / mu0
    T = np.exp(-tau / mu0) + s.flux_dn(md[0][1][-1]) / mu0
    return R, T


def test_energy_conservation():
    for tau in (0.05, 0.3, 1.0):
        for mu0 in (1.0, 0.6):
            R, T = fluxes(tau, 1.0, RAY, mu0)
            assert abs(R + T - 1) < 2e-3, (tau, mu0, R, T)


def test_absorption_reduces_energy():
    R, T = fluxes(0.3, 0.9, RAY, 0.8)
    assert R + T < 0.97


def test_single_scattering_limit():
    """Very thin layer: reflectance -> ssa * tau * P(Theta) / (4 mu0 muv)."""
    tau, mu0, muv = 1e-4, 0.8, 0.9
    a = sos.Atmosphere.from_layers([tau], [1.0], [RAY], 1e-4)
    s = sos.SOS(a, 8, extra_mu=[muv])
    md = s.solar(mu0)
    for dphi in (0.0, 90.0, 180.0):
        I = sos.radiance(md, 0, s.n_gauss, dphi, "up")
        cos_t = -mu0 * muv + np.sqrt(1 - mu0 ** 2) * np.sqrt(1 - muv ** 2) * np.cos(np.radians(dphi))
        P = 1 + RAY[2] * 0.5 * (3 * cos_t ** 2 - 1)
        expect = tau * P / (4 * np.pi * muv)
        assert abs(I / expect - 1) < 2e-3, (dphi, I, expect)


def test_converges_with_streams_and_layers():
    a = fluxes(0.3, 1.0, RAY, 0.7, n=8, dt=0.01)[0]
    b = fluxes(0.3, 1.0, RAY, 0.7, n=16, dt=0.002)[0]
    assert abs(a / b - 1) < 5e-3


def test_spherical_albedo_thin_limit():
    """Thin conservative Rayleigh layer: S ~ tau (half of the scattered light goes back down,
    and the mean slant path of isotropic light is 2x the vertical)."""
    tau = 0.01
    a = sos.Atmosphere.from_layers([tau], [1.0], [RAY], 0.001)
    s = sos.SOS(a, 16)
    up, dn = s.isotropic_bottom()
    S = s.flux_dn(dn[-1]) / np.pi
    assert 0.8 * tau < S < 1.1 * tau, S


def test_mie_published_cases():
    from dustran.mie import mie
    x = 2 * np.pi * 0.525 / 0.6328                       # Bohren & Huffman appendix
    assert abs(mie(x, 1.55 + 0j, [1.0])[1] - 3.10543) < 1e-4
    assert abs(mie(100.0, 1.5 + 0j, [1.0])[0] - 2.0944) < 1e-4          # Wiscombe 1979
    qe, qs, _ = mie(10.0, 1.5 + 0.1j, [1.0])
    assert abs(qe - 2.459791) < 1e-5 and abs(qs - 1.235144) < 1e-5


def test_energy_conservation_forward_peaked():
    """Non-absorbing, strongly forward-peaked (Henyey-Greenstein g=0.7) with delta-M."""
    g, L = 0.7, 64
    beta = (2 * np.arange(L) + 1) * g ** np.arange(L)
    t, w, b = sos.delta_m(np.array([0.5]), np.array([1.0]), beta[None, :], 32)
    a = sos.Atmosphere.from_layers(t, w, b, 0.01)
    s = sos.SOS(a, 16)
    md = s.solar(0.8, n_modes=1)
    R = s.flux_up(md[0][0][0]) / 0.8
    T = np.exp(-a.tau[-1] / 0.8) + s.flux_dn(md[0][1][-1]) / 0.8
    assert abs(R + T - 1) < 3e-3, (R, T)
