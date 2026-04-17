#!/usr/bin/env python3
"""
evaporation_rate.py
-------------------
Compute the DM evaporation rate E_⊙(m_χ) from the Sun using Eq. (5.1)
from Garani & Palomares-Ruiz (2017) [arXiv:1702.02768].

Constant (velocity-independent, isotropic) cross section.
Targets: electrons (DM-e), hydrogen (DM-nucleon SD), all nuclei (DM-nucleon SI).
Produces a plot of E_⊙ vs m_χ analogous to Fig. 3 (top-left panel).
"""

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import erf, hyp0f1, gammaln, ive
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================
# Physical Constants (CGS)
# =============================================================
G_N   = 6.674e-8       # gravitational constant [cm³/(g s²)]
k_B   = 1.381e-16      # Boltzmann constant [erg/K]
m_u   = 1.661e-24      # atomic mass unit [g]
m_p_g = 1.672e-24      # proton mass [g]
m_e_g = 9.109e-28      # electron mass [g]
R_sun = 6.957e10       # solar radius [cm]
M_sun = 1.989e33       # solar mass [g]
GeV2g = 1.783e-24      # 1 GeV/c² → grams

# Elements for SI interactions: (name, Z, A, column_index_0based)
# Column indices in AGSS09: 6=H1, 7=He4, 8=He3, 9=C12, 10=C13,
# 11=N14, 12=N15, 13=O16, 14=O17, 15=O18, 16=Ne, 17=Na, 18=Mg,
# 19=Al, 20=Si, 21=P, 22=S, 23=Cl, 24=Ar, 25=K, 26=Ca, ..., 33=Fe
SI_ELEMENTS = [
    ("H1",  1,  1,   6),
    ("He4", 2,  4,   7),
    ("C12", 6,  12,  9),
    ("N14", 7,  14, 11),
    ("O16", 8,  16, 13),
    ("Ne",  10, 20, 16),
    ("Mg",  12, 24, 18),
    ("Si",  14, 28, 20),
    ("S",   16, 32, 22),
    ("Fe",  26, 56, 33),
]

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR  = os.path.dirname(BASE_DIR)
DATA_DIR  = os.path.join(ROOT_DIR, "data")
SSM_PATH  = os.path.join(DATA_DIR, "model_agss09.dat")
FIG_DIR   = os.path.join(ROOT_DIR, "figure")
os.makedirs(FIG_DIR, exist_ok=True)


# =============================================================
# Load AGSS09 Solar Model
# =============================================================
def load_solar_model():
    """Return raw AGSS09 data: array of shape (N, 35)."""
    rows = []
    with open(SSM_PATH) as f:
        for line in f:
            s = line.strip()
            if not s or not s[0].isdigit():
                continue
            vals = s.split()
            if len(vals) < 20:
                continue
            try:
                row = [float(x) for x in vals]
                rows.append(row)
            except ValueError:
                continue
    return np.array(rows)


def build_solar_profiles(data):
    """
    Compute derived solar quantities from AGSS09 data.
    Returns dict with keys:
      r, T, rho, n_e, n_H, phi, v_esc, X_frac (dict of element arrays)
    All in CGS.
    """
    r_frac  = data[:, 1]       # R/R_sun
    T_K     = data[:, 2]       # K
    rho     = data[:, 3]       # g/cm³
    M_frac  = data[:, 0]       # M/M_sun
    X_H     = data[:, 6]       # hydrogen mass fraction

    N = len(r_frac)
    r   = r_frac * R_sun                          # cm
    T   = k_B * T_K                               # erg
    M_r = M_frac * M_sun                          # enclosed mass [g]

    # Electron number density (fully ionized plasma)
    n_e = rho * (1.0 + X_H) / (2.0 * m_u)        # cm⁻³

    # Hydrogen number density
    n_H = rho * X_H / m_p_g                       # cm⁻³

    # Element mass-fraction arrays (for SI)
    X_frac = {}
    for name, Z, A, col in SI_ELEMENTS:
        if col < data.shape[1]:
            X_frac[name] = data[:, col]
        else:
            X_frac[name] = np.zeros(N)

    # Gravitational potential: φ(r) = ∫₀ʳ G M(r') / r'² dr'
    phi = np.zeros(N)
    for k in range(1, N):
        if r[k] <= 0:
            continue
        f_k   = G_N * M_r[k]   / r[k]**2
        f_km1 = G_N * M_r[k-1] / max(r[k-1], 1e-5)**2
        phi[k] = phi[k-1] + 0.5 * (f_k + f_km1) * (r[k] - r[k-1])

    # Escape velocity: v²_e(r) = 2[φ(R☉) - φ(r) + G M☉/R☉]
    phi_surf = phi[-1]
    v_esc_sq = 2.0 * (phi_surf - phi + G_N * M_sun / R_sun)
    v_esc    = np.sqrt(np.maximum(v_esc_sq, 0.0))

    return dict(r=r, T=T, rho=rho, n_e=n_e, n_H=n_H,
                phi=phi, v_esc=v_esc, X_frac=X_frac)


# =============================================================
# DM Temperature  (Knudsen limit, Eq. B.11, no cutoff)
# =============================================================
def solve_T_chi(m_chi_g, targets, sol):
    """
    Solve energy-balance equation (B.11) for DM temperature.

    targets : list of (m_target_g, n_target_array, weight)
              where weight = σ_{i,0} relative factor (cancels, set to 1
              unless multiple targets with different σ).
    sol     : solar profile dict

    Returns T_chi [erg].
    """
    r, T, phi = sol['r'], sol['T'], sol['phi']

    def balance(T_chi):
        total = 0.0
        boltz = np.exp(-m_chi_g * phi / T_chi)
        for m_i, n_i, w_i in targets:
            kin = m_chi_g * m_i / (m_i + m_chi_g)**2
            vel = np.sqrt((m_i * T_chi + m_chi_g * T) / (m_chi_g * m_i))
            integrand = w_i * n_i * kin * vel * (T - T_chi) * boltz * r**2
            total += trapezoid(integrand, r)
        return total

    T_lo = 0.3 * T.min()
    T_hi = 1.5 * T.max()
    try:
        return brentq(balance, T_lo, T_hi, xtol=1e-20, rtol=1e-8)
    except ValueError:
        # fallback: weighted average
        boltz = np.exp(-m_chi_g * phi / T[0])
        w = sol['n_e'] * boltz * r**2
        return trapezoid(w * T, r) / trapezoid(w, r)


# =============================================================
# R⁺  for constant cross section  (Eq. A.23)
# =============================================================
def chi_func(a, b):
    """χ(a,b) = √π/2 (erf(b) - erf(a))."""
    return 0.5 * np.sqrt(np.pi) * (erf(b) - erf(a))


def R_plus_const_vec(W, V, n_i, sigma_0, m_chi_g, m_i_g, T_sun):
    """
    Vectorised R⁺_const(w→v) for v>w, Eq. (A.23).

    W, V : arrays of DM velocities w and v [cm/s]  (broadcastable)
    Returns R⁺ array [cm⁻¹].
    """
    mu  = m_chi_g / m_i_g
    mup = (mu + 1.0) / 2.0
    mum = (mu - 1.0) / 2.0
    u_i = np.sqrt(2.0 * T_sun / m_i_g)

    alpha_m = (mup * V - mum * W) / u_i          # α₋
    alpha_p = (mup * V + mum * W) / u_i          # α₊
    beta_m  = (mum * V - mup * W) / u_i          # β₋
    beta_p  = (mum * V + mup * W) / u_i          # β₊

    c1 = chi_func(alpha_m, alpha_p)

    # Exponential factor (clip to avoid overflow)
    exponent = mu * (W**2 - V**2) / u_i**2
    exponent = np.clip(exponent, -500.0, 500.0)
    c2 = chi_func(beta_m, beta_p) * np.exp(exponent)

    prefactor = (2.0 / np.sqrt(np.pi)) * mup**2 / mu * (V / np.maximum(W, 1e-30)) * n_i * sigma_0
    R = prefactor * (c1 + c2)

    # Enforce v > w condition
    R = np.where(V > W, R, 0.0)
    return np.maximum(R, 0.0)


# =============================================================
# Column density for SD (hydrogen)
# =============================================================
def compute_column_density_H(sol):
    """
    Compute radial hydrogen column density from radius r to R_⊙:
    Σ_H(r) = ∫_r^{R_⊙} n_H(r') dr'   [cm⁻²]
    """
    r   = sol['r']
    n_H = sol['n_H']
    N   = len(r)
    col = np.zeros(N)
    for k in range(N - 2, -1, -1):
        col[k] = col[k + 1] + 0.5 * (n_H[k] + n_H[k + 1]) * (r[k + 1] - r[k])
    return col


# =============================================================
# Suppression factor  s(r)  —  Eq. (5.2) / (C.8)–(C.13)
# =============================================================
def compute_suppression(tau_arr, phi_hat_arr):
    """
    s(r) = η_ang(r) × η_mult(r) × exp(−τ(r))

    η_ang  : Eq. (C.11) — non-radial trajectories
    η_mult : Eq. (C.13) — multiple scatterings  ( ₀F₁(;1+2φ̂/3; τ) )

    Parameters
    ----------
    tau_arr     : array  – optical depth at each r
    phi_hat_arr : array  – φ̂(r) = m_χ v_e²(r) / (2 T⊙(r))

    Returns
    -------
    s : array,  0 ≤ s ≤ 1
    """
    tau = np.asarray(tau_arr, dtype=float)
    phi_hat = np.asarray(phi_hat_arr, dtype=float)
    s = np.zeros_like(tau)

    # ---- regime masks ----
    thin = tau < 0.01                          # optically thin  → s ≈ 1
    med  = (tau >= 0.01) & (tau < 500)         # moderate τ      → direct
    # thick (τ ≥ 500):  s < 10⁻²⁰⁰ → leave as 0

    s[thin] = 1.0

    if np.any(med):
        t = tau[med]
        p = phi_hat[med]
        b = 1.0 + 2.0 / 3.0 * p               # parameter of ₀F₁

        # η_ang (C.11)
        eta_ang = (7.0 / 10.0) * (1.0 - np.exp(-10.0 * t / 7.0)) / t

        # η_mult (C.13):  ₀F₁(; b; τ)
        eta_mult = np.ones_like(t)
        for j in range(len(t)):
            try:
                eta_mult[j] = float(hyp0f1(b[j], t[j]))
            except (OverflowError, FloatingPointError):
                # fallback: Bessel-ive representation
                nu  = b[j] - 1.0
                sq  = np.sqrt(t[j])
                log_val = (gammaln(b[j])
                           + (1.0 - b[j]) / 2.0 * np.log(t[j])
                           + np.log(max(float(ive(nu, 2.0 * sq)), 1e-300))
                           + 2.0 * sq - t[j])
                eta_mult[j] = np.exp(min(log_val, 700))

        vals = eta_ang * eta_mult * np.exp(-t)
        s[med] = np.where(np.isfinite(vals), vals, 0.0)

    return np.clip(s, 0.0, 1.0)


# =============================================================
# Thermal diffusivity α₀(μ)  — Gould & Raffelt 1990, Table I
# (constant, velocity-independent cross section)
# =============================================================
_alpha0_mu_table = np.array([
    0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0,
    5.0, 7.0, 10.0, 30.0, 100.0, 1000.0])
_alpha0_val_table = np.array([
    1.00, 1.003, 1.04, 1.38, 1.77, 2.22, 3.01, 4.37, 5.09, 5.51,
    5.45, 5.34, 5.24, 5.07, 5.01, 5.00])
_alpha0_interp = interp1d(np.log10(_alpha0_mu_table),
                          _alpha0_val_table,
                          kind='cubic', fill_value=(1.0, 5.0),
                          bounds_error=False)

def alpha0_const(mu):
    """Thermal diffusivity for constant cross section (Gould & Raffelt 1990)."""
    return float(_alpha0_interp(np.log10(max(mu, 1e-3))))


# =============================================================
# LTE distribution  n_χ,LTE(r)  — Eq. (4.7)
# =============================================================
def compute_n_LTE(m_chi_g, sol, alpha_val):
    """
    Compute n_χ,LTE(r) / N_χ   (normalised so that ∫ n 4πr² dr = 1).

    n_χ,LTE ∝ (T⊙/T⊙(0))^{3/2} exp(-∫₀ʳ [α dT/dr + m_χ dφ/dr] / T dr)

    Parameters
    ----------
    alpha_val : float – α₀(μ) for SD (single hydrogen target)
    """
    r   = sol['r']
    T   = sol['T']
    phi = sol['phi']
    N   = len(r)

    # Build the exponent via trapezoidal integration
    # integrand_k = [α dT/dr + m_χ dφ/dr] / T
    log_profile = np.zeros(N)
    for k in range(1, N):
        dT = T[k] - T[k - 1]
        dphi = phi[k] - phi[k - 1]
        dr = r[k] - r[k - 1]
        if dr <= 0 or T[k] <= 0:
            continue
        T_mid = 0.5 * (T[k] + T[k - 1])
        integrand = (alpha_val * dT + m_chi_g * dphi) / T_mid
        log_profile[k] = log_profile[k - 1] + integrand

    # n_LTE ∝ (T/T₀)^{3/2} exp(-log_profile)
    n_raw = (T / T[0])**1.5 * np.exp(-log_profile)
    n_raw = np.where(np.isfinite(n_raw), n_raw, 0.0)
    norm = trapezoid(n_raw * 4.0 * np.pi * r**2, r)
    if norm <= 0:
        return np.zeros(N)
    return n_raw / norm


# =============================================================
# Knudsen number & transition  — Eq. (4.10)–(4.12)
# =============================================================
def knudsen_number(sigma_0, n_H_center, rho_center, T_center, m_chi_g):
    """
    K = ℓ(0) / r_χ
    For SD (constant σ): ℓ(0) = 1/(σ₀ n_H(0))
    """
    ell_0 = 1.0 / (sigma_0 * n_H_center)
    r_chi = np.sqrt(3.0 * T_center /
                    (2.0 * np.pi * G_N * rho_center * m_chi_g))
    return ell_0 / r_chi


def f_K(K, K0=0.4):
    """Transition function f(K) = 1/(1 + (K/K₀)²), Eq. (4.12)."""
    return 1.0 / (1.0 + (K / K0)**2)


# =============================================================
# Full evaporation rate  E_⊙(σ)  with LTE/iso transition
# =============================================================
def compute_evap_vs_sigma_SD_full(m_chi_GeV, log_sigma_arr, sol,
                                  n_w=80, n_v=150):
    """
    Compute E_⊙(σ) for SD interactions with the full Knudsen↔LTE
    transition (Eq. 4.11) and suppression factor s(r).

    For each σ, the DM distribution is:
      n_χ f_χ = f(K) n_LTE f_LTE + (1-f(K)) n_iso f_iso
    and the evaporation integrand is recomputed accordingly.
    """
    m_chi_g = m_chi_GeV * GeV2g
    m_i_g   = m_p_g
    mu      = m_chi_g / m_i_g

    r   = sol['r']
    T   = sol['T']
    phi = sol['phi']
    v_e = sol['v_esc']
    n_H = sol['n_H']
    rho = sol['rho']
    N   = len(r)

    # ── α₀ for SD (single H target) ──
    alpha_val = alpha0_const(mu)

    # ── Knudsen-limit DM temperature ──
    target_list = [(m_p_g, n_H, 1.0)]
    T_chi = solve_T_chi(m_chi_g, target_list, sol)
    v_chi = np.sqrt(2.0 * T_chi / m_chi_g)

    # ── isothermal distribution g_iso(r) ──
    boltz_arg = m_chi_g * phi / T_chi
    boltz_arg = boltz_arg - boltz_arg.min()
    boltz_iso = np.exp(-boltz_arg)
    norm_iso  = trapezoid(boltz_iso * 4.0 * np.pi * r**2, r)
    g_iso     = boltz_iso / norm_iso

    # ── LTE distribution g_LTE(r) ──
    g_LTE = compute_n_LTE(m_chi_g, sol, alpha_val)

    # ── column density & φ̂ (σ-independent) ──
    col_H = compute_column_density_H(sol)
    phi_hat = m_chi_g * v_e**2 / (2.0 * T)
    phi_hat = np.where(T > 0, phi_hat, 1e30)

    # ── Gauss-Legendre for w integral ──
    gl_x, gl_w = np.polynomial.legendre.leggauss(n_w)

    # ── scale height & radial range ──
    r_chi = np.sqrt(3.0 * T_chi /
                    (2.0 * np.pi * G_N * rho[0] * m_chi_g))

    # ── precompute R⁺ integral at each shell (σ-independent part) ──
    # For iso: use f_χ with T_χ global
    # For LTE: use f_χ with T_⊙(r) local
    r_max = 0.95 * R_sun

    # We store the velocity-integral result per shell for both regimes
    Iwv_iso = np.zeros(N)   # ∫ f_iso(w) 4πw² dw ∫ R⁺(σ=1) dv
    Iwv_LTE = np.zeros(N)   # ∫ f_LTE(w) 4πw² dw ∫ R⁺(σ=1) dv

    # Masks for shells worth computing
    g_max_iso = np.max(g_iso * r**2)
    g_max_LTE = np.max(g_LTE * r**2) if np.max(g_LTE) > 0 else 1.0

    for k in range(N):
        if r[k] < 1e4 or r[k] > r_max:
            continue

        ve_k = v_e[k]
        T_k  = T[k]
        n_H_k = n_H[k]
        if ve_k < 1e3 or T_k < 1e-20 or n_H_k <= 0:
            continue

        vc_k = ve_k
        u_i  = np.sqrt(2.0 * T_k / m_i_g)
        mup  = (mu + 1.0) / 2.0

        v_max_k = max(4.0 * ve_k,
                      ve_k + 8.0 * u_i / np.sqrt(max(mu, 1.0)))
        v_pts = np.linspace(ve_k, v_max_k, n_v)

        w_pts = (gl_x + 1.0) / 2.0 * vc_k
        w_wts = gl_w * vc_k / 2.0

        # R⁺ grid (σ=1)
        WW, VV = np.meshgrid(w_pts, v_pts, indexing='ij')
        Rgrid = R_plus_const_vec(WW, VV, n_H_k, 1.0,
                                 m_chi_g, m_i_g, T_k)
        Iv = trapezoid(Rgrid, v_pts, axis=1)    # (n_w,)

        # --- isothermal f_χ(w) with global T_χ ---
        if g_iso[k] * r[k]**2 > 1e-30 * g_max_iso:
            xc_iso = vc_k / v_chi
            norm_f_iso = erf(xc_iso) - (2.0/np.sqrt(np.pi)) * xc_iso * np.exp(-xc_iso**2)
            if norm_f_iso > 1e-300:
                fw_iso = np.exp(-w_pts**2 / v_chi**2) / (np.pi**1.5 * v_chi**3 * norm_f_iso)
                Iwv_iso[k] = np.sum(fw_iso * 4.0*np.pi * w_pts**2 * Iv * w_wts)

        # --- LTE f_χ(w) with local T_⊙(r) ---
        if g_LTE[k] * r[k]**2 > 1e-30 * g_max_LTE:
            v_chi_local = np.sqrt(2.0 * T_k / m_chi_g)
            xc_lte = vc_k / v_chi_local
            norm_f_lte = erf(xc_lte) - (2.0/np.sqrt(np.pi)) * xc_lte * np.exp(-xc_lte**2)
            if norm_f_lte > 1e-300:
                fw_lte = np.exp(-w_pts**2 / v_chi_local**2) / (np.pi**1.5 * v_chi_local**3 * norm_f_lte)
                Iwv_LTE[k] = np.sum(fw_lte * 4.0*np.pi * w_pts**2 * Iv * w_wts)

    # ── sweep over σ ──
    E_arr = np.zeros(len(log_sigma_arr))
    for j, log_sig in enumerate(log_sigma_arr):
        sigma = 10.0**log_sig

        # Knudsen number
        K = knudsen_number(sigma, n_H[0], rho[0], T[0], m_chi_g)
        fK = f_K(K)

        # suppression factor
        tau = sigma * col_H
        s_r = compute_suppression(tau, phi_hat)

        # combined integrand: σ × s(r) × [fK g_LTE Iwv_LTE + (1-fK) g_iso Iwv_iso] × 4πr²
        integrand = sigma * s_r * (
            fK * g_LTE * Iwv_LTE + (1.0 - fK) * g_iso * Iwv_iso
        ) * 4.0 * np.pi * r**2

        E_arr[j] = trapezoid(integrand, r)

    return E_arr, T_chi


# =============================================================
# Capture Rate  C_⊙  (SD, constant σ)   — Eq. (3.1), (3.7), (3.9)
# =============================================================
# Halo parameters
rho_chi = 0.3 * GeV2g        # local DM density [g/cm³]  (0.3 GeV/cm³)
v_sun   = 220e5              # Sun speed w.r.t. DM rest frame [cm/s]
v_d     = 270e5              # velocity dispersion [cm/s]


def capture_rate_weak_SD(m_chi_GeV, sigma_0, sol, n_u=200):
    """
    Cweak_⊙ for SD (hydrogen) with constant cross section.
    C_weak = sigma_0 × _capture_rate_unit_SD(...)
    """
    return sigma_0 * _capture_rate_unit_SD(m_chi_GeV, sol, n_u)


def _capture_rate_unit_SD(m_chi_GeV, sol, n_u=200):
    """
    Compute C_weak / σ₀  (σ-independent part).
    Since R⁻ ∝ σ₀, we factor it out and compute once per mass.
    """
    m_chi_g = m_chi_GeV * GeV2g
    m_i_g   = m_p_g
    r   = sol['r']
    T   = sol['T']
    v_e = sol['v_esc']
    n_H = sol['n_H']
    N   = len(r)

    mu  = m_chi_g / m_i_g
    mup = (mu + 1.0) / 2.0
    mum = (mu - 1.0) / 2.0

    # DM halo velocity grid
    u_max = v_sun + 4.0 * v_d
    u_pts = np.linspace(1e3, u_max, n_u)

    # f_{v⊙}(u) — Eq. (3.6)
    coeff = np.sqrt(3.0 / (2.0 * np.pi))
    f_halo = coeff * u_pts / (v_sun * v_d) * (
        np.exp(-3.0 * (u_pts - v_sun)**2 / (2.0 * v_d**2))
      - np.exp(-3.0 * (u_pts + v_sun)**2 / (2.0 * v_d**2))
    )

    # Gauss-Legendre for v integral (down-scatter)
    gl_x_v, gl_w_v = np.polynomial.legendre.leggauss(60)

    # Radial + halo velocity double integral  (vectorized over u)
    C_total = 0.0
    for k in range(N - 1):
        if r[k] < 1e4 or n_H[k] <= 0:
            continue
        ve_k  = v_e[k]
        T_k   = T[k]
        n_H_k = n_H[k]
        u_i   = np.sqrt(2.0 * T_k / m_i_g)

        # v grid for R⁻: [0, ve_k]  (n_gl,)
        v_pts_k = (gl_x_v + 1.0) / 2.0 * ve_k
        v_wts_k = gl_w_v * ve_k / 2.0

        # w(r, u) for all halo velocities  (n_u,)
        w_all = np.sqrt(u_pts**2 + ve_k**2)

        # Broadcast: V (n_gl, 1),  W (1, n_u)
        V = v_pts_k[:, None]
        W = w_all[None, :]

        alpha_m = (mup * V - mum * W) / u_i
        alpha_p = (mup * V + mum * W) / u_i
        beta_m  = (mum * V - mup * W) / u_i
        beta_p  = (mum * V + mup * W) / u_i

        c1 = chi_func(-alpha_m, alpha_p)
        exponent = np.clip(mu * (W**2 - V**2) / u_i**2, -500, 500)
        c2 = chi_func(-beta_m, beta_p) * np.exp(exponent)

        # σ=1 here (factored out)
        prefactor = (2.0/np.sqrt(np.pi)) * mup**2/mu * (V/W) * n_H_k
        R_minus = np.maximum(prefactor * (c1 + c2), 0.0)   # (n_gl, n_u)
        R_minus = np.where(V < W, R_minus, 0.0)

        # ∫ R⁻ dv for each u  → (n_u,)
        Iv = np.sum(R_minus * v_wts_k[:, None], axis=0)

        inner_u = (rho_chi / m_chi_g) * f_halo / u_pts * w_all * Iv

        shell_contrib = trapezoid(inner_u, u_pts)
        dr = r[k+1] - r[k]
        C_total += shell_contrib * 4.0 * np.pi * r[k]**2 * dr

    return C_total


def capture_rate_geom(m_chi_GeV, sol):
    """
    Geometric capture rate limit — Eq. (3.7).
    """
    m_chi_g = m_chi_GeV * GeV2g
    ve_surf = sol['v_esc'][-1]

    v_avg_0 = np.sqrt(8.0 / (3.0 * np.pi)) * v_d

    # ξ(v⊙, v_d) — Eq. (3.8)
    x = np.sqrt(1.5) * v_sun / v_d
    xi = (v_d**2 * np.exp(-1.5 * v_sun**2 / v_d**2)
          + np.sqrt(np.pi / 6.0) * (v_d / v_sun)
          * (v_d**2 + 3.0 * ve_surf**2 + 3.0 * v_sun**2) * erf(x)) \
         / (2.0 * v_d**2 + 3.0 * ve_surf**2)

    C_geom = np.pi * R_sun**2 * (rho_chi / m_chi_g) * v_avg_0 \
             * (1.0 + 1.5 * ve_surf**2 / v_d**2) * xi
    return C_geom


def capture_rate_SD(m_chi_GeV, sigma_0, sol):
    """
    Full capture rate with saturation — Eq. (3.9).
    C = Cweak × (1 - exp(-Cgeom/Cweak))
    """
    C_weak = capture_rate_weak_SD(m_chi_GeV, sigma_0, sol)
    C_geom = capture_rate_geom(m_chi_GeV, sol)
    if C_weak > 0:
        ratio = C_geom / C_weak
        return C_weak * (1.0 - np.exp(-ratio))
    return 0.0


# =============================================================
# Annihilation Rate  A_⊙  — Eq. (4.15)
# =============================================================
sigma_A_v = 3e-26            # ⟨σ_A v⟩ [cm³/s]  (thermal relic)


def annihilation_rate(m_chi_GeV, sigma_0, sol):
    """
    A_⊙ = ⟨σ_A v⟩ × ∫ n²_χ 4πr² dr / (∫ n_χ 4πr² dr)²
    Uses the Knudsen↔LTE transition distribution.
    """
    m_chi_g = m_chi_GeV * GeV2g
    r   = sol['r']
    T   = sol['T']
    phi = sol['phi']
    n_H = sol['n_H']
    rho = sol['rho']

    mu = m_chi_g / m_p_g
    alpha_val = alpha0_const(mu)

    # isothermal distribution
    target_list = [(m_p_g, n_H, 1.0)]
    T_chi = solve_T_chi(m_chi_g, target_list, sol)

    boltz_arg = m_chi_g * phi / T_chi
    boltz_arg = boltz_arg - boltz_arg.min()
    g_iso = np.exp(-boltz_arg)
    norm_iso = trapezoid(g_iso * 4.0 * np.pi * r**2, r)
    g_iso /= norm_iso

    # LTE distribution
    g_LTE = compute_n_LTE(m_chi_g, sol, alpha_val)

    # Knudsen transition
    K = knudsen_number(sigma_0, n_H[0], rho[0], T[0], m_chi_g)
    fK = f_K(K)

    g_total = fK * g_LTE + (1.0 - fK) * g_iso

    int_n2 = trapezoid(g_total**2 * 4.0 * np.pi * r**2, r)
    int_n  = trapezoid(g_total * 4.0 * np.pi * r**2, r)

    return sigma_A_v * int_n2 / int_n**2


# =============================================================
# N_χ(t_⊙)  — Eq. (6.1)–(6.2)
# =============================================================
t_sun = 4.57e9 * 3.156e7     # 4.57 Gyr in seconds


def N_total(C_sun, A_sun, E_sun):
    """
    Eq. (6.2): N_χ(t⊙) = √(C/A) × tanh(κ t⊙/τ_eq) / (κ + E τ_eq tanh/2)
    """
    if C_sun <= 0 or A_sun <= 0:
        return 0.0
    tau_eq = 1.0 / np.sqrt(A_sun * C_sun)
    kappa  = np.sqrt(1.0 + (E_sun * tau_eq / 2.0)**2)
    arg    = kappa * t_sun / tau_eq
    th     = np.tanh(min(arg, 500.0))
    N = np.sqrt(C_sun / A_sun) * th / (kappa + 0.5 * E_sun * tau_eq * th)
    return max(N, 0.0)


# =============================================================
# Batch sweep:  C_⊙(σ), A_⊙(σ), N_χ(σ) for one mass
# =============================================================
def compute_N_vs_sigma_SD(m_chi_GeV, log_sigma_arr, E_arr, sol):
    """
    Compute N_χ(t_⊙) for an array of σ values, given pre-computed E_⊙(σ).

    Key optimisation: C_weak ∝ σ  →  compute σ-independent unit once.
    g_iso, g_LTE are σ-independent → precompute once, vary f(K) only.

    Returns
    -------
    N_arr, N_no_evap, C_arr, A_arr : arrays of length len(log_sigma_arr)
    """
    m_chi_g = m_chi_GeV * GeV2g
    r   = sol['r']
    T   = sol['T']
    phi = sol['phi']
    n_H = sol['n_H']
    rho = sol['rho']

    mu = m_chi_g / m_p_g

    # ── Capture rate: σ-independent unit (computed once) ──
    C_unit = _capture_rate_unit_SD(m_chi_GeV, sol)
    C_geom = capture_rate_geom(m_chi_GeV, sol)

    # ── Annihilation: precompute distributions (σ-independent) ──
    alpha_val = alpha0_const(mu)
    target_list = [(m_p_g, n_H, 1.0)]
    T_chi = solve_T_chi(m_chi_g, target_list, sol)

    boltz_arg = m_chi_g * phi / T_chi
    boltz_arg = boltz_arg - boltz_arg.min()
    g_iso = np.exp(-boltz_arg)
    norm_iso = trapezoid(g_iso * 4.0 * np.pi * r**2, r)
    g_iso /= norm_iso

    g_LTE = compute_n_LTE(m_chi_g, sol, alpha_val)

    # Precompute integrands for annihilation
    r2_4pi = 4.0 * np.pi * r**2
    iso2_integrand = g_iso**2 * r2_4pi
    lte2_integrand = g_LTE**2 * r2_4pi
    iso_lte_integrand = g_iso * g_LTE * r2_4pi   # cross term for n²
    iso_int  = trapezoid(g_iso * r2_4pi, r)      # should be 1.0
    lte_int  = trapezoid(g_LTE * r2_4pi, r)      # should be 1.0
    iso2_int = trapezoid(iso2_integrand, r)
    lte2_int = trapezoid(lte2_integrand, r)
    iso_lte_int = trapezoid(iso_lte_integrand, r)

    n_sigma = len(log_sigma_arr)
    N_arr      = np.zeros(n_sigma)
    N_no_evap  = np.zeros(n_sigma)
    C_arr      = np.zeros(n_sigma)
    A_arr      = np.zeros(n_sigma)

    for j, ls in enumerate(log_sigma_arr):
        sigma = 10.0**ls

        # Capture rate
        C_weak = sigma * C_unit
        if C_weak > 0:
            C_val = C_weak * (1.0 - np.exp(-C_geom / C_weak))
        else:
            C_val = 0.0

        # Annihilation rate with f(K) mixing
        K = knudsen_number(sigma, n_H[0], rho[0], T[0], m_chi_g)
        fK_val = f_K(K)
        # ∫ n² dV = fK² ∫g_LTE² + 2fK(1-fK) ∫g_iso g_LTE + (1-fK)² ∫g_iso²
        int_n2 = (fK_val**2 * lte2_int
                  + 2.0 * fK_val * (1.0 - fK_val) * iso_lte_int
                  + (1.0 - fK_val)**2 * iso2_int)
        # ∫ n dV = fK ∫g_LTE + (1-fK) ∫g_iso = fK + (1-fK) = 1
        int_n = fK_val * lte_int + (1.0 - fK_val) * iso_int
        A_val = sigma_A_v * int_n2 / int_n**2

        C_arr[j] = C_val
        A_arr[j] = A_val
        N_arr[j] = N_total(C_val, A_val, E_arr[j])
        N_no_evap[j] = N_total(C_val, A_val, 0.0)

    return N_arr, N_no_evap, C_arr, A_arr
