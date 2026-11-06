#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diseño de dique en talud — Manto principal de bloques de hormigón o escollera
Autores: R. Bergillos, M. Ortega, M. Clavero  |  GDFA · IISTA · UGR · 2018-2025
Versión web (Streamlit): M. Clavero Gilabert — 2026
"""

from __future__ import division

import math

import numpy as np
import streamlit as st

# ═══════════════════════════════════════════════════════════════════
#  FUNCIONES DE CÁLCULO  (idénticas al GUI de escritorio)
# ═══════════════════════════════════════════════════════════════════


def wavnum(t, h):
    """Número de onda — resuelve σ²/g = k·tanh(k·h) con fsolve."""
    from scipy.optimize import fsolve

    g = 9.81
    gamma = (2 * np.pi / t) ** 2.0 * h / g
    sigma2 = (2 * np.pi / t) ** 2

    def func(k):
        return sigma2 / g - k * np.tanh(k * h)

    ko = np.sqrt(gamma / h) if gamma > np.pi**2 else gamma / h
    k = fsolve(func, ko)
    return float(abs(k[0]))


def calcular_ir(params):
    """Calcula Ir a partir de las condiciones de oleaje (sin necesitar φ).
    Propaga desde SIMAR si oleaje_origen == 'simar'; usa Hs directa si es 'pie_dique'."""
    T = params["T"]
    h = params["h"]
    Hs0 = params["Hs0"]
    tanb = params["tanb"]
    theta0 = params.get("theta0", 0.0)
    oleaje_origen = params.get("oleaje_origen", "simar")
    k1 = wavnum(T, h)
    L1 = 2 * math.pi / k1
    Lo = 9.81 * T**2 / (2 * math.pi)
    if oleaje_origen == "simar":
        sigma = 2 * math.pi / T
        L0 = 9.81 * T**2 / (2 * math.pi)
        k0 = 2 * math.pi / L0
        cg0 = 0.5 * L0 / T
        cg1 = 0.5 * sigma / k1 * (1 + 2 * k1 * h / math.sinh(2 * k1 * h))
        theta1 = math.asin(k0 * math.sin(math.radians(theta0)) / k1)
        Kr = math.sqrt(math.cos(math.radians(theta0)) / math.cos(theta1))
        Ks = math.sqrt(cg0 / cg1)
        Hs = Hs0 * Ks * Kr
    else:
        Hs = Hs0
    R = 0.0
    Hmax = 1.6 * Hs
    Hrot = (0.11 + 0.03 * (1 - R) / (1 + R)) * math.tanh(k1 * h) * L1
    Hp = min(Hmax, Hrot)
    Ir = tanb / math.sqrt(Hp / Lo)
    return Ir, Hp, L1


def calcular(params):
    """Cálculo completo del dique. Devuelve dict de resultados o lanza ValueError."""
    gammaw = params["gammaw"]
    gammas = params["gammas"]
    gammae = params["gammae"]
    tanb = params["tanb"]
    Hs0 = params["Hs0"]
    T = params["T"]
    theta0 = params["theta0"]
    h = params["h"]
    phi = params["phi"]
    aux_fc = params["aux_fc"]
    qmax = params["qmax"]
    prob = params["prob"]
    con_berma = params["con_berma"]
    hb = params["hb"]
    Fm = params["Fm"]
    CS = params["CS"]
    mu = params["mu"]
    botaolas = int(params["botaolas"])
    Ae = params["Ae"]

    Sr = gammas / gammaw
    Sr_e = gammae / gammaw
    Delta_e = Sr_e - 1
    alpha = math.atan(tanb)

    # ── Propagación ──
    oleaje_origen = params.get("oleaje_origen", "simar")
    sigma = 2 * math.pi / T
    L0 = 9.81 * T**2 / (2 * math.pi)
    k0 = 2 * math.pi / L0
    k1 = wavnum(T, h)
    L1 = 2 * math.pi / k1
    cg0 = 0.5 * L0 / T
    cg1 = 0.5 * sigma / k1 * (1 + 2 * k1 * h / math.sinh(2 * k1 * h))
    if oleaje_origen == "simar":
        theta1 = math.asin(k0 * math.sin(math.radians(theta0)) / k1)
        Kr = math.sqrt(math.cos(math.radians(theta0)) / math.cos(theta1))
        Ks = math.sqrt(cg0 / cg1)
        Hs = Hs0 * Ks * Kr
    else:
        Ks, Kr = 1.0, 1.0
        Hs = Hs0
    L = L1
    k = k1

    # ── Rotura ──
    R = 0.0
    Hmax = 1.6 * Hs
    Hrot = (0.11 + 0.03 * (1 - R) / (1 + R)) * math.tanh(k * h) * L
    Hp = min(Hmax, Hrot)
    rotura = Hp == Hrot

    # ── Iribarren ──
    Lo = 9.81 * T**2 / 2 / math.pi
    Ir = tanb / math.sqrt(Hp / Lo)
    Ir_s = tanb / math.sqrt(Hs / Lo)
    Ad, Bd = -0.72, -0.42
    Rd = Ad * Hp * (1 - math.exp(Bd * Ir))

    # ── Manto principal ──
    tipo_pieza = params.get("tipo_pieza", "hormigon")
    if tipo_pieza == "escollera":
        Sr_manto = Sr_e
        gamma_manto = gammae
    else:
        Sr_manto = Sr
        gamma_manto = gammas
    Wp = phi * gammaw * Hp**3 * Sr_manto / (Sr_manto - 1) ** 3
    ep = 2 * (Wp / gamma_manto) ** (1 / 3)

    # ── Mantos secundarios ──
    Ws = np.zeros(4)
    Ws[0] = Wp / 15
    for i in range(3):
        capa = Ws[i] / 15
        if capa >= 0.150:
            Ws[i + 1] = capa
    es = np.zeros(4)
    for i in range(4):
        es[i] = 2 * (Ws[i] / gammae) ** (1 / 3)
    if Ws[0] / 15 <= 0.150:
        Ws[1] = 2.5 * 0.150
        es[1] = (Ws[1] / gammas) ** (1 / 3)

    # ── Berma de coronación ──
    Bc = ep
    Ac1 = 3 * Hs / 4
    Ac2 = Ae + np.sum(es) + ep
    Ac = max(Ac1, Ac2)

    # ── Berma de pie ──
    Wb, Db = None, None
    if con_berma:
        if hb <= 0:
            raise ValueError("hb debe ser > 0 para diseñar la berma de pie.")
        Nod = 0.5
        Ns = (2.4 * Nod ** (1 / 3)) / (1.4 - 0.4 * (hb / Hs))
        Db = Hs / (Delta_e * Ns)
        Wb = gammae * Db**3

    # ── Francobordo ──
    Au, Bu = 1.05, -0.67
    if aux_fc == 1:
        Hrms = Hs / math.sqrt(2)
        Hc = Hrms * (-math.log(prob)) ** 0.5
        Irc = tanb / math.sqrt(Hc / Lo)
        Ruf = Au * Hc * (1 - math.exp(Bu * Irc))
        Fc = math.ceil(Ruf * 10) / 10
    elif aux_fc == 2:
        Ru = Au * Hs * (1 - math.exp(Bu * Ir_s))
        Fc = math.ceil(Ru * 10) / 10
    else:
        aq, bq, cq = 0.09, 1.5 / 0.47, 1.3
        Fc_ = (Hs / bq) * (-math.log(qmax / (aq * math.sqrt(9.81 * Hs**3)))) ** (1 / cq)
        Fc = math.ceil(Fc_ * 10) / 10

    # ── Comprobación geométrica ──
    Rumax = Au * Hp * (1 - math.exp(Bu * Ir))
    So = Hp * (1 - Ac / Rumax)

    adv_espaldon = None
    if So <= 0:
        adv_espaldon = (
            f"So = {So:.3f} m (negativo): el run-up no alcanza la cota de coronación "
            f"(Ac = {Ac:.3f} m > Rumax = {Rumax:.3f} m). "
            "No se generan fuerzas dinámicas sobre el espaldón con estos parámetros."
        )

    # ── Espaldón (solo si So > 0) ──
    B_p_desliz = B_p_vuelco = B_r_desliz = B_r_vuelco = None
    coef_mu = Bmax = None

    if So > 0:
        lamb = 0.8 * math.exp(-10.9 * Bc / L)
        cf = 1.45
        alphao = 2 * cf * ((Rumax / Hp) * math.cos(alpha)) ** 2
        pso = alphao * gammaw * So
        pu = lamb * pso

        coef_a = 0.362
        coef_b = 0.069
        coef_c = 357.1
        coef_mu = coef_a * math.exp(coef_c * ((Hp / L - coef_b) ** 2))

        FH_p = (Fc - Ac) * pso + (Ac - Ae) * pu
        B1_p_desliz = FH_p * CS / mu
        B2_p_desliz = ((Fc - Fm) / botaolas) * gammas + (Fm - Ae) * gammas - 0.5 * pu
        if abs(B2_p_desliz) < 1e-10:
            raise ValueError("Denominador nulo en deslizamiento dinámico — revisa Fm, Ae y Fc.")
        B_p_desliz = B1_p_desliz / B2_p_desliz

        Mh_p = (Fc - Ac) * pso * (0.5 * (Fc - Ac) + (Ac - Ae)) + (Ac - Ae) * pu * (Ac - Ae) / 2
        B1_p_vuelco = Mh_p * CS
        B2_p_vuelco = (
            ((Fc - Fm) / botaolas) * gammas * (1 - 1 / (2 * botaolas))
            + (Fm - Ae) * gammas * 0.5
            - CS * pu / 3
        )
        ratio_p = B1_p_vuelco / B2_p_vuelco if abs(B2_p_vuelco) > 1e-10 else None
        if ratio_p is None or ratio_p < 0:
            raise ValueError(
                f"No se puede calcular vuelco dinámico (B1/B2 = {ratio_p}). "
                "Comprueba los valores de Fm, CS y los parámetros del espaldón."
            )
        B_p_vuelco = math.sqrt(ratio_p)
        Bmax_p = max(B_p_desliz, B_p_vuelco)

        Pr_abajo = coef_mu * gammaw * (So + Ac - Ae)
        Pr_arriba = coef_mu * gammaw * (So + Ac - Fc)
        FH_r = (Fc - Ae) * Pr_arriba + 0.5 * (Pr_abajo - Pr_arriba) * (Fc - Ae)
        B1_r_desliz = FH_r * CS / mu
        B2_r_desliz = (
            ((Fc - Fm) / botaolas) * gammas + (Fm - Ae) * gammas - 0.5 * Pr_abajo
        )
        if abs(B2_r_desliz) < 1e-10:
            raise ValueError("Denominador nulo en deslizamiento reflectivo.")
        B_r_desliz = B1_r_desliz / B2_r_desliz

        B1_r_vuelco = Mh_p * CS
        B2_r_vuelco = (
            ((Fc - Fm) / botaolas) * gammas * (1 - 1 / (2 * botaolas))
            + (Fm - Ae) * gammas * 0.5
            - 0.5 * Pr_abajo * 2 / 3
        )
        ratio_r = B1_r_vuelco / B2_r_vuelco if abs(B2_r_vuelco) > 1e-10 else None
        if ratio_r is None or ratio_r < 0:
            raise ValueError(
                f"No se puede calcular vuelco reflectivo (B1/B2 = {ratio_r}). "
                "Comprueba los parámetros del espaldón."
            )
        B_r_vuelco = math.sqrt(ratio_r)
        Bmax_r = max(B_r_desliz, B_r_vuelco)
        Bmax = max(Bmax_p, Bmax_r)

    return {
        "Hs": Hs, "L": L, "Hp": Hp, "rotura": rotura, "Ir": Ir,
        "Rd": Rd, "Rumax": Rumax, "So": So, "Wp": Wp, "ep": ep,
        "Ws": Ws, "es": es, "Bc": Bc, "Ac": Ac, "Fc": Fc,
        "B_p_desliz": B_p_desliz, "B_p_vuelco": B_p_vuelco,
        "B_r_desliz": B_r_desliz, "B_r_vuelco": B_r_vuelco,
        "coef_mu": coef_mu, "Bmax": Bmax, "Wb": Wb, "Db": Db,
        "adv_espaldon": adv_espaldon,
        "Ks": Ks if oleaje_origen == "simar" else 1.0,
        "Kr": Kr if oleaje_origen == "simar" else 1.0,
    }


# ═══════════════════════════════════════════════════════════════════
#  CONFIGURACIÓN DE PÁGINA
# ═══════════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="Dique en Talud — GDFA·UGR",
    page_icon="🌊",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background-color: #0d2b55;
}
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] p,
[data-testid="stSidebar"] span,
[data-testid="stSidebar"] div,
[data-testid="stSidebar"] .stMarkdown {
    color: #e8edf4 !important;
}
[data-testid="stSidebar"] .stExpander {
    background-color: #142f5c !important;
    border: 1px solid #1a4a8a !important;
    border-radius: 8px;
    margin-bottom: 6px;
}
[data-testid="stSidebar"] .stExpander summary {
    color: #ffffff !important;
    font-weight: 600;
}
/* Estado abierto: mismo tono que el hover nativo de Streamlit */
[data-testid="stSidebar"] details[open] > summary {
    border-radius: 6px 6px 0 0;
}
[data-testid="stSidebar"] details[open],
[data-testid="stSidebar"] details[open] > summary,
[data-testid="stSidebar"] details[open] > div,
[data-testid="stSidebar"] details[open] .streamlit-expanderContent,
[data-testid="stSidebar"] details[open] .st-expander-content {
    background-color: #162d58 !important;
}
[data-testid="stSidebar"] details[open] > div * {
    background-color: transparent !important;
}
[data-testid="stSidebar"] details[open] input {
    background-color: #122648 !important;
}
[data-testid="stSidebar"] input,
[data-testid="stSidebar"] select {
    background-color: #1a3a6b !important;
    color: #ffffff !important;
    border: 1px solid #3a6aab !important;
}
/* ── Cabecera principal ── */
.app-header {
    background: linear-gradient(135deg, #0d2b55 0%, #1a4a8a 100%);
    color: white;
    padding: 1.5rem 2rem;
    border-radius: 12px;
    margin-bottom: 1.5rem;
}
.app-header h1 { color: white !important; margin: 0; font-size: 1.6rem; }
.app-header p  { color: #d0d8e4 !important; margin: 0.3rem 0 0 0; font-size: 0.9rem; }
/* ── Tarjetas de resultados ── */
.metric-card {
    background: white;
    border: 1px solid #d0d8e4;
    border-radius: 10px;
    padding: 1rem;
    text-align: center;
    box-shadow: 0 2px 6px rgba(13,43,85,.08);
}
/* ── Tabla de resultados ── */
table { width: 100%; border-collapse: collapse; }
th { background: #0d2b55; color: white; padding: 8px 12px; text-align: left; }
td { padding: 6px 12px; border-bottom: 1px solid #e8edf4; }
tr:hover td { background: #f2f4f7; }
/* ── Tabla resumen — título expander y celdas ── */
[data-testid="stMain"] [data-testid="stExpander"] summary p,
[data-testid="stMain"] [data-testid="stExpander"] summary span {
    font-size: 20px !important;
    font-weight: 600 !important;
}
/* ── Botón principal ── */
div[data-testid="stButton"] > button {
    background-color: #0d2b55;
    color: white;
    border: none;
    border-radius: 8px;
    padding: 0.6rem 2.5rem;
    font-size: 1rem;
    font-weight: 600;
    letter-spacing: .05em;
    width: 100%;
    transition: background .2s;
}
div[data-testid="stButton"] > button:hover {
    background-color: #1a4a8a;
}
/* ── Ocultar botones +/- de los campos numéricos ── */
button[data-testid="stNumberInputStepUp"],
button[data-testid="stNumberInputStepDown"] {
    display: none;
}
</style>
""", unsafe_allow_html=True)

# ═══════════════════════════════════════════════════════════════════
#  CABECERA
# ═══════════════════════════════════════════════════════════════════

st.markdown("""
<div class="app-header">
  <h1>🌊 DIMENSIONAMIENTO DE DIQUE EN TALUD</h1>
  <p>M. Clavero &nbsp;·&nbsp; GDFA · IISTA · UGR</p>
</div>
""", unsafe_allow_html=True)

# ═══════════════════════════════════════════════════════════════════
#  BARRA LATERAL — DATOS DE ENTRADA
# ═══════════════════════════════════════════════════════════════════

with st.sidebar:
    st.markdown("## 📋 DATOS DE ENTRADA")

    # ── 1. Geometría y materiales ──────────────────────────────────
    with st.expander("Geometría y materiales", expanded=False):
        gammaw = st.number_input("γ agua (t/m³)", value=1.025, step=0.001, format="%.3f", key="gammaw")
        gammas = st.number_input("γ hormigón (t/m³)", value=2.35, step=0.01, format="%.3f", key="gammas")
        gammae = st.number_input("γ piedra escollera (t/m³)", value=2.65, step=0.01, format="%.3f", key="gammae")
        tanb   = st.number_input("Talud V:H  (V/H)", value=0.5, step=0.05, format="%.3f", key="tanb")
        Ae     = st.number_input("Cota base espaldón Ae (m)", value=0.5, step=0.1, format="%.2f", key="Ae")

    # ── 2. Oleaje y emplazamiento ──────────────────────────────────
    with st.expander("Condiciones de oleaje y emplazamiento", expanded=False):
        oleaje_origen = st.radio(
            "Origen del oleaje",
            options=["simar", "pie_dique"],
            format_func=lambda x: (
                "Hs₀ en aguas profundas (SIMAR → propaga a pie de dique)"
                if x == "simar"
                else "Hs a pie de dique (directo, sin propagación)"
            ),
            index=0,
            key="oleaje_origen",
        )
        hs0_label = "Hs₀ en SIMAR (m)" if oleaje_origen == "simar" else "Hs a pie de dique (m)"
        Hs0 = st.number_input(hs0_label, value=9.11, step=0.01, format="%.3f", key="Hs0")
        T   = st.number_input("Periodo T (s)", value=15.76, step=0.01, format="%.3f", key="T")
        if oleaje_origen == "simar":
            theta0 = st.number_input("Ángulo incidencia θ₀ (°)", value=0.0, step=1.0, format="%.1f", key="theta0")
        else:
            theta0 = 0.0
        h = st.number_input("Calado pie de dique h (m)", value=25.0, step=0.5, format="%.2f", key="h")

    # ── 3. Tipo de pieza del manto ─────────────────────────────────
    with st.expander("Tipo de pieza del manto principal", expanded=False):
        tipo_pieza = st.radio(
            "Material del manto principal",
            options=["hormigon", "escollera"],
            format_func=lambda x: (
                "Bloque de hormigón"
                if x == "hormigon"
                else "Escollera natural"
            ),
            index=0,
            key="tipo_pieza",
        )

    # ── 4. Iribarren y función de estabilidad ──────────────────────
    with st.expander("Índice de Iribarren y función de estabilidad", expanded=False):
        # Cálculo automático de Ir (sin necesitar φ)
        try:
            _ir_p = {
                "T": T, "h": h, "Hs0": Hs0, "tanb": tanb,
                "theta0": theta0, "oleaje_origen": oleaje_origen,
            }
            _Ir_prev, _Hp_prev, _L_prev = calcular_ir(_ir_p)
            st.info(
                f"**Ir = {_Ir_prev:.4f}**  ·  Hp = {_Hp_prev:.3f} m  ·  L = {_L_prev:.2f} m"
            )
            st.caption("Con este Ir, lee φ del ábaco de Losada et al. (1986).")
        except Exception as _e:
            st.warning(f"Ir no calculable: {_e}")

        phi = st.number_input(
            "φ — función de estabilidad",
            value=None,
            min_value=0.001,
            max_value=1.0,
            step=0.001,
            format="%.4f",
            placeholder="Introduce φ del ábaco",
            key="phi",
        )

    # ── 5. Berma de pie ───────────────────────────────────────────
    with st.expander("Berma de pie", expanded=False):
        con_berma = st.checkbox("Diseñar berma de pie", value=False, key="con_berma")
        hb = 5.0
        if con_berma:
            hb = st.number_input("Calado sobre berma hb (m)", value=5.0, step=0.1, format="%.2f", key="hb")

    # ── 6. Francobordo ────────────────────────────────────────────
    with st.expander("Cota de francobordo", expanded=False):
        aux_fc = st.radio(
            "Método de cálculo",
            options=[1, 2, 3],
            format_func=lambda x: {
                1: "(1) Probabilidad de ocurrencia de run-up",
                2: "(2) Sin rebase — run-up al 100 %",
                3: "(3) Caudal de rebase admisible q_adm",
            }[x],
            index=2,
            key="aux_fc",
        )
        qmax = 0.001
        prob = 0.068
        if aux_fc == 3:
            qmax = st.number_input(
                "q_adm (m³/s·ml)", value=0.001, step=0.0001, format="%.4f", key="qmax"
            )
        elif aux_fc == 1:
            prob = st.number_input(
                "Probabilidad de ocurrencia", value=0.068, step=0.001, format="%.4f", key="prob"
            )

    # ── 7. Espaldón ───────────────────────────────────────────────
    with st.expander("Espaldón", expanded=False):
        Fm       = st.number_input("Cota muelle Fm (m)", value=3.0, step=0.1, format="%.2f", key="Fm")
        CS       = st.number_input("Coef. seguridad CS", value=1.2, step=0.05, format="%.2f", key="CS")
        mu       = st.number_input("Coef. rozamiento μ", value=0.6, step=0.05, format="%.3f", key="mu")
        botaolas = st.number_input("Divisor base espaldón", value=4, step=1, min_value=1, key="botaolas")

# ═══════════════════════════════════════════════════════════════════
#  ZONA PRINCIPAL — BOTÓN Y RESULTADOS
# ═══════════════════════════════════════════════════════════════════

col_btn, col_sp = st.columns([1, 3])
with col_btn:
    calcular_btn = st.button("⚙️  CALCULAR", use_container_width=True)

# ── Botón "Acerca de" fijo en esquina inferior derecha ──
st.markdown("""
<div style="position:fixed; bottom:24px; right:28px; z-index:9999;">
  <details>
    <summary style="
        background:#0d2b55; color:white;
        padding:8px 18px; border-radius:20px;
        cursor:pointer; list-style:none;
        font-weight:600; font-size:0.88rem;
        box-shadow:0 3px 10px rgba(0,0,0,.25);
        user-select:none;">
      ℹ️&nbsp; Acerca de
    </summary>
    <div style="
        background:white; border:1px solid #d0d8e4;
        border-radius:10px; padding:1rem 1.2rem;
        margin-top:8px; min-width:230px;
        box-shadow:0 6px 20px rgba(0,0,0,.15);
        font-size:0.9rem; line-height:1.7;
        position:relative; bottom:0; right:0;">
      <b style="color:#0d2b55;">Autores base</b><br>
      R. Bergillos, M. Ortega, M. Clavero<br>
      2018 – 2025<br><br>
      <b style="color:#0d2b55;">Integración web</b><br>
      M. Clavero Gilabert<br>
      2026<br><br>
      <span style="color:#6b7a8d;">GDFA · IISTA · UGR</span>
    </div>
  </details>
</div>
""", unsafe_allow_html=True)

if calcular_btn:
    # ── Validación φ ──
    if phi is None:
        st.error(
            "**Dato obligatorio — φ no introducido.**\n\n"
            "Calcula el Índice de Iribarren (mostrado en la sección correspondiente "
            "de la barra lateral) y consulta el ábaco de Losada et al. (1986) para "
            "obtener el valor de φ (rango típico: 0.04 – 0.09)."
        )
        st.stop()

    # ── Construcción del dict de parámetros ──
    params = {
        "gammaw": gammaw, "gammas": gammas, "gammae": gammae,
        "tanb": tanb, "Ae": Ae,
        "Hs0": Hs0, "T": T, "theta0": theta0, "h": h,
        "oleaje_origen": oleaje_origen,
        "phi": phi,
        "tipo_pieza": tipo_pieza,
        "con_berma": con_berma, "hb": hb,
        "aux_fc": aux_fc, "qmax": qmax, "prob": prob,
        "Fm": Fm, "CS": CS, "mu": mu, "botaolas": int(botaolas),
    }

    # ── Cálculo ──
    try:
        r = calcular(params)
    except ValueError as e:
        st.error(f"**Error de cálculo:** {e}")
        st.stop()
    except Exception as e:
        st.error(f"**Error inesperado:** {e}")
        st.stop()

    st.success("✅ Cálculo completado")
    st.markdown("---")

    # ══════════════════════════════════════════
    #  BLOQUE 1 — OLEAJE A PIE DE DIQUE
    # ══════════════════════════════════════════
    st.markdown("### Oleaje a pie de dique")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Hs pie de dique", f"{r['Hs']:.3f} m")
    c2.metric("Hp diseño", f"{r['Hp']:.3f} m",
              help="Mínimo de Hmax = 1.6·Hs y altura de rotura Hrot")
    c3.metric("Longitud de onda L", f"{r['L']:.2f} m")
    c4.metric("Rotura", "Sí ⚠" if r["rotura"] else "No ✓",
              help="Hp limitada por rotura si 'Sí'")

    if oleaje_origen == "simar":
        st.caption(
            f"Propagación desde SIMAR: Ks = {r.get('Ks', 1.0):.4f}  ·  "
            f"Kr = {r.get('Kr', 1.0):.4f}  →  Hs = Hs₀·Ks·Kr = {r['Hs']:.3f} m"
        )

    # ══════════════════════════════════════════
    #  BLOQUE 2 — IRIBARREN Y REMONTE
    # ══════════════════════════════════════════
    st.markdown("### Iribarren y remonte")
    c1, c2, c3 = st.columns(3)
    c1.metric("Índice de Iribarren Ir", f"{r['Ir']:.4f}")
    c2.metric("Descenso de remonte Rd", f"{r['Rd']:.3f} m")
    c3.metric("Remonte máximo Rumax", f"{r['Rumax']:.3f} m")

    # ══════════════════════════════════════════
    #  BLOQUE 3 — MANTO PRINCIPAL
    # ══════════════════════════════════════════
    st.markdown("### Manto principal y secundarios")
    c1, c2 = st.columns(2)
    _unidad = "hormigón" if tipo_pieza == "hormigon" else "escollera"
    c1.metric(f"Peso pieza manto Wp", f"{r['Wp']:.4f} t",
              help=f"Piezas de {_unidad}")
    c2.metric("Espesor manto ep", f"{r['ep']:.3f} m")

    # Mantos secundarios — mismo formato que manto principal
    for i in range(4):
        if r["Ws"][i] > 0:
            cs1, cs2 = st.columns(2)
            cs1.metric(f"Peso manto secundario {i+1} Ws", f"{r['Ws'][i]:.4f} t")
            cs2.metric(f"Espesor manto secundario {i+1} es", f"{r['es'][i]:.3f} m")

    # ══════════════════════════════════════════
    #  BLOQUE 4 — GEOMETRÍA DE CORONACIÓN
    # ══════════════════════════════════════════
    st.markdown("### Geometría de coronación")
    c1, c2, c3 = st.columns(3)
    c1.metric("Cota coronación Ac", f"{r['Ac']:.3f} m")
    c2.metric("Francobordo Fc", f"{r['Fc']:.1f} m")
    c3.metric("Berma coronación Bc", f"{r['Bc']:.3f} m")
    c_so, _ = st.columns([1, 2])
    c_so.metric("So (run-up sobre Ac)", f"{r['So']:.3f} m",
                help="Si So ≤ 0 el run-up no alcanza la coronación")

    # ══════════════════════════════════════════
    #  BLOQUE 5 — BERMA DE PIE (opcional)
    # ══════════════════════════════════════════
    if con_berma and r["Wb"] is not None:
        st.markdown("### Berma de pie")
        c1, c2 = st.columns(2)
        c1.metric("Peso pieza berma Wb", f"{r['Wb']:.4f} t")
        c2.metric("Diámetro pieza Db", f"{r['Db']:.3f} m")

    # ══════════════════════════════════════════
    #  BLOQUE 6 — ESPALDÓN
    # ══════════════════════════════════════════
    st.markdown("### Espaldón")

    if r["adv_espaldon"]:
        st.warning(f"⚠ {r['adv_espaldon']}")
    else:
        c1, c2, c3 = st.columns(3)
        c1.metric("B desliz. dinámico", f"{r['B_p_desliz']:.3f} m" if r["B_p_desliz"] else "—")
        c2.metric("B vuelco dinámico",  f"{r['B_p_vuelco']:.3f} m" if r["B_p_vuelco"] else "—")
        c3.metric("Anchura mínima Bmax", f"{r['Bmax']:.3f} m"      if r["Bmax"]       else "—")
        c1, c2, c3 = st.columns(3)
        c1.metric("B desliz. reflectivo", f"{r['B_r_desliz']:.3f} m" if r["B_r_desliz"] else "—")
        c2.metric("B vuelco reflectivo",  f"{r['B_r_vuelco']:.3f} m" if r["B_r_vuelco"] else "—")
        c3.metric("Coef. presión μ",      f"{r['coef_mu']:.4f}"      if r["coef_mu"]     else "—")

    # ══════════════════════════════════════════
    #  RESUMEN COMPLETO (tabla)
    # ══════════════════════════════════════════
    with st.expander("Tabla resumen de todos los resultados", expanded=False):
        import pandas as pd

        filas = [
            ("Oleaje", "Hs a pie de dique (m)",         f"{r['Hs']:.4f}"),
            ("Oleaje", "Hp de diseño (m)",               f"{r['Hp']:.4f}"),
            ("Oleaje", "Longitud de onda L (m)",         f"{r['L']:.4f}"),
            ("Oleaje", "Rotura",                         "Sí" if r["rotura"] else "No"),
            ("Iribarren", "Ir",                          f"{r['Ir']:.6f}"),
            ("Iribarren", "Remonte máx. Rumax (m)",      f"{r['Rumax']:.4f}"),
            ("Iribarren", "Descenso remonte Rd (m)",     f"{r['Rd']:.4f}"),
            ("Manto", "Peso Wp (t)",                     f"{r['Wp']:.6f}"),
            ("Manto", "Espesor ep (m)",                  f"{r['ep']:.4f}"),
            ("Coronación", "Cota Ac (m)",                f"{r['Ac']:.4f}"),
            ("Coronación", "Francobordo Fc (m)",         f"{r['Fc']:.2f}"),
            ("Coronación", "Berma coronación Bc (m)",    f"{r['Bc']:.4f}"),
            ("Coronación", "So (m)",                     f"{r['So']:.4f}"),
        ]
        if r["Wb"]:
            filas += [
                ("Berma pie", "Peso Wb (t)",  f"{r['Wb']:.6f}"),
                ("Berma pie", "Diám. Db (m)", f"{r['Db']:.4f}"),
            ]
        if r["Bmax"]:
            filas += [
                ("Espaldón", "B desliz. dinámico (m)", f"{r['B_p_desliz']:.4f}"),
                ("Espaldón", "B vuelco dinámico (m)",  f"{r['B_p_vuelco']:.4f}"),
                ("Espaldón", "B desliz. reflectivo (m)", f"{r['B_r_desliz']:.4f}"),
                ("Espaldón", "B vuelco reflectivo (m)",  f"{r['B_r_vuelco']:.4f}"),
                ("Espaldón", "Bmax — anchura espaldón (m)", f"{r['Bmax']:.4f}"),
            ]

        filas_html = "".join(
            f"<tr><td>{b}</td><td>{p}</td><td>{v}</td></tr>"
            for b, p, v in filas
        )
        st.markdown(f"""
<style>
.tabla-resumen {{
    width: 100%;
    border-collapse: collapse;
    font-size: 20px;
}}
.tabla-resumen th {{
    background-color: #0d2b55;
    color: white;
    padding: 8px 12px;
    text-align: left;
    font-size: 20px;
}}
.tabla-resumen td {{
    padding: 6px 12px;
    border-bottom: 1px solid #d0d8e4;
    font-size: 20px;
}}
.tabla-resumen tr:nth-child(even) td {{
    background-color: #f2f4f7;
}}
.tabla-resumen tr:hover td {{
    background-color: #dce4f0;
}}
</style>
<table class="tabla-resumen">
  <thead><tr><th>Bloque</th><th>Parámetro</th><th>Valor</th></tr></thead>
  <tbody>{filas_html}</tbody>
</table>
""", unsafe_allow_html=True)

    # ══════════════════════════════════════════
    #  REFERENCIAS
    # ══════════════════════════════════════════
    st.markdown("---")
    st.markdown("### Referencias empleadas en el dimensionamiento")

    st.markdown("""
**Ascenso y remonte**
- Losada, M.A. & Giménez-Curto, L.A. (1981). Flow characteristics on rough, permeable slopes under wave action. *Coastal Engineering*, 4, 187–206.
- Losada, M.A., Desire, J.M. & Alejo, L.M. (1986). Stability of blocks as breakwater armoring. *Journal of Structural Engineering*, ASCE, 112(11), 2392–2401.
- EurOtop (2018). *Manual on wave overtopping of sea defences and related structures.* Van der Meer, J.W. et al.

**Manto principal**
- Losada, M.A. & Giménez-Curto, L.A. (1979). The joint effect of the wave height and period on the stability of rubble mound breakwaters using Ir. *Coastal Engineering*, 3, 77–96.
- CIRIA, CUR, CETMEF (2007). *The Rock Manual. The use of rock in hydraulic engineering* (2nd ed.). C683, CIRIA, London.

**Berma de pie**
- Muttray, M. (2013). A pragmatic approach to rock toe design. *Coastal Engineering*, 82, 56–63.

**Cota de francobordo**
- Losada, M.A., Desire, J.M. & Alejo, L.M. (1986). Stability of blocks as breakwater armoring. *Journal of Structural Engineering*, ASCE, 112(11), 2392–2401.
- EurOtop (2018). *Manual on wave overtopping of sea defences and related structures.* Van der Meer, J.W. et al.

**Estabilidad del espaldón**
- Martín, F.L., Losada, M.A. & Medina, R. (1999). Wave loads on rubble mound breakwater crown walls. *Coastal Engineering*, 37, 149–174.
""")

else:
    # ── Estado inicial — instrucciones ──
    st.info(
        "👈 **Introduce los datos en la barra lateral** y pulsa **⚙️ CALCULAR**.\n\n"
        "Las secciones están plegadas; haz clic en cada una para expandirla. "
        "El Índice de Iribarren se calcula automáticamente al introducir los datos "
        "de oleaje — úsalo para leer φ del ábaco de Losada et al. (1986)."
    )
    st.markdown("""
---
#### Flujo de trabajo recomendado

1. **Geometría y materiales** — pesos específicos y talud.
2. **Condiciones de oleaje** — elige si propagas desde SIMAR o das Hs a pie de dique.
3. **Tipo de pieza** — hormigón o escollera.
4. **Iribarren / φ** — el Ir se muestra automáticamente; lee φ del ábaco y anótalo.
5. **Berma de pie**, **Francobordo**, **Espaldón** — completa o deja valores por defecto.
6. Pulsa **⚙️ CALCULAR**.
""")
