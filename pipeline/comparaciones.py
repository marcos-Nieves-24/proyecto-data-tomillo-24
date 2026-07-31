"""Fase 8: comparaciones multiples entre tecnicas de extraccion.

Dependiendo del modelo seleccionado: Tukey HSD (via parametrica) o test de
Dunn con correccion FDR (via no parametrica), mas Wilcoxon (Mann-Whitney)
con FDR como robustez. Genera grupos homogeneos con letras compactas (CLD)
y la figura correspondiente.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as st
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

from pipeline.config import (
    METODOS,
    METODO_LABEL,
    PALETA_METODOS,
    TAMANO_FIG,
    VARIABLE_LABEL,
    guardar_tabla,
    save_figure_pub,
)

ALPHA = 0.05


# ---------------------------------------------------------------------------
# Tests no parametricos
# ---------------------------------------------------------------------------


def test_dunn(df: pd.DataFrame, variable: str) -> tuple[pd.DataFrame, dict]:
    """Test de Dunn (no parametrico, basado en rangos) con correccion FDR."""
    d = df.copy()
    d["rango"] = d[variable].rank(method="average")
    n = len(d)
    rangos = {m: d.loc[d["metodo_extraccion"] == m, "rango"].values for m in METODOS}
    tamanos = {m: len(v) for m, v in rangos.items()}
    rangos_medios = {m: float(v.mean()) for m, v in rangos.items()}

    conteos = d[variable].value_counts()
    ties = float(((conteos ** 3) - conteos).sum())
    varianza = (n * (n + 1)) / 12.0 - ties / (12.0 * (n - 1))

    pares = []
    for i, m1 in enumerate(METODOS):
        for m2 in METODOS[i + 1:]:
            diff = rangos_medios[m1] - rangos_medios[m2]
            se = np.sqrt(varianza * (1.0 / tamanos[m1] + 1.0 / tamanos[m2]))
            z = diff / se if se > 0 else np.nan
            p = 2 * (1 - st.norm.cdf(abs(z))) if np.isfinite(z) else np.nan
            pares.append({"par": f"{m1} vs {m2}", "estadistico_z": z, "p_valor": p})
    out = pd.DataFrame(pares)
    if out["p_valor"].notna().any():
        _, p_fdr, _, _ = multipletests(out["p_valor"].fillna(1.0), method="fdr_bh")
        out["p_valor_ajustado"] = p_fdr
    else:
        out["p_valor_ajustado"] = np.nan
    return out, rangos_medios


def _wilcoxon_fdr(df: pd.DataFrame, variable: str) -> pd.DataFrame:
    """Wilcoxon (Mann-Whitney U) por pares con correccion FDR como robustez."""
    pares = []
    for i, m1 in enumerate(METODOS):
        for m2 in METODOS[i + 1:]:
            a = df.loc[df["metodo_extraccion"] == m1, variable].values
            b = df.loc[df["metodo_extraccion"] == m2, variable].values
            try:
                p = st.mannwhitneyu(a, b, alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
            pares.append({"par": f"{m1} vs {m2}", "p_valor": p})
    out = pd.DataFrame(pares)
    if out["p_valor"].notna().any():
        _, p_fdr, _, _ = multipletests(out["p_valor"].fillna(1.0), method="fdr_bh")
        out["p_valor_ajustado"] = p_fdr
    else:
        out["p_valor_ajustado"] = np.nan
    return out


def _tukey_hsd(df: pd.DataFrame, variable: str) -> pd.DataFrame:
    """Tukey HSD entre metodos (via parametrica)."""
    tukey = pairwise_tukeyhsd(df[variable], df["metodo_extraccion"], alpha=ALPHA)
    grupos = list(tukey.groupsunique)
    pares = [(a, b) for i, a in enumerate(grupos) for b in grupos[i + 1:]]
    filas = []
    for i, (a, b) in enumerate(pares):
        filas.append({
            "par": f"{a} vs {b}",
            "diferencia_medias": float(tukey.meandiffs[i]),
            "p_valor_ajustado": float(tukey.pvalues[i]),
            "ic95_inferior": float(tukey.confint[i, 0]),
            "ic95_superior": float(tukey.confint[i, 1]),
            "significativo": bool(tukey.reject[i]),
        })
    return pd.DataFrame(filas)


# ---------------------------------------------------------------------------
# Letras compactas (CLD)
# ---------------------------------------------------------------------------


def compact_letter_display(
    df: pd.DataFrame, variable: str, tabla_pares: pd.DataFrame
) -> dict:
    """Asigna letras compactas (CLD) a los metodos.

    Dos metodos comparten al menos una letra si y solo si NO difieren
    significativamente (p ajustada >= alpha). Algoritmo greedy sobre grupos
    ordenados por media descendente, con correccion final para pares
    significativos.
    """
    medias = df.groupby("metodo_extraccion")[variable].mean().sort_values(ascending=False)
    grupos = list(medias.index)
    pmat = {}
    for _, fila in tabla_pares.iterrows():
        m1, m2 = fila["par"].split(" vs ")
        pmat[(m1, m2)] = float(fila["p_valor_ajustado"])
        pmat[(m2, m1)] = float(fila["p_valor_ajustado"])

    letras: dict[str, list] = {g: [] for g in grupos}
    abc = "abcdefghijklmnopqrstuvwxyz"
    letra_actual = 0

    for i, g in enumerate(grupos):
        compartida = False
        for j in range(i):
            anterior = grupos[j]
            if pmat.get((g, anterior), 1.0) >= ALPHA and letras[anterior]:
                letras[g].append(letras[anterior][0])
                compartida = True
                break
        if not compartida:
            letras[g].append(abc[letra_actual])
            letra_actual += 1

    # Correccion: pares significativos no deben compartir ninguna letra.
    for i in range(len(grupos)):
        for j in range(i):
            m1, m2 = grupos[i], grupos[j]
            if pmat.get((m1, m2), 1.0) < ALPHA:
                comunes = set(letras[m1]) & set(letras[m2])
                if comunes:
                    letras[m2].append(abc[letra_actual])
                    letra_actual += 1

    return {g: "".join(sorted(set(l))) for g, l in letras.items()}


# ---------------------------------------------------------------------------
# Orquestador
# ---------------------------------------------------------------------------


def comparaciones_posthoc(df: pd.DataFrame, variable: str, tipo_modelo: str) -> dict:
    """Ejecuta comparaciones multiples entre metodos segun el modelo.

    Para modelos OLS usa Tukey HSD; para la via no parametrica usa Dunn (FDR)
    con Wilcoxon (FDR) como robustez. Genera letras CLD, guarda tablas y
    la figura ``posthoc_<variable>_letras``.
    """
    if tipo_modelo in ("factorial_ols", "anova_ols"):
        tabla_pares = _tukey_hsd(df, variable)
        metodo = "Tukey HSD"
        robustez = None
    else:
        tabla_pares, rangos_medios = test_dunn(df, variable)
        robustez = _wilcoxon_fdr(df, variable)
        metodo = "Dunn (FDR)"

    letras = compact_letter_display(df, variable, tabla_pares)

    tabla_pares = tabla_pares.copy()
    tabla_pares["letras"] = tabla_pares["par"].map(
        lambda par: _letras_del_par(par, letras)
    )
    guardar_tabla(tabla_pares, f"posthoc_{variable}", index=False)
    guardar_tabla(
        pd.DataFrame([
            {"metodo_extraccion": m, "media": round(float(df.loc[df['metodo_extraccion'] == m, variable].mean()), 3),
             "letras": letras[m]}
            for m in METODOS
        ]),
        f"posthoc_{variable}_letras", index=False,
    )
    if robustez is not None:
        guardar_tabla(robustez, f"posthoc_{variable}_wilcoxon", index=False)

    # Figura con letras CLD
    medias = df.groupby("metodo_extraccion")[variable].mean()
    err = df.groupby("metodo_extraccion")[variable].sem()
    orden = [m for m in METODOS if m in medias.index]
    x = np.arange(len(orden))
    colores = [PALETA_METODOS[m] for m in orden]

    fig, ax = plt.subplots(figsize=TAMANO_FIG)
    ax.bar(x, medias[orden], yerr=err[orden], color=colores, capsize=4, alpha=0.85)
    ax.set_xticks(x, [METODO_LABEL[m] for m in orden])
    ax.set_ylabel(VARIABLE_LABEL[variable])
    rango = float(medias.max() - medias.min()) if medias.max() != medias.min() else 1.0
    for xi, m in enumerate(orden):
        ax.text(
            xi, float(medias[m]) + float(err[m]) + 0.03 * rango,
            letras[m], ha="center", fontweight="bold", fontsize=13,
        )
    save_figure_pub(
        fig, f"posthoc_{variable}_letras",
        titulo=f"Comparaciones multiples - {VARIABLE_LABEL[variable]}",
    )

    print("=" * 72)
    print(f"FASE 8 - Post-hoc ({metodo}): {VARIABLE_LABEL[variable]}")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print(tabla_pares.to_string(index=False))
    print(f"Letras CLD: {letras}")

    return {
        "variable": variable,
        "metodo": metodo,
        "tabla_pares": tabla_pares,
        "letras": letras,
        "robustez": robustez,
    }


def _letras_del_par(par: str, letras: dict) -> str:
    m1, m2 = par.split(" vs ")
    comunes = sorted(set(letras[m1]) & set(letras[m2]))
    if comunes:
        return "".join(comunes)
    return "sin letra comun"
