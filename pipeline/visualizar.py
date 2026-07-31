"""Fase 9: figuras de resultados principales.

Genera, por variable de respuesta: medias con SD/SE/IC95%, grafico de
interaccion metodo x aislado, efectos principales y comparacion con letras
CLD. Guarda ademas la tabla ``medias_<variable>.csv``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt

from pipeline.config import (
    METODOS,
    METODO_LABEL,
    PALETA_METODOS,
    VARIABLES_RESPUESTA,
    VARIABLE_LABEL,
    guardar_tabla,
    save_figure_pub,
)


def _medias_resumen(df: pd.DataFrame, variable: str) -> pd.DataFrame:
    """Resumen de medias por metodo con SD, SE e IC95%."""
    filas = []
    for m in METODOS:
        s = df.loc[df["metodo_extraccion"] == m, variable].dropna()
        n = s.count()
        media = float(s.mean())
        desv = float(s.std(ddof=1)) if n > 1 else np.nan
        err = desv / np.sqrt(n) if n > 1 else np.nan
        if n > 1:
            inf, sup = st.t.interval(0.95, df=n - 1, loc=media, scale=err)
        else:
            inf = sup = np.nan
        filas.append({
            "metodo_extraccion": m,
            "media": round(media, 3),
            "desviacion_estandar": round(desv, 3) if np.isfinite(desv) else np.nan,
            "error_estandar": round(err, 3) if np.isfinite(err) else np.nan,
            "ic95_inferior": round(inf, 3) if np.isfinite(inf) else np.nan,
            "ic95_superior": round(sup, 3) if np.isfinite(sup) else np.nan,
            "n": int(n),
        })
    return pd.DataFrame(filas)


def _panel_medias(df: pd.DataFrame, variable: str):
    """Figura con tres paneles: media +/- SD, media +/- SE, media con IC95%."""
    resumen = _medias_resumen(df, variable)
    medias = resumen.set_index("metodo_extraccion")
    orden = METODOS
    x = np.arange(len(orden))
    colores = [PALETA_METODOS[m] for m in orden]

    fig, ejes = plt.subplots(1, 3, figsize=(15, 4.8))
    paneles = [
        ("SD", "desviacion_estandar", "Media +/- DE"),
        ("SE", "error_estandar", "Media +/- EE"),
        ("IC95", None, "Media con IC95%"),
    ]
    for ax, (nombre, col, titulo) in zip(ejes, paneles):
        barras = ax.bar(
            x, medias["media"], color=colores, alpha=0.8,
            yerr=medias[col] if col else None,
            capsize=4, error_kw={"elinewidth": 1},
        )
        if col is None:
            ax.errorbar(
                x, medias["media"],
                yerr=[medias["media"] - medias["ic95_inferior"],
                      medias["ic95_superior"] - medias["media"]],
                fmt="none", ecolor="black", capsize=4,
            )
        ax.set_xticks(x, [METODO_LABEL[m] for m in orden])
        ax.set_ylabel(VARIABLE_LABEL[variable])
        ax.set_title(titulo)
        for b, m in zip(barras, orden):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                    f"{medias.loc[m, 'media']:.1f}", ha="center", va="bottom", fontsize=9)
    save_figure_pub(
        fig, f"resultados_{variable}_medias",
        titulo=f"Medias por metodo - {VARIABLE_LABEL[variable]}",
    )


def _figura_interaccion(df: pd.DataFrame, variable: str):
    """Grafico de interaccion metodo x aislado (medias por celda)."""
    orden_aislados = sorted(
        df["aislamiento"].unique(),
        key=lambda a: float(df.loc[df["aislamiento"] == a, variable].mean()),
        reverse=True,
    )
    x = np.arange(len(orden_aislados))
    fig, ax = plt.subplots(figsize=(13, 5.5))
    for m in METODOS:
        sub = df[df["metodo_extraccion"] == m]
        medias_celda = sub.groupby("aislamiento")[variable].mean().reindex(orden_aislados)
        ax.plot(x, medias_celda.values, marker="o", ms=4, lw=1.4,
                color=PALETA_METODOS[m], label=METODO_LABEL[m])
    ax.set_xticks(x, orden_aislados, rotation=90, fontsize=8)
    ax.set_xlabel("Aislado")
    ax.set_ylabel(VARIABLE_LABEL[variable])
    ax.legend(title="Metodo", fontsize=9)
    save_figure_pub(
        fig, f"resultados_{variable}_interaccion",
        titulo=f"Interacción método × aislado - {VARIABLE_LABEL[variable]}",
    )


def _figura_efectos(df: pd.DataFrame, variable: str):
    """Efectos principales: media marginal por metodo y por aislado."""
    fig, ejes = plt.subplots(1, 2, figsize=(14, 5))

    medias_metodo = df.groupby("metodo_extraccion")[variable].mean()
    err_metodo = df.groupby("metodo_extraccion")[variable].sem()
    x = np.arange(len(METODOS))
    ejes[0].bar(x, medias_metodo[METODOS], yerr=err_metodo[METODOS],
                color=[PALETA_METODOS[m] for m in METODOS], capsize=4, alpha=0.85)
    ejes[0].set_xticks(x, [METODO_LABEL[m] for m in METODOS])
    ejes[0].set_title("Efecto principal del metodo")
    ejes[0].set_ylabel(VARIABLE_LABEL[variable])

    medias_aislado = df.groupby("aislamiento")[variable].mean().sort_values(ascending=False)
    err_aislado = df.groupby("aislamiento")[variable].sem().reindex(medias_aislado.index)
    xa = np.arange(len(medias_aislado))
    ejes[1].bar(xa, medias_aislado.values, yerr=err_aislado.values,
                color="#6A737D", capsize=2, alpha=0.8)
    ejes[1].set_xticks(xa, medias_aislado.index, rotation=90, fontsize=8)
    ejes[1].set_title("Efecto principal del aislado")
    save_figure_pub(
        fig, f"resultados_{variable}_efectos_principales",
        titulo=f"Efectos principales - {VARIABLE_LABEL[variable]}",
    )


def _figura_letras(df: pd.DataFrame, variable: str, letras: dict):
    """Comparacion de medias con letras CLD superpuestas."""
    medias = df.groupby("metodo_extraccion")[variable].mean()
    err = df.groupby("metodo_extraccion")[variable].sem()
    orden = METODOS
    x = np.arange(len(orden))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(x, medias[orden], yerr=err[orden],
           color=[PALETA_METODOS[m] for m in orden], capsize=4, alpha=0.85)
    ax.set_xticks(x, [METODO_LABEL[m] for m in orden])
    ax.set_ylabel(VARIABLE_LABEL[variable])
    rango = float(medias.max() - medias.min()) if medias.max() != medias.min() else 1.0
    for xi, m in enumerate(orden):
        ax.text(xi, float(medias[m]) + float(err[m]) + 0.03 * rango,
                letras.get(m, ""), ha="center", fontweight="bold", fontsize=13)
    save_figure_pub(
        fig, f"resultados_{variable}_letras",
        titulo=f"Comparación con letras CLD - {VARIABLE_LABEL[variable]}",
    )


def figuras_resultados(df_bio: pd.DataFrame, posthoc: dict) -> list[str]:
    """Genera las figuras de resultados (prefijo 'resultados_').

    ``posthoc`` es un dict variable -> resultados de comparaciones_posthoc.
    """
    generadas: list[str] = []
    for variable in VARIABLES_RESPUESTA:
        resumen = _medias_resumen(df_bio, variable)
        guardar_tabla(resumen, f"medias_{variable}", index=False)

        _panel_medias(df_bio, variable)
        generadas.append(f"resultados_{variable}_medias")
        _figura_interaccion(df_bio, variable)
        generadas.append(f"resultados_{variable}_interaccion")
        _figura_efectos(df_bio, variable)
        generadas.append(f"resultados_{variable}_efectos_principales")

        if variable in posthoc:
            letras = posthoc[variable].get("letras", {})
            if letras:
                _figura_letras(df_bio, variable, letras)
                generadas.append(f"resultados_{variable}_letras")

    print(f"  Figuras de resultados generadas: {len(generadas)}")
    return generadas
