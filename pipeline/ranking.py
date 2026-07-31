"""Fase 11: ranking de tecnicas de extraccion.

Combina el rendimiento medio, la inhibicion micelial media y la inhibicion
de conidias media (cada una normalizada min-max 0-1) en un score compuesto
(promedio simple) y genera el ranking con una tabla comparativa y un grafico
radar de tres ejes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

from pipeline.config import (
    METODOS,
    METODO_LABEL,
    PALETA_METODOS,
    guardar_tabla,
    save_figure_pub,
)


def _normalizar_minmax(serie: pd.Series) -> pd.Series:
    """Normaliza min-max a 0-1 (mayor valor -> 1)."""
    rango = serie.max() - serie.min()
    if rango == 0:
        return pd.Series(0.5, index=serie.index)
    return (serie - serie.min()) / rango


def ranking_tecnicas(df_bio: pd.DataFrame, df_rend: pd.DataFrame) -> dict:
    """Construye el ranking de tecnicas por rendimiento y actividad antifungica.

    Metrica por tecnica:
      - rendimiento medio (%) sobre las 3 replicas,
      - %INH micelial medio sobre todos los aislados,
      - %INH conidias medio sobre todos los aislados.
    Normaliza cada metrica a 0-1, promedia y ordena.
    """
    rend_medio = df_rend.groupby("metodo_extraccion")["rendimiento_pct"].mean()
    inh_mic_medio = df_bio.groupby("metodo_extraccion")["porcentaje_inhibicion_micelial"].mean()
    inh_con_medio = df_bio.groupby("metodo_extraccion")["porcentaje_inhibicion_conidias"].mean()

    tabla = pd.DataFrame({
        "metodo_extraccion": METODOS,
        "rendimiento_medio_pct": rend_medio.reindex(METODOS).round(3).values,
        "inhib_micelial_medio_pct": inh_mic_medio.reindex(METODOS).round(3).values,
        "inhib_conidias_medio_pct": inh_con_medio.reindex(METODOS).round(3).values,
    })
    tabla["rendimiento_norm"] = _normalizar_minmax(tabla["rendimiento_medio_pct"]).round(3).values
    tabla["inhib_micelial_norm"] = _normalizar_minmax(tabla["inhib_micelial_medio_pct"]).round(3).values
    tabla["inhib_conidias_norm"] = _normalizar_minmax(tabla["inhib_conidias_medio_pct"]).round(3).values
    tabla["score_compuesto"] = (
        tabla[["rendimiento_norm", "inhib_micelial_norm", "inhib_conidias_norm"]]
        .mean(axis=1)
        .round(3)
    )
    tabla = tabla.sort_values("score_compuesto", ascending=False).reset_index(drop=True)
    tabla.insert(0, "ranking", range(1, len(tabla) + 1))
    guardar_tabla(tabla, "ranking_tecnicas", index=False)

    # Radar de tres ejes
    metricas = ["rendimiento_norm", "inhib_micelial_norm", "inhib_conidias_norm"]
    etiquetas_ejes = ["Rendimiento", "INH micelial", "INH conidias"]
    angulos = np.linspace(0, 2 * np.pi, len(metricas), endpoint=False).tolist()
    angulos += angulos[:1]

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
    for m in METODOS:
        valores = [tabla.loc[tabla["metodo_extraccion"] == m, met].iloc[0] for met in metricas]
        valores += valores[:1]
        ax.plot(angulos, valores, lw=2, color=PALETA_METODOS[m],
                label=METODO_LABEL[m], marker="o", ms=5)
        ax.fill(angulos, valores, color=PALETA_METODOS[m], alpha=0.12)
    ax.set_xticks(angulos[:-1])
    ax.set_xticklabels(etiquetas_ejes)
    ax.set_ylim(0, 1)
    ax.set_title("Ranking de tecnicas (metricas normalizadas 0-1)")
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.05))
    save_figure_pub(fig, "ranking_radar", titulo="Ranking de técnicas de extracción (radar)")

    print("=" * 72)
    print("FASE 11 - Ranking de técnicas de extracción")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print(tabla.to_string(index=False))

    return {"tabla": tabla, "metricas": metricas}
