"""Fase 4: analisis exploratorio de datos (EDA).

Resumenes descriptivos por método para cada variable de respuesta y figuras
exploratorias de calidad de publicacion (histogramas, boxplots, violines,
densidades con rug, QQ-plots, matriz de correlacion y scatter).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns

from pipeline.config import (
    METODOS,
    METODO_LABEL,
    PALETA_METODOS,
    TAMANO_FIG,
    VARIABLES_RESPUESTA,
    VARIABLE_LABEL,
    guardar_tabla,
    save_figure_pub,
)


def _ic95(serie: pd.Series) -> tuple[float, float]:
    """Intervalo de confianza del 95 % para la media (t de Student)."""
    n = serie.count()
    if n < 2:
        return (np.nan, np.nan)
    media = serie.mean()
    se = serie.std(ddof=1) / np.sqrt(n)
    inf, sup = st.t.interval(0.95, df=n - 1, loc=media, scale=se)
    return inf, sup


def resumen_descriptivo(dfs: pd.DataFrame) -> pd.DataFrame:
    """Calcula estadisticos descriptivos por método para cada respuesta.

    Si el dataset incluye los controles C4 (columnas ``control_*``), agrega
    filas de línea de base con su media y desviación. Devuelve y guarda
    ``eda_descriptivos.csv``.
    """
    filas = []
    for metodo in METODOS:
        grupo = dfs[dfs["metodo_extraccion"] == metodo]
        for var in VARIABLES_RESPUESTA:
            s = grupo[var].dropna()
            n = s.count()
            if n == 0:
                continue
            media = float(s.mean())
            desv = float(s.std(ddof=1)) if n > 1 else np.nan
            err = desv / np.sqrt(n) if n > 1 else np.nan
            inf, sup = _ic95(s)
            filas.append({
                "metodo_extraccion": METODO_LABEL[metodo],
                "variable": VARIABLE_LABEL[var],
                "n": int(n),
                "media": round(media, 3),
                "desviacion_estandar": round(desv, 3) if np.isfinite(desv) else np.nan,
                "error_estandar": round(err, 3) if np.isfinite(err) else np.nan,
                "ic95_inferior": round(inf, 3) if np.isfinite(inf) else np.nan,
                "ic95_superior": round(sup, 3) if np.isfinite(sup) else np.nan,
                "minimo": round(float(s.min()), 3),
                "maximo": round(float(s.max()), 3),
            })
        # Línea de base: controles C4 del aislado (compartidos por las 3 réplicas).
        for var in ("control_crecimiento_mm", "control_conidias_log10"):
            if var not in dfs.columns:
                continue
            s = grupo[var].dropna()
            n = s.count()
            if n == 0:
                continue
            etiqueta = "Control C4: crecimiento micelial (mm)" if var == "control_crecimiento_mm" else "Control C4: conidias (log10/mL)"
            media = float(s.mean())
            desv = float(s.std(ddof=1)) if n > 1 else np.nan
            err = desv / np.sqrt(n) if n > 1 else np.nan
            inf, sup = _ic95(s)
            filas.append({
                "metodo_extraccion": METODO_LABEL[metodo],
                "variable": etiqueta,
                "n": int(n),
                "media": round(media, 3),
                "desviacion_estandar": round(desv, 3) if np.isfinite(desv) else np.nan,
                "error_estandar": round(err, 3) if np.isfinite(err) else np.nan,
                "ic95_inferior": round(inf, 3) if np.isfinite(inf) else np.nan,
                "ic95_superior": round(sup, 3) if np.isfinite(sup) else np.nan,
                "minimo": round(float(s.min()), 3),
                "maximo": round(float(s.max()), 3),
            })
    tabla = pd.DataFrame(filas)
    guardar_tabla(tabla, "eda_descriptivos", index=False)

    print("=" * 72)
    print("FASE 4 - Resumen descriptivo por método")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(tabla.to_string(index=False))

    return tabla


def figuras_eda(dfs: pd.DataFrame) -> list[str]:
    """Genera las figuras exploratorias con prefijo 'eda_' (PNG y PDF)."""
    generadas: list[str] = []
    sns.set_theme(style="whitegrid")

    df = dfs.copy()
    df["metodo_label"] = df["metodo_extraccion"].map(METODO_LABEL)
    orden = [METODO_LABEL[m] for m in METODOS]
    paleta = {METODO_LABEL[m]: PALETA_METODOS[m] for m in METODOS}

    # 1. Histogramas por método y variable ---------------------------------
    for var in VARIABLES_RESPUESTA:
        fig, ejes = plt.subplots(1, 3, figsize=(16, 4.5))
        for ax, metodo in zip(ejes, orden):
            datos = df.loc[df["metodo_label"] == metodo, var].dropna()
            ax.hist(datos, bins=15, color=paleta[metodo], edgecolor="white", alpha=0.85)
            ax.set_title(metodo)
            ax.set_xlabel(VARIABLE_LABEL[var])
            ax.set_ylabel("Frecuencia")
        save_figure_pub(
            fig, f"eda_histogramas_{var}",
            titulo=f"Histogramas por método - {VARIABLE_LABEL[var]}",
        )
        generadas.append(f"eda_histogramas_{var}")

    # 2. Boxplots por método -----------------------------------------------
    fig, ejes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, var in zip(ejes.ravel(), VARIABLES_RESPUESTA):
        sns.boxplot(
            data=df, x="metodo_label", y=var, order=orden, palette=paleta, ax=ax,
            hue="metodo_label", hue_order=orden, legend=False,
        )
        ax.set_title(VARIABLE_LABEL[var])
        ax.set_xlabel("")
    save_figure_pub(fig, "eda_boxplots_por_metodo", titulo="Boxplots por método")
    generadas.append("eda_boxplots_por_metodo")

    # 3. Violin plots por método -------------------------------------------
    fig, ejes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, var in zip(ejes.ravel(), VARIABLES_RESPUESTA):
        sns.violinplot(
            data=df, x="metodo_label", y=var, order=orden, palette=paleta, ax=ax,
            hue="metodo_label", hue_order=orden, legend=False, cut=0,
        )
        ax.set_title(VARIABLE_LABEL[var])
        ax.set_xlabel("")
    save_figure_pub(fig, "eda_violines_por_metodo", titulo="Violin plots por método")
    generadas.append("eda_violines_por_metodo")

    # 4. Densidad con rug por variable --------------------------------------
    for var in VARIABLES_RESPUESTA:
        fig, ax = plt.subplots(figsize=TAMANO_FIG)
        for metodo in METODOS:
            datos = df.loc[df["metodo_extraccion"] == metodo, var].dropna()
            sns.kdeplot(
                datos, ax=ax, label=METODO_LABEL[metodo],
                color=PALETA_METODOS[metodo], fill=True, alpha=0.25,
            )
            ax.plot(datos, np.full_like(datos, -0.02), "|", color=PALETA_METODOS[metodo], alpha=0.4)
        ax.set_xlabel(VARIABLE_LABEL[var])
        ax.set_ylabel("Densidad")
        ax.legend(title="Método")
        save_figure_pub(
            fig, f"eda_densidad_{var}", titulo=f"Densidad por método - {VARIABLE_LABEL[var]}"
        )
        generadas.append(f"eda_densidad_{var}")

    # 5. QQ-plot por variable y metodo ---------------------------------------
    for var in VARIABLES_RESPUESTA:
        fig, ejes = plt.subplots(1, 3, figsize=(16, 4.5))
        for ax, metodo in zip(ejes, orden):
            datos = df.loc[df["metodo_label"] == metodo, var].dropna()
            st.probplot(datos, dist="norm", plot=ax)
            ax.set_title(f"{metodo} - {VARIABLE_LABEL[var]}")
            ax.set_ylabel("Cuantiles muestrales")
        save_figure_pub(
            fig, f"eda_qqplot_{var}", titulo=f"QQ-plot por método - {VARIABLE_LABEL[var]}"
        )
        generadas.append(f"eda_qqplot_{var}")

    # 6. Matriz de correlacion entre respuestas ------------------------------
    vars_short = VARIABLES_RESPUESTA
    corr = df[vars_short].corr(method="pearson")
    pmat = pd.DataFrame(np.ones((4, 4)), index=vars_short, columns=vars_short)
    for i in range(4):
        for j in range(4):
            if i != j:
                pmat.iloc[i, j] = st.pearsonr(df[vars_short[i]], df[vars_short[j]]).pvalue

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    mascara = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(
        corr, mask=mascara, annot=True, fmt=".2f", cmap="RdBu_r", vmin=-1, vmax=1,
        square=True, cbar_kws={"label": "Pearson r"}, ax=ax,
        xticklabels=[VARIABLE_LABEL[v] for v in vars_short],
        yticklabels=[VARIABLE_LABEL[v] for v in vars_short],
    )
    for i in range(4):
        for j in range(i):
            ax.text(j + 0.5, i + 0.5, f"p={pmat.iloc[i, j]:.2g}", ha="center", va="center", fontsize=8, color="black")
    ax.set_title("Correlacion entre variables de respuesta (Pearson, con p)")
    fig.tight_layout()
    save_figure_pub(fig, "eda_correlacion", titulo="Matriz de correlación entre variables de respuesta")
    generadas.append("eda_correlacion")

    # 7. Scatter crecimiento vs conidias coloreado por método ------------------
    fig, ax = plt.subplots(figsize=TAMANO_FIG)
    for metodo in METODOS:
        sub = df[df["metodo_extraccion"] == metodo]
        ax.scatter(
            sub["crecimiento_micelial_mm"], sub["conidias_log10_ml"],
            s=30, alpha=0.65, color=PALETA_METODOS[metodo], label=METODO_LABEL[metodo],
        )
    ax.set_xlabel(VARIABLE_LABEL["crecimiento_micelial_mm"])
    ax.set_ylabel(VARIABLE_LABEL["conidias_log10_ml"])
    ax.legend(title="Método")
    save_figure_pub(
        fig, "eda_scatter_crecimiento_conidias",
        titulo="Crecimiento micelial vs conidias por método",
    )
    generadas.append("eda_scatter_crecimiento_conidias")

    print(f"  Figuras EDA generadas: {len(generadas)}")
    return generadas
