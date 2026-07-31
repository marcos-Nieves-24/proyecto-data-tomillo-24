"""Fase 3.3: Exploración descriptiva y distribución del rendimiento BDCA.

Genera resumenes descriptivos por tratamiento, medias e intervalos de confianza del 95 %
(reutilizando el patrón _ic95 del pipeline DCA cuando existe) y figuras exploratorias.
"""

from __future__ import annotations

import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns

from pipeline.config import (
    guardar_tabla,
    save_figure_pub,
    VARIABLE_LABEL,
)


def _ic95(serie: pd.Series) -> tuple[float, float]:
    """Intervalo de confianza del 95 % para la media (t de Student).

    Este es el mismo helper usado por el pipeline DCA y heredado aquí.
    """
    n = serie.count()
    if n < 2:
        return (np.nan, np.nan)
    media = serie.mean()
    se = serie.std(ddof=1) / np.sqrt(n)
    inf, sup = st.t.interval(0.95, df=n - 1, loc=media, scale=se)
    return inf, sup


def resumen_descriptivo(df: pd.DataFrame) -> pd.DataFrame:
    """Estadísticos descriptivos por tratamiento + medias + IC95%.

    Para cada valor de ``trt`` (R, T0, T1, T2) calcula:
        • n (observaciones)
        • media, desviación estándar, error estándar
        • IC95 (inferior/superior)
        • mínimo, máximo

    Devuelve un DataFrame con una fila por tratamiento y columna, y guarda como
    ``bdca/resutados/tablas/eda_descriptivos.csv``.
    """
    # Normalizar trt a mayúsculas para consistencia
    df = df.copy()
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()

    filas = []
    for trt in sorted(df["trt"].unique()):
        sub = df[df["trt"] == trt]["yield"].dropna()
        n = int(sub.count())
        if n == 0:
            continue
        media = float(sub.mean())
        desv = float(sub.std(ddof=1)) if n > 1 else np.nan
        err = desv / np.sqrt(n) if n > 1 else np.nan
        inf, sup = _ic95(sub)
        filas.append({
            "trt": trt,
            "n": n,
            "media": round(media, 3),
            "desviacion_estandar": round(desv, 3) if np.isfinite(desv) else np.nan,
            "error_estandar": round(err, 3) if np.isfinite(err) else np.nan,
            "ic95_inferior": round(inf, 3) if np.isfinite(inf) else np.nan,
            "ic95_superior": round(sup, 3) if np.isfinite(sup) else np.nan,
            "minimo": round(float(sub.min()), 3),
            "maximo": round(float(sub.max()), 3),
        })

    tabla = pd.DataFrame(filas)
    guardar_tabla(tabla, "eda_descriptivos", index=False)

    print("=" * 72)
    print("FASE 3.3 - Resumen descriptivo por tratamiento (BDCA)")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print(tabla.to_string(index=False))

    return tabla


def figuras_eda(df: pd.DataFrame) -> list[str]:
    """Genera figuras exploratorias con prefijo ``eda_`` (PNG/PDF).

    Crea:
        1. Boxplot por tratamiento del rendimiento.
        2. Histogramas por tratamiento.
        3. QQ-plot por tratamiento (normalidad).
    Todas las figuras se guardan usando ``save_figure_pub``.
    """
    generadas: list[str] = []
    sns.set_theme(style="whitegrid")

    df = df.copy()
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()

    # 1. Boxplot por tratamiento
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sns.boxplot(data=df, x="trt", y="yield", palette="Set2")
    ax.set_title("Distribución del rendimiento por tratamiento")
    ax.set_xlabel("Tratamiento")
    ax.set_ylabel("Rendimiento")
    save_figure_pub(fig, "eda_boxplot_por_trt", titulo="Boxplot de rendimiento por tratamiento (BDCA)")
    generadas.append("eda_boxplot_por_trt")

    # 2. Histogramas por tratamiento
    fig, ejes = plt.subplots(2, 2, figsize=(16, 12), sharey=True)
    tratamientos = sorted(df["trt"].unique())
    for ax, trt in zip(ejes.ravel(), tratamientos):
        datos = df.loc[df["trt"] == trt, "yield"].dropna()
        ax.hist(datos, bins=12, color="steelblue", edgecolor="white", alpha=0.85)
        ax.set_title(f"Tratamiento {trt}")
        ax.set_xlabel("Rendimiento")
        ax.set_ylabel("Frecuencia")
    fig.suptitle("Histogramas de rendimiento por tratamiento (BDCA)", y=1.02)
    fig.tight_layout()
    save_figure_pub(fig, "eda_histogramas_por_trt", titulo="Histogramas por tratamiento (BDCA)")
    generadas.append("eda_histogramas_por_trt")

    # 3. QQ-plot por tratamiento
    fig, ejes = plt.subplots(2, 2, figsize=(16, 12))
    for ax, trt in zip(ejes.ravel(), tratamientos):
        datos = df.loc[df["trt"] == trt, "yield"].dropna()
        if len(datos) >= 3:
            st.probplot(datos, dist="norm", plot=ax)
            ax.set_title(f"QQ-plot: Tratamiento {trt}")
        else:
            ax.text(0.5, 0.5, f"Datos insuficientes para QQ-plot ({len(datos)})",
                    ha="center", va="center", transform=ax.transAxes)
    fig.suptitle("QQ-plot de normalidad por tratamiento (BDCA)", y=1.02)
    fig.tight_layout()
    save_figure_pub(fig, "eda_qqplot_por_trt", titulo="QQ-plot por tratamiento (BDCA)")
    generadas.append("eda_qqplot_por_trt")

    print(f"  Figuras EDA generadas: {len(generadas)}")
    return generadas


if __name__ == "__main__":
    # Permite ejecución directa para pruebas rápidas
    from pipeline.bdca.cargar import cargar
    df, _ = cargar()
    resumen_descriptivo(df)
    figuras_eda(df)
