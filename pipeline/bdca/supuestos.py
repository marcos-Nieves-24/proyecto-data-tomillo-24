"""Fase 3.4: Verificación de supuestos para el modelo de bloques RCBD (BDCA).

Evalúa normalidad de residuos (Shapiro-Wilk), homocedasticidad (Levene) e
independencia (Durbin-Watson) sobre los residuos del ANOVA clásico de bloques,
y documenta la ruta elegida para la inferencia principal según AGENTS.md §7.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns

from pipeline.config import (
    guardar_tabla,
    save_figure_pub,
    VARIABLE_LABEL,
)

# Durbin-Watson vive en statsmodels.stats.stattools, no en scipy.stats
from statsmodels.stats.stattools import durbin_watson


def analisis_supuestos(df: pd.DataFrame) -> dict:
    """Evalúa supuestos y decide la ruta de inferencia.

    Realiza:
        1. ANOVA de bloques RCBD clásico (yield ~ C(trt) + C(block)) para obtener residuos.
        2. Shapiro-Wilk (normalidad).
        3. Levene (homocedasticidad entre trt).
        4. Durbin-Watson (independencia serial).

    Documenta violaciones, su magnitud (p-valor) y recomienda:
        • Si los supuestos se cumplen (p > 0.05 para todos): tabla F (paramétrica).
        • Si fallan: ruta no paramétrica (Kruskal-Wallis + Scheirer-Ray-Hare) como inferencia principal,
          conservando la tabla ANOVA como referencia descriptiva.

    Guarda tabla de supuestos, figuras de residuos y una descripción justificada.
    """
    # ---------------------------------------------------------------
    # 1. Modelo de bloques RCBD clásico para residuos
    # ---------------------------------------------------------------
    df = df.copy()
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()
    df["block"] = df["block"].astype(str).str.strip().str.upper()

    # Modelo de bloques: yield ~ C(trt) + C(block).
    # `yield` es palabra reservada de Python: se cita con Q("yield").
    modelo = smf.ols("Q('yield') ~ C(trt) + C(block)", data=df).fit()
    residuos = modelo.resid
    df["residuo"] = residuos

    # ---------------------------------------------------------------
    # 2. Normalidad: Shapiro-Wilk (muestra ≤ 5000)
    # ---------------------------------------------------------------
    if len(residuos) <= 5000:
        shapiro_stat, shapiro_p = st.shapiro(residuos)
    else:
        # Para muestras muy grandes, usar KS o omitir; aquí usamos KS
        ks_stat, shapiro_p = st.kstest(residuos, "norm")
        shapiro_stat = ks_stat

    # ---------------------------------------------------------------
    # 3. Homocedasticidad: Levene entre tratamientos
    # ---------------------------------------------------------------
    grupos = [df.loc[df["trt"] == trt, "yield"].values for trt in sorted(df["trt"].unique())]
    if len(grupos) >= 2:
        levene_stat, levene_p = st.levene(*grupos)
    else:
        levene_stat, levene_p = np.nan, np.nan

    # ---------------------------------------------------------------
    # 4. Independencia serial: Durbin-Watson (ordenado por bloque)
    # ---------------------------------------------------------------
    # Ordena por block para examinar posibles efectos de orden
    df_ordenado = df.sort_values("block")
    residuos_ordenados = df_ordenado["residuo"].values
    dw_stat = durbin_watson(residuos_ordenados) if len(residuos_ordenados) > 1 else np.nan

    # ---------------------------------------------------------------
    # Decisiones de ruta y justificación
    # ---------------------------------------------------------------
    normal_ok = shapiro_p > 0.05 if not np.isnan(shapiro_p) else False
    homo_ok = levene_p > 0.05 if not np.isnan(levene_p) else False
    # Durbin-Watson va de 0 a 4; 2 = sin autocorrelación serial. Se considera
    # independencia aceptable la banda 1.5-2.5 (regla práctica estándar);
    # fuera de ella hay indicio de autocorrelación positiva (<1.5) o negativa (>2.5).
    independiente_ok = not np.isnan(dw_stat) and (1.5 <= dw_stat <= 2.5)

    asuman_az_afirmadas = normal_ok and homo_ok and independiente_ok

    if asuman_az_afirmadas:
        tipo_modelo = "parametrica"
        justificacion = (
            "Los residuos son normales (Shapiro-Wilk p={:.4f}), homocedásticos "
            "(Levene p={:.4f}) y muestran baja autocorrelación serial "
            "(Durbin-Watson = {:.4f}). El ANOVA RCBD como tabla F es adecuado."
        ).format(shapiro_p, levene_p, dw_stat)
    else:
        tipo_modelo = "no_parametrica"
        justificacion = (
            "Uno o más supuestos fallan: normalidad (Shapiro-Wilk p={:.4f}), "
            "homocedasticidad (Levene p={:.4f}) o independencia serial "
            "(Durbin-Watson = {:.4f}). Se usa la ruta no paramétrica "
            "(Kruskal-Wallis por trt y Scheirer-Ray-Hare para interacción) como "
            "inferencia principal, conservando el ANOVA como referencia descriptiva."
        ).format(shapiro_p, levene_p, dw_stat)

    # ---------------------------------------------------------------
    # 5. Tablas y figuras
    # ---------------------------------------------------------------
    # Tabla de supuestos
    tabla_supuestos = pd.DataFrame([
        {"supuesto": "Normalidad (Shapiro-Wilk)", "estadistico": shapiro_stat, "p_valor": shapiro_p, "cumple": normal_ok},
        {"supuesto": "Homocedasticidad (Levene)", "estadistico": levene_stat, "p_valor": levene_p, "cumple": homo_ok},
        {"supuesto": "Independencia serial (Durbin-Watson)", "estadistico": dw_stat, "p_valor": np.nan, "cumple": independiente_ok},
    ])
    guardar_tabla(tabla_supuestos, "supuestos_modelo", index=False)

    # Figura: histograma de residuos con PDF superpuesta
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sns.histplot(residuos, kde=True, stat="density", ax=ax)
    x = np.linspace(residuos.min(), residuos.max(), 100)
    ax.plot(x, st.norm.pdf(x, residuos.mean(), residuos.std(ddof=1)), "r-", label="PDF normal")
    ax.set_title("Histograma de residuos con PDF normal")
    ax.set_xlabel("Residuos")
    ax.set_ylabel("Densidad")
    ax.legend()
    save_figure_pub(fig, "supuestos_residuos_histograma", titulo="Histograma de residuos con PDF normal (BDCA)")

    # Figura: residuos vs ajuste para detectar heterocedasticidad
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ajustado = modelo.fittedvalues
    ax.scatter(ajustado, residuos, alpha=0.7)
    ax.axhline(y=0, color="r", linestyle="-")
    ax.set_xlabel("Valores ajustados")
    ax.set_ylabel("Residuos")
    ax.set_title("Residuos vs valores ajustados")
    save_figure_pub(fig, "supuestos_residuos_vs_ajuste", titulo="Residuos vs valores ajustados (BDCA)")

    print("=" * 72)
    print("FASE 3.4 - Verificación de supuestos (BDCA)")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(tabla_supuestos.to_string(index=False))
    print(f"\nDecisión: {tipo_modelo.upper()} -> {justificacion}")

    return {
        "tipo_modelo": tipo_modelo,
        "tabla_supuestos": tabla_supuestos,
        "justificacion": justificacion,
        "residuos": residuos,
        "modelo_ols": modelo,
    }


if __name__ == "__main__":
    # Permite ejecución directa para pruebas rápidas
    from pipeline.bdca.cargar import cargar
    df, _ = cargar()
    analisis_supuestos(df)
