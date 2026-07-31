"""Fase 6: verificacion de supuestos de los modelos factoriales.

Para cada variable de respuesta evalúa normalidad de residuos (Shapiro-Wilk),
homocedasticidad (Levene y Bartlett), independencia (Durbin-Watson) y genera
las figuras de diagnostico de residuos. Los resultados se interpretan en
espanol y se guardan como CSV y figuras.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.formula.api as smf
from statsmodels.stats.stattools import durbin_watson
import matplotlib.pyplot as plt

from pipeline.config import (
    METODOS,
    METODO_LABEL,
    PALETA_METODOS,
    VARIABLE_LABEL,
    guardar_tabla,
    save_figure_pub,
)


def _bartlett_seguro(grupos):
    """Bartlett con proteccion ante grupos de varianza nula.

    scipy 1.18 lanza un error interno si algun grupo es constante (varianza 0),
    situacion que ocurre en este dataset (crecimiento 0 = inhibicion completa).
    Se devuelve NaN y una nota en ese caso.
    """
    varianzas = [float(np.var(g)) for g in grupos]
    if any(v == 0 for v in varianzas):
        return float("nan"), float("nan"), "No calculable (grupo con varianza nula)"
    # scipy 1.18 falla con arrays enteros (bug de la API array); se castea a float.
    grupos_float = [np.asarray(g, dtype=float) for g in grupos]
    stat, p = st.bartlett(*grupos_float)
    return float(stat), float(p), "OK"


def verificar_supuestos(df: pd.DataFrame, variable: str) -> dict:
    """Verifica supuestos del modelo OLS factorial para una variable.

    Ajusta ``variable ~ C(metodo_extraccion) * C(aislamiento)`` y evalúa
    normalidad, homocedasticidad e independencia de los residuos, generando
    tabla CSV y figuras de diagnostico.
    """
    modelo = smf.ols(f"{variable} ~ C(metodo_extraccion) * C(aislamiento)", data=df).fit()
    residuos = modelo.resid

    shapiro_w, shapiro_p = st.shapiro(residuos)
    grupos = [df.loc[df["metodo_extraccion"] == m, variable].values for m in METODOS]
    levene_stat, levene_p = st.levene(*grupos)
    bartlett_stat, bartlett_p, bartlett_nota = _bartlett_seguro(grupos)
    dw = durbin_watson(residuos)

    interpretaciones = {
        "shapiro": (
            "No se rechaza normalidad de los residuos (p>0.05)."
            if shapiro_p > 0.05 else
            "Se rechaza normalidad de los residuos (p<0.05)."
        ),
        "levene": (
            "Homocedasticidad aceptada (p>0.05)."
            if levene_p > 0.05 else
            "Heterocedasticidad detectada (p<0.05)."
        ),
        "bartlett": (
            "Homocedasticidad aceptada (p>0.05)."
            if (not np.isnan(bartlett_p)) and bartlett_p > 0.05 else
            ("No calculable (grupo con varianza nula)." if np.isnan(bartlett_p) else
             "Heterocedasticidad detectada (p<0.05).")
        ),
        "dw": (
            "Independencia razonable (DW cercano a 2)."
            if 1.5 <= dw <= 2.5 else
            "Posible autocorrelacion (DW fuera de 1.5-2.5)."
        ),
    }

    tabla = pd.DataFrame([
        {"test": "Shapiro-Wilk (residuos)", "estadistico": round(float(shapiro_w), 4),
         "p_valor": round(float(shapiro_p), 4), "interpretacion": interpretaciones["shapiro"]},
        {"test": "Levene (por metodo)", "estadistico": round(float(levene_stat), 4),
         "p_valor": round(float(levene_p), 4), "interpretacion": interpretaciones["levene"]},
        {"test": "Bartlett (por metodo)", "estadistico": round(bartlett_stat, 4) if not np.isnan(bartlett_stat) else float("nan"),
         "p_valor": round(bartlett_p, 4) if not np.isnan(bartlett_p) else float("nan"),
         "interpretacion": interpretaciones["bartlett"]},
        {"test": "Durbin-Watson (independencia)", "estadistico": round(float(dw), 4),
         "p_valor": float("nan"), "interpretacion": interpretaciones["dw"]},
    ])
    guardar_tabla(tabla, f"supuestos_{variable}", index=False)

    label = VARIABLE_LABEL[variable]
    fig, ejes = plt.subplots(2, 2, figsize=(11, 8))

    ejes[0, 0].hist(residuos, bins=25, color="#4C72B0", edgecolor="white", alpha=0.85)
    ejes[0, 0].set_title("Histograma de residuos")
    ejes[0, 0].set_xlabel("Residuo")

    st.probplot(residuos, dist="norm", plot=ejes[0, 1])
    ejes[0, 1].set_title("QQ-plot de residuos")
    ejes[0, 1].set_ylabel("Cuantiles de los residuos")

    ejes[1, 0].scatter(modelo.fittedvalues, residuos, s=18, alpha=0.6, color="#4C72B0")
    ejes[1, 0].axhline(0, color="grey", lw=0.8)
    ejes[1, 0].set_title("Residuos vs ajustados")
    ejes[1, 0].set_xlabel("Valores ajustados")
    ejes[1, 0].set_ylabel("Residuo")

    for m in METODOS:
        r_m = residuos[df["metodo_extraccion"] == m]
        ejes[1, 1].scatter(
            np.full(len(r_m), METODO_LABEL[m]), r_m, s=18, alpha=0.6,
            color=PALETA_METODOS[m], label=METODO_LABEL[m],
        )
    ejes[1, 1].axhline(0, color="grey", lw=0.8)
    ejes[1, 1].set_title("Residuos por metodo")
    ejes[1, 1].legend(fontsize=8)
    save_figure_pub(
        fig, f"supuestos_{variable}_residuos",
        titulo=f"Diagnostico de supuestos - {label}",
    )

    interpretacion = (
        f"Normalidad: {interpretaciones['shapiro']} | Homocedasticidad (Levene): "
        f"{interpretaciones['levene']} | Independencia (DW={dw:.2f}): {interpretaciones['dw']}"
    )

    print(f"[{label}] Shapiro-Wilk p={shapiro_p:.4f}; Levene p={levene_p:.4f}; "
          f"Bartlett p={'NA' if np.isnan(bartlett_p) else f'{bartlett_p:.4f}'}; "
          f"Durbin-Watson={dw:.2f}")

    return {
        "variable": variable,
        "modelo": modelo,
        "tabla_supuestos": tabla,
        "shapiro_w": float(shapiro_w),
        "shapiro_p": float(shapiro_p),
        "levene_p": float(levene_p),
        "bartlett_p": float(bartlett_p),
        "durbin_watson": float(dw),
        "interpretacion": interpretacion,
    }
