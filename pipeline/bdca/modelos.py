"""Fase 3.5: Modelos estadísticos para análisis RCBD de rendimiento BDCA.

Implementa:
  • ANOVA clásico de bloques RCBD (yield ~ C(trt) + C(block)) — educativo,
    complemento tabular tipo II (según statsmodels.anova_lm typ=2).
  • Modelo mixto lineal (LMM) como modelo principal: yield ~ C(trt) con bloque aleatorio,
    reportando efectos fijos, varianzas (bloque, residual) e ICC.
  • Nota sobre limitación de aditividad (sin término de interacción porque solo hay una observación por celda trt×block).

Tablas generadas:
  • ANOVA con ANOVA clásico y tamaños de efecto (eta2 parcial).
  • LMM: efectos fijos de tratamiento, varianzas, ICC.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.weightstats import DescrStatsW

from pipeline.config import (
    guardar_tabla,
    VARIABLE_LABEL,
)


def _renombrar_anova_tabla(tabla_anova: pd.DataFrame) -> pd.DataFrame:
    """Renombra índices de tabla ANOVA a nombres legibles para reportes.

    Conversa niveles de factores (C(trt), C(block)) a etiquetas legibles.
    """
    renombres = {
        "C(trt)": "tratamiento",
        "C(block)": "bloque",
    }
    # Renombra índices si están presentes
    if not tabla_anova.index.name:
        tabla_anova.index.name = "fuente"
    # Normaliza nombres de filas
    tabla_anova.index = tabla_anova.index.map(lambda x: renombres.get(x, x))
    tabla_anova = tabla_anova.reset_index().rename(columns={"index": "fuente"})
    return tabla_anova


def _calcular_eta2_parcial(tabla_anova: pd.DataFrame) -> dict:
    """Eta cuadrado parcial = SS_efecto / (SS_efecto + SS_residual).

    Sirve para evaluar magnitud de efecto para cada término.
    """
    # Encuentra SS_residual
    resid_ss = float(tabla_anova.loc[tabla_anova["fuente"] == "Residual", "sum_sq"].iloc[0])
    eta2 = {}
    for _, fila in tabla_anova.iterrows():
        if fila["fuente"] == "Residual":
            continue
        ss = float(fila["sum_sq"])
        eta2[fila["fuente"]] = ss / (ss + resid_ss)
    return eta2


def anova_bloques(df: pd.DataFrame) -> dict:
    """ANOVA clásico de bloques RCBD (complemento educativo).

    Ajusta ``yield ~ C(trt) + C(block)`` y reporta tabla ANOVA tipo II.
    Nota: no hay término de interacción porque solo hay una observación por celda.

    Devuelve un dict con la tabla ANOVA y el modelo ajustado.
    """
    df = df.copy()
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()
    df["block"] = df["block"].astype(str).str.strip().str.upper()

    # `yield` es palabra reservada de Python: patsy no puede evaluarla en la
    # fórmula, por lo que se cita con Q("yield") (misma columna, sin renombrar).
    modelo = smf.ols("Q('yield') ~ C(trt) + C(block)", data=df).fit()
    tabla = anova_lm(modelo, typ=2)  # ANOVA tipo II
    tabla_ren = _renombrar_anova_tabla(tabla)
    eta2 = _calcular_eta2_parcial(tabla_ren)

    # Agregar eta2 parcial a la tabla
    tabla_ren["eta2_parcial"] = tabla_ren["fuente"].map(eta2)
    tabla_ren["eta2_parcial"] = tabla_ren["eta2_parcial"].round(4)

    # Guardar tabla
    guardar_tabla(tabla_ren, "anova_bloques", index=False)

    print("=" * 72)
    print("FASE 3.5 - ANOVA clásico de bloques RCBD (BDCA)")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(tabla_ren.to_string(index=False))

    return {
        "modelo": modelo,
        "tabla_anova": tabla_ren,
        "eta2_parcial": eta2,
    }


def lmm_bloques(df: pd.DataFrame) -> dict:
    """Modelo mixto lineal primario para diseño RCBD BDCA.

    Ajusta ``yield ~ C(trt)`` con bloque aleatorio ``(1|block)`` (REML).
    Reporta:
      • tabla de efectos fijos (coeficientes, error estándar, p-valor, IC95)
      • varianza de bloque, varianza residual y ICC (varianza_bloque / (varianza_bloque + residual))
      • justificación: no se puede probar aditividad porque una observación por celda.

    Devuelve un dict con los resultados del modelo, varianzas e ICC.
    """
    df = df.copy()
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()
    df["block"] = df["block"].astype(str).str.strip().str.upper()

    # LMM: efecto fijo = tratamiento, aleatorio = bloque
    # LMM: efecto fijo = tratamiento, aleatorio = bloque.
    # `yield` es palabra reservada de Python: se cita con Q("yield").
    # `reml` es parámetro de .fit(), no del constructor.
    modelo = smf.mixedlm("Q('yield') ~ C(trt)", groups=df["block"], data=df)
    ajuste = modelo.fit(reml=True)

    # Tabla de efectos fijos
    filas_fijos = []
    for indice in ajuste.params.index:
        coef = float(ajuste.params[indice])
        se = float(ajuste.bse[indice])
        t = float(ajuste.tvalues[indice])
        p = float(ajuste.pvalues[indice])
        filas_fijos.append({
            "efecto": indice,
            "coeficiente": round(coef, 4),
            "error_estandar": round(se, 4),
            "t": round(t, 4),
            "p_valor": round(p, 4),
            "ic95_inferior": round(coef - 1.96 * se, 4),
            "ic95_superior": round(coef + 1.96 * se, 4),
        })
    tabla_fija = pd.DataFrame(filas_fijos)

    # Varianzas: bloque (cov_re) y residual (scale)
    var_bloque = float(ajuste.cov_re.iloc[0, 0]) if ajuste.cov_re.shape[0] > 0 else 0.0
    var_residual = float(ajuste.scale)
    icc = var_bloque / (var_bloque + var_residual) if (var_bloque + var_residual) > 0 else 0.0

    # Guardar tablas
    guardar_tabla(tabla_fija, "lmm_bloques_fijos", index=False)
    guardar_tabla(pd.DataFrame([
        {"parametro": "varianza_bloque", "valor": round(var_bloque, 4)},
        {"parametro": "varianza_residual", "valor": round(var_residual, 4)},
        {"parametro": "ICC", "valor": round(icc, 4)},
    ]), "lmm_bloques_varianzas", index=False)

    print("=" * 72)
    print("FASE 3.5 - Modelo mixto lineal de bloques RCBD (BDCA)")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(tabla_fija.to_string(index=False))
    print(f"\nVarianza de bloque = {var_bloque:.4f}")
    print(f"Varianza residual = {var_residual:.4f}")
    print(f"ICC = {icc:.4f}")

    # Limite de aditividad
    limitacion_aditividad = (
        "Dado un solo cultivo por celda trt × bloque, el término de interacción "
        "no puede ser estimado. Por lo tanto, la aditividad (efecto aditivo puro) "
        "es una suposición no testable; la inferencia se basa en el modelo de bloques "
        "RCBD sin interacción."
    )

    return {
        "modelo": ajuste,
        "tabla_fija": tabla_fija,
        "var_bloque": var_bloque,
        "var_residual": var_residual,
        "icc": icc,
        "limitacion_aditividad": limitacion_aditividad,
    }


if __name__ == "__main__":
    # Permite ejecución directa para pruebas rápidas
    from pipeline.bdca.cargar import cargar
    df, _ = cargar()

    print("=== ANOVA clásico ===")
    anova_bloques(df)
    print("\n=== Modelo mixto lineal ===")
    lmm_bloques(df)
