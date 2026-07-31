"""Fase 7: analisis estadisticos principales con seleccion automatica.

Evalúa el rendimiento de extraccion (ANOVA de una via), el efecto factorial
metodo x aislado sobre cada variable de respuesta (con seleccion automatica
de la ruta parametrica/no parametrica/GLM), analisis de sensibilidad con
modelos mixtos lineales (LMM con aislado aleatorio) y el diagnostico
especifico de conidias (por que no aplica Poisson).

Filosofia del proyecto: no se elige un test solo por producir significancia;
la seleccion se justifica con supuestos y diagnostico de datos.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

from pipeline.config import (
    METODOS,
    VARIABLE_LABEL,
    guardar_tabla,
)

# Mediciones continuas que, aunque se registren como enteros, no son conteos.
VARIABLES_MEDICION_CONTINUA = {"crecimiento_micelial_mm"}

NOMBRE_FACTOR_METODO = "C(metodo_extraccion)"
NOMBRE_FACTOR_AISLADO = "C(aislamiento)"
NOMBRE_INTERACCION = "C(metodo_extraccion):C(aislamiento)"


# ---------------------------------------------------------------------------
# Utilidades compartidas
# ---------------------------------------------------------------------------


def _tipo_dato(serie: pd.Series, variable: str) -> str:
    """Clasifica la variable como 'conteo' o 'continua' con justificacion.

    Un conteo exige valores enteros no negativos, mas de 5 valores distintos
    y sobredispersion (varianza > media). Las mediciones en mm se excluyen.
    """
    s = serie.dropna()
    es_entero = bool(np.all(s == np.floor(s)))
    no_negativo = bool((s >= 0).all())
    if variable in VARIABLES_MEDICION_CONTINUA:
        return "continua"
    if es_entero and no_negativo and s.nunique() > 5 and s.var() > s.mean():
        return "conteo"
    return "continua"


def _renombrar_anova(tabla: pd.DataFrame) -> pd.DataFrame:
    """Renombra los indices de la tabla ANOVA a nombres legibles."""
    renombres = {
        NOMBRE_FACTOR_METODO: "metodo_extraccion",
        NOMBRE_FACTOR_AISLADO: "aislamiento",
        NOMBRE_INTERACCION: "metodo_extraccion:aislamiento",
    }
    tabla = tabla.rename(index=renombres)
    tabla = tabla.reset_index().rename(columns={"index": "fuente"})
    return tabla


def _eta_cuadrado_parcial(tabla: pd.DataFrame, residuo_ss: float) -> dict:
    """Eta cuadrado parcial = SS_efecto / (SS_efecto + SS_residual)."""
    out = {}
    for _, fila in tabla.iterrows():
        if fila["fuente"] == "Residual":
            continue
        ss = float(fila["sum_sq"])
        out[fila["fuente"]] = ss / (ss + residuo_ss)
    return out


def _omega_cuadrado(tabla: pd.DataFrame) -> dict:
    """Omega cuadrado parcial por termino."""
    out = {}
    residuo_ss = float(tabla.loc[tabla["fuente"] == "Residual", "sum_sq"].iloc[0])
    residuo_df = float(tabla.loc[tabla["fuente"] == "Residual", "df"].iloc[0])
    ms_resid = residuo_ss / residuo_df if residuo_df else np.nan
    for _, fila in tabla.iterrows():
        if fila["fuente"] == "Residual":
            continue
        ss = float(fila["sum_sq"])
        df = float(fila["df"])
        ms = ss / df if df else np.nan
        out[fila["fuente"]] = (df * (ms - ms_resid)) / (ss + residuo_ss + ms_resid)
    return out


def _tabla_anotada(tabla: pd.DataFrame, eta_p: dict, omega_p: dict) -> pd.DataFrame:
    """Anade columnas de tamano de efecto a la tabla ANOVA."""
    out = tabla.copy()
    out["eta2_parcial"] = out["fuente"].map(eta_p)
    out["omega2_parcial"] = out["fuente"].map(omega_p)
    out["eta2_parcial"] = out["eta2_parcial"].round(4)
    out["omega2_parcial"] = out["omega2_parcial"].round(4)
    return out


# ---------------------------------------------------------------------------
# Objetivo 1: rendimiento de extraccion
# ---------------------------------------------------------------------------


def analisis_rendimiento(df_rend: pd.DataFrame) -> dict:
    """ANOVA de una via del rendimiento de extraccion segun metodo.

    Reporta tabla ANOVA tipo II, eta cuadrado, omega cuadrado y supuestos.
    Si los supuestos fallan, complementa con Kruskal-Wallis.
    """
    df = df_rend.copy()
    modelo = smf.ols("rendimiento_pct ~ C(metodo_extraccion)", data=df).fit()
    tabla = anova_lm(modelo, typ=2)

    ss_total = float(tabla["sum_sq"].sum())
    ss_efecto = float(tabla.loc[NOMBRE_FACTOR_METODO, "sum_sq"])
    ss_resid = float(tabla.loc["Residual", "sum_sq"])
    df_efecto = float(tabla.loc[NOMBRE_FACTOR_METODO, "df"])
    ms_error = ss_resid / float(tabla.loc["Residual", "df"])
    eta2 = ss_efecto / ss_total
    omega2 = (ss_efecto - df_efecto * ms_error) / (ss_total + ms_error)

    shapiro_w, shapiro_p = st.shapiro(modelo.resid)
    grupos = [df.loc[df["metodo_extraccion"] == m, "rendimiento_pct"].values for m in METODOS]
    levene_stat, levene_p = st.levene(*grupos)
    kw_stat, kw_p = st.kruskal(*grupos)

    tabla_anotada = _tabla_anotada(
        _renombrar_anova(tabla),
        _eta_cuadrado_parcial(_renombrar_anova(tabla), ss_resid),
        _omega_cuadrado(_renombrar_anova(tabla)),
    )
    guardar_tabla(tabla_anotada, "modelos_rendimiento", index=False)

    guardar_tabla(pd.DataFrame([
        {"estadistico": "Shapiro-Wilk (residuos)", "valor": round(float(shapiro_w), 4),
         "p_valor": round(float(shapiro_p), 4)},
        {"estadistico": "Levene", "valor": round(float(levene_stat), 4),
         "p_valor": round(float(levene_p), 4)},
        {"estadistico": "Kruskal-Wallis", "valor": round(float(kw_stat), 4),
         "p_valor": round(float(kw_p), 4)},
    ]), "supuestos_rendimiento", index=False)

    tipo_modelo = "anova_ols"
    justificacion = (
        "El ANOVA de una via es adecuado si los residuos son normales y "
        "homocedasticos. "
    )
    if shapiro_p > 0.05 and levene_p > 0.05:
        justificacion += "Ambos supuestos se cumplen; se usa la tabla F como inferencia principal."
    else:
        justificacion += (
            f"Supuestos no plenamente cumplidos (Shapiro p={shapiro_p:.3f}, "
            f"Levene p={levene_p:.3f}); se complementa con Kruskal-Wallis."
        )
    if kw_p < 0.05:
        justificacion += " Kruskal-Wallis confirma diferencias entre metodos."
    else:
        justificacion += " Kruskal-Wallis no confirma diferencias."

    medias = df.groupby("metodo_extraccion")["rendimiento_pct"].mean().round(3)

    print("=" * 72)
    print("FASE 7 - Rendimiento de extracción (ANOVA de una vía)")
    print("=" * 72)
    print("Medias de rendimiento (%):")
    print(medias.to_string())
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print("\nTabla ANOVA:")
        print(tabla_anotada.to_string(index=False))
    print(f"eta2={eta2:.4f} | omega2={omega2:.4f} | Kruskal-Wallis p={kw_p:.4f}")

    return {
        "tipo_modelo": tipo_modelo,
        "tabla_anova": tabla_anotada,
        "eta2": float(eta2),
        "omega2": float(omega2),
        "shapiro_p": float(shapiro_p),
        "levene_p": float(levene_p),
        "kruskal_wallis": {"H": float(kw_stat), "p_valor": float(kw_p)},
        "medias": medias,
        "justificacion": justificacion,
    }


# ---------------------------------------------------------------------------
# Rama GLM para conteos (condicional; inactiva en este dataset)
# ---------------------------------------------------------------------------


def glm_conteos(df_bio: pd.DataFrame, variable: str) -> dict:
    """Rama condicional GLM Poisson/NB para variables de conteo sobredispersas.

    Solo se activa si la variable es un conteo entero >=0 con sobredispersion.
    Para este dataset ninguna variable la activa (las conidias son log10
    continuas, el %INH conidias incluye negativos y el crecimiento en mm es una
    medicion). La funcion queda documentada para datasets futuros.
    """
    modelo_p = smf.glm(
        f"{variable} ~ C(metodo_extraccion) * C(aislamiento)",
        data=df_bio, family=sm.families.Poisson(),
    ).fit()
    dispersion = float(modelo_p.pearson_chi2 / modelo_p.df_resid)
    familia = "Poisson"
    modelo = modelo_p
    if dispersion > 1.5:
        modelo = smf.glm(
            f"{variable} ~ C(metodo_extraccion) * C(aislamiento)",
            data=df_bio, family=sm.families.NegativeBinomial(),
        ).fit()
        familia = "BinomialNegativa"

    print(f"[{VARIABLE_LABEL[variable]}] Rama GLM {familia} activada "
          f"(dispersion={dispersion:.2f}).")

    return {
        "variable": variable,
        "tipo_modelo": "glm_conteos",
        "familia": familia,
        "dispersion": dispersion,
        "modelo": modelo,
        "justificacion": (
            f"Conteo entero >=0 sobredisperso (dispersion {dispersion:.2f}); "
            f"se usa GLM {familia}."
        ),
    }


# ---------------------------------------------------------------------------
# Ruta no parametrica
# ---------------------------------------------------------------------------


def _kruskal_por_metodo(df: pd.DataFrame, variable: str) -> dict:
    """Kruskal-Wallis de la variable entre metodos."""
    grupos = [df.loc[df["metodo_extraccion"] == m, variable].dropna().values for m in METODOS]
    h, p = st.kruskal(*grupos)
    return {"H": float(h), "p_valor": float(p), "df": len(grupos) - 1}


def _scheirer_ray_hare(df: pd.DataFrame, variable: str) -> dict:
    """Test no parametrico de dos vias (Scheirer-Ray-Hare) sobre rangos.

    Se calcula como ANOVA tipo II sobre los rangos promedio; el estadistico
    H = SS_factor / (SS_total/(N-1)) sigue una chi-cuadrado.
    """
    d = df.copy()
    d["rango"] = d[variable].rank(method="average")
    modelo = smf.ols("rango ~ C(metodo_extraccion) * C(aislamiento)", data=d).fit()
    anova = anova_lm(modelo, typ=2)
    n = len(d)
    ss_total = float(((d["rango"] - d["rango"].mean()) ** 2).sum())
    ms_total = ss_total / (n - 1)

    resultado = {}
    for termino in (NOMBRE_FACTOR_METODO, NOMBRE_FACTOR_AISLADO, NOMBRE_INTERACCION):
        ss = float(anova.loc[termino, "sum_sq"])
        dfi = int(anova.loc[termino, "df"])
        h = ss / ms_total
        p = st.chi2.sf(h, dfi)
        resultado[termino] = {"H": float(h), "df": dfi, "p_valor": float(p)}

    tabla = pd.DataFrame([
        {"fuente": _renombrar_fuente(k), "H": round(v["H"], 4),
         "df": v["df"], "p_valor": round(v["p_valor"], 4)}
        for k, v in resultado.items()
    ])
    guardar_tabla(tabla, f"no_parametrico_{variable}", index=False)
    return resultado


def _renombrar_fuente(termino: str) -> str:
    renombres = {
        NOMBRE_FACTOR_METODO: "metodo_extraccion",
        NOMBRE_FACTOR_AISLADO: "aislamiento",
        NOMBRE_INTERACCION: "metodo_extraccion:aislamiento",
    }
    return renombres.get(termino, termino)


# ---------------------------------------------------------------------------
# Objetivos 2 y 3: analisis factorial
# ---------------------------------------------------------------------------


def analisis_factorial(df_bio: pd.DataFrame, variable: str) -> dict:
    """Ajusta el modelo factorial y selecciona la ruta de inferencia.

    Si la variable es un conteo entero >=0 sobredisperso usa la rama GLM
    Poisson/NB. Si es continua ajusta el OLS factorial
    ``variable ~ metodo * aislado`` y, segun los supuestos, decide entre la
    tabla F (parametrica) o la ruta no parametrica (Kruskal-Wallis +
    Scheirer-Ray-Hare). Devuelve tipo_modelo, resultados y justificacion.
    """
    tipo = _tipo_dato(df_bio[variable], variable)
    if tipo == "conteo":
        return glm_conteos(df_bio, variable)

    df = df_bio.copy()
    modelo = smf.ols(f"{variable} ~ C(metodo_extraccion) * C(aislamiento)", data=df).fit()
    tabla = _renombrar_anova(anova_lm(modelo, typ=2))
    ss_resid = float(tabla.loc[tabla["fuente"] == "Residual", "sum_sq"].iloc[0])
    eta_p = _eta_cuadrado_parcial(tabla, ss_resid)
    omega_p = _omega_cuadrado(tabla)

    shapiro_w, shapiro_p = st.shapiro(modelo.resid)
    grupos = [df.loc[df["metodo_extraccion"] == m, variable].values for m in METODOS]
    levene_stat, levene_p = st.levene(*grupos)

    resultados: dict = {
        "variable": variable,
        "tipo_modelo": "factorial_ols",
        "modelo_ols": modelo,
        "tabla_anova": _tabla_anotada(tabla, eta_p, omega_p),
        "eta_cuadrado_parcial": eta_p,
        "omega_cuadrado_parcial": omega_p,
        "shapiro_w": float(shapiro_w),
        "shapiro_p": float(shapiro_p),
        "levene_p": float(levene_p),
    }
    guardar_tabla(resultados["tabla_anova"], f"modelos_{variable}", index=False)

    cumplen = (shapiro_p > 0.05) and (levene_p > 0.05)
    if not cumplen:
        resultados["tipo_modelo"] = "factorial_no_parametrico"
        resultados["kruskal_wallis"] = _kruskal_por_metodo(df, variable)
        resultados["scheirer_ray_hare"] = _scheirer_ray_hare(df, variable)
        resultados["justificacion"] = (
            f"Los supuestos del ANOVA factorial no se cumplen "
            f"(Shapiro-Wilk p={shapiro_p:.4f}; Levene p={levene_p:.4f}). "
            "Se utiliza la via no parametrica (Kruskal-Wallis por metodo y "
            "Scheirer-Ray-Hare para la interaccion) como inferencia principal, "
            "conservando la tabla ANOVA y los tamanos de efecto como referencia "
            "descriptiva. Para %INH micelial, el efecto techo (muchos valores "
            "= 100) explica la violacion de normalidad."
        )
    else:
        resultados["justificacion"] = (
            "Los supuestos de normalidad y homocedasticidad se cumplen "
            "(Shapiro-Wilk y Levene p>0.05); el ANOVA factorial es la via "
            "adecuada. Nota: las replicas de %INH comparten el control C4 "
            "(pseudorreplicacion), lo que puede inflar levemente la precision."
        )

    print("=" * 72)
    print(f"FASE 7 - Análisis factorial: {VARIABLE_LABEL[variable]}")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(resultados["tabla_anova"].to_string(index=False))
    print(f"Tipo de modelo: {resultados['tipo_modelo']}")
    print(f"Justificacion: {resultados['justificacion']}")

    return resultados


# ---------------------------------------------------------------------------
# Analisis de sensibilidad con LMM
# ---------------------------------------------------------------------------


def analisis_sensibilidad_lmm(df_bio: pd.DataFrame, variable: str) -> dict:
    """Modelo mixto lineal con aislado aleatorio como analisis de sensibilidad.

    Ajusta ``variable ~ C(metodo_extraccion) + (1|aislamiento)``, reporta
    efectos fijos, ICC (varianza del aislado / varianza total) y compara la
    conclusion con el modelo factorial.
    """
    df = df_bio.copy()
    modelo = smf.mixedlm(
        f"{variable} ~ C(metodo_extraccion)", df, groups=df["aislamiento"]
    )
    ajuste = modelo.fit(reml=True)

    filas = []
    for indice in ajuste.params.index:
        coef = float(ajuste.params[indice])
        se = float(ajuste.bse[indice])
        t = float(ajuste.tvalues[indice])
        p = float(ajuste.pvalues[indice])
        filas.append({
            "efecto": str(indice),
            "coeficiente": round(coef, 4),
            "error_estandar": round(se, 4),
            "t": round(t, 4),
            "p_valor": round(p, 4),
            "ic95_inferior": round(coef - 1.96 * se, 4),
            "ic95_superior": round(coef + 1.96 * se, 4),
        })
    tabla_fija = pd.DataFrame(filas)

    var_aislado = float(ajuste.cov_re.iloc[0, 0])
    var_residual = float(ajuste.scale)
    icc = var_aislado / (var_aislado + var_residual)

    # p del efecto global del metodo con prueba de razon de verosimilitud
    modelo_reducido = smf.mixedlm(f"{variable} ~ 1", df, groups=df["aislamiento"])
    ajuste_reducido = modelo_reducido.fit(reml=False)
    lrt = 2 * (ajuste.llf - ajuste_reducido.llf)
    p_metodo_lrt = float(st.chi2.sf(lrt, len(METODOS) - 1))

    guardar_tabla(tabla_fija, f"lmm_{variable}", index=False)
    guardar_tabla(pd.DataFrame([
        {"parametro": "varianza_aislado", "valor": round(var_aislado, 4)},
        {"parametro": "varianza_residual", "valor": round(var_residual, 4)},
        {"parametro": "ICC", "valor": round(icc, 4)},
        {"parametro": "LRT_metodo", "valor": round(lrt, 4)},
        {"parametro": "p_valor_metodo_LRT", "valor": round(p_metodo_lrt, 4)},
    ]), f"lmm_{variable}_varianzas", index=False)

    # Significancia del metodo en el modelo mixto
    print("=" * 72)
    print(f"FASE 7 - LMM sensibilidad: {VARIABLE_LABEL[variable]}")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(tabla_fija.to_string(index=False))
    print(f"ICC (aislado)={icc:.4f}; p_metodo (LRT)={p_metodo_lrt:.4f}")

    return {
        "variable": variable,
        "modelo": ajuste,
        "tabla_efectos_fijos": tabla_fija,
        "icc": icc,
        "var_aislado": var_aislado,
        "var_residual": var_residual,
        "lrt": float(lrt),
        "p_valor_metodo_lrt": p_metodo_lrt,
    }


# ---------------------------------------------------------------------------
# Diagnostico de conidias (Objetivo 3)
# ---------------------------------------------------------------------------


def analisis_conidias(df_bio: pd.DataFrame) -> dict:
    """Diagnostico de la variable de conidias y modelo lineal sobre log10.

    Muestra que conidias_log10_ml es continua (no entera) y que por lo tanto
    no aplica la rama Poisson/NB; luego ajusta el modelo factorial lineal
    sobre la escala log10.
    """
    s = df_bio["conidias_log10_ml"].dropna()
    pct_enteros = float((s == np.floor(s)).mean() * 100)
    diag = pd.DataFrame({
        "metrica": ["n", "media", "varianza", "minimo", "maximo", "pct_valores_enteros"],
        "valor": [
            int(s.count()),
            round(float(s.mean()), 4),
            round(float(s.var(ddof=1)), 4),
            float(s.min()),
            float(s.max()),
            round(pct_enteros, 2),
        ],
    })
    guardar_tabla(diag, "conidias_diagnostico", index=False)

    explicacion = (
        "La variable conidias_log10_ml es una concentracion de conidias ya "
        "transformada en log10 por el laboratorio (valores continuos, no "
        "enteros). No se trata de conteos crudos, por lo que la rama "
        "Poisson/Binomial negativa del pipeline NO aplica: dicha rama solo se "
        "activa si se detectaran conteos enteros >=0 sobredispersos. Ademas, "
        "el %INH conidias esta expresado en escala log10. Por ello se utiliza "
        "un modelo lineal sobre la escala log10, que es la escala natural de "
        "medicion reportada por el laboratorio."
    )

    factorial = analisis_factorial(df_bio, "conidias_log10_ml")

    print("=" * 72)
    print("FASE 7 - Diagnóstico de conidias")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print(diag.to_string(index=False))
    print(f"\n{explicacion}")

    return {"diagnostico": diag, "explicacion": explicacion, "factorial": factorial}
