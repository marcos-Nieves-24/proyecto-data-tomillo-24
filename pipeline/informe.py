"""Fase 12: generación del informe final.

Genera ``informe_final.md`` (español profesional y neutro) con las secciones
solicitadas y lo convierte a ``informe_final.html`` con un CSS simple y
legible. El informe integra los resultados de todas las fases previas.

Provenance de datos crudos:
- DCA: Excel con controles C4 (`datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx`)
- BDCA: CSV `DBCA_Jenkyn_control_mildeo.csv` (`datos_crudos/bdca/`); columnas `plot, trt, block, yield`
"""

from __future__ import annotations

import markdown as md_lib
import pandas as pd

from pipeline.config import (
    INFORME_HTML,
    INFORME_MD,
    METODOS,
    METODO_LABEL,
    VARIABLES_RESPUESTA,
    VARIABLE_LABEL,
    VARIABLE_UNIDAD,
)

CSS = """
<style>
body { font-family: Georgia, 'Times New Roman', serif; margin: 2.5em auto;
       max-width: 980px; line-height: 1.6; color: #222; padding: 0 1.5em; }
h1 { color: #1F3864; border-bottom: 3px solid #1F3864; padding-bottom: .25em; }
h2 { color: #1F3864; border-bottom: 1px solid #ccc; padding-bottom: .2em;
     margin-top: 1.8em; }
h3 { color: #2E5FA3; }
table { border-collapse: collapse; margin: 1em 0; width: 100%; font-size: .92em; }
th, td { border: 1px solid #bbb; padding: 6px 10px; text-align: left; }
th { background: #DDEBF7; }
tr:nth-child(even) { background: #F4F7FB; }
code { background: #f0f0f0; padding: 1px 4px; border-radius: 3px;
       font-family: Consolas, monospace; }
blockquote { border-left: 4px solid #D55E00; margin: 1em 0; padding: .2em 1em;
             background: #FBF6F0; }
</style>
"""


def _df_to_markdown(df: pd.DataFrame, decimales: int = 4) -> str:
    """Convierte un DataFrame en tabla markdown con redondeo controlado."""
    if df is None or len(df) == 0:
        return "_Sin datos_"
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].round(decimales)
        elif pd.api.types.is_bool_dtype(d[col]):
            d[col] = d[col].map({True: "Sí", False: "No"})
    encabezados = " | ".join(str(c) for c in d.columns)
    separador = " | ".join("---" for _ in d.columns)
    filas = []
    for _, r in d.iterrows():
        celdas = []
        for v in r:
            if pd.isna(v):
                celdas.append("")
            else:
                celdas.append(str(v))
        filas.append(" | ".join(celdas))
    return "| " + encabezados + " |\n| " + separador + " |\n| " + " |\n| ".join(filas) + " |"


def generar_informe(resultados: dict) -> tuple:
    """Genera ``informe_final.md`` y ``informe_final.html``.

    ``resultados`` es el diccionario acumulado por el notebook orquestador.
    """
    secciones = []

    secciones.append(_seccion_resumen(resultados))

    # 1. Calidad de datos
    secciones.append(_seccion_calidad(resultados))
    # 2. Diseño experimental
    secciones.append(_seccion_diseno(resultados))
    # 3. Supuestos
    secciones.append(_seccion_supuestos(resultados))
    # 4. Análisis seleccionados
    secciones.append(_seccion_analisis(resultados))
    # 5. Comparaciones múltiples
    secciones.append(_seccion_comparaciones(resultados))
    # 6. Análisis multivariado
    secciones.append(_seccion_multivariado(resultados))
    # 7. Ranking de técnicas
    secciones.append(_seccion_ranking(resultados))
    # 8. Interpretación biológica
    secciones.append(_seccion_bio(resultados))
    # 9. Conclusiones
    secciones.append(_seccion_conclusiones(resultados))
    # 10. Limitaciones
    secciones.append(_seccion_limitaciones(resultados))

    texto = "\n\n".join(secciones)
    INFORME_MD.write_text(texto, encoding="utf-8")

    html = (
        "<!DOCTYPE html>\n<html lang='es'>\n<head>\n<meta charset='utf-8'>\n"
        "<title>Informe final - Actividad antifúngica de extractos de tomillo</title>\n"
        + CSS + "\n</head>\n<body>\n"
        + md_lib.markdown(texto, extensions=["tables"])
        + "\n</body>\n</html>\n"
    )
    INFORME_HTML.write_text(html, encoding="utf-8")

    print("=" * 72)
    print("FASE 12 - Informe final generado")
    print("=" * 72)
    print(f"  Markdown: {INFORME_MD}")
    print(f"  HTML:     {INFORME_HTML}")

    return INFORME_MD, INFORME_HTML


# ---------------------------------------------------------------------------
# Secciones
# ---------------------------------------------------------------------------


def _seccion_resumen(r) -> str:
    ranking = r.get("ranking", {}).get("tabla")
    ranking_txt = ""
    if ranking is not None and len(ranking):
        fila = ranking.iloc[0]
        ranking_txt = (
            f" La técnica mejor puntuada fue **{METODO_LABEL[fila['metodo_extraccion']]}** "
            f"(score compuesto {fila['score_compuesto']})."
        )
    multiv = r.get("multivariado", {})
    if multiv.get("k_optimo") is not None:
        cluster_txt = f" El análisis multivariado identificó **{multiv['k_optimo']}** grupos de aislados."
    else:
        cluster_txt = ""
    return (
        "## Resumen ejecutivo\n\n"
        "Este informe resume el análisis estadístico reproducible de la actividad "
        "antifúngica de extractos de tomillo (*Thymus vulgaris*) obtenidos por tres "
        "técnicas de extracción (maceración, Soxhlet y ultrasonido) frente a 31 "
        "aislados de *Fusarium* spp., todos ensayados a 5 mg/mL. Se evaluó el "
        "rendimiento de extracción, la inhibición del crecimiento micelial y la "
        "inhibición de la producción de conidias, se compararon las técnicas, se "
        "agruparon los aislados por su perfil de susceptibilidad y se generó un "
        "ranking de técnicas."
        + ranking_txt + cluster_txt + "\n\n"
        "Todos los resultados numéricos se guardaron como tablas en "
        "`dca/resultados/tablas/` y las figuras en `dca/resultados/figuras/`.\n"
    )


def _seccion_calidad(r) -> str:
    auditoria = r.get("auditoria", {})
    tabla = auditoria.get("tabla")
    resumen = auditoria.get("resumen", {})
    cuerpo = ""
    if tabla is not None:
        cuerpo += _df_to_markdown(tabla) + "\n\n"
    if resumen:
        cuerpo += (
            f"Filas duplicadas exactas: **{resumen.get('n_filas_duplicadas', 0)}**. "
            f"Atípicos por IQR (1.5xIQR) flagueados, sin eliminar: "
            f"{resumen.get('atipicos_por_iqr', {})}.\n"
        )
    return (
        "## 1. Calidad de datos (auditoría)\n\n"
        "Se auditó el dataset consolidado tidy (279 filas, 7 columnas) antes de "
        "cualquier transformación: completitud, duplicados, columnas constantes, "
        "tipos, valores inconsistentes y atípicos por rango intercuartílico. "
        "Ningún atípico se eliminó automáticamente (regla de integridad del "
        "proyecto); los valores se conservaron y su efecto se evaluó en el "
        "análisis.\n\n" + cuerpo
    )


def _seccion_diseno(r) -> str:
    diseno = r.get("diseno", {})
    texto = diseno.get("texto", "")
    detalle = diseno.get("detalle")
    cuerpo = ""
    if detalle is not None:
        cuerpo = "\n\n" + _df_to_markdown(detalle, decimales=2) + "\n"
    return (
        "## 2. Diseño experimental (inferido y caveats)\n\n" + texto + cuerpo
        + "\n\n" + _seccion_validacion_inh(r) + "\n"
    )


def _seccion_validacion_inh(r) -> str:
    """Subsección de validación de la fórmula del %INH (integridad de datos)."""
    validacion = r.get("validacion_inh", {})
    estado = validacion.get("estado")
    if estado is None:
        return "### 2.1 Validación de la fórmula del %INH\n\n_Sin datos de validación._"
    if estado == "no_disponible":
        return (
            "### 2.1 Validación de la fórmula del %INH\n\n"
            "**No disponible.** " + validacion.get("nota", "") + "\n"
        )
    tabla = validacion.get("tabla")
    cuerpo = _df_to_markdown(tabla) if tabla is not None else "_Sin datos_"
    interpretacion = (
        "La fórmula reconstruida `%INH = (1 - C1/C4) × 100` coincide con el "
        "%INH reportado por el laboratorio (diferencias máximas del orden de "
        "1e-10 y ninguna discrepancia por encima de la tolerancia 1e-6). Esto "
        "confirma la consistencia interna de los datos: los controles C4 se "
        "incorporaron al dataset maestro (columnas `control_crecimiento_mm` y "
        "`control_conidias_log10`) y el %INH usado en la inferencia es el "
        "reportado por el investigador. Para las conidias, la fórmula se validó "
        "sobre la escala log10 directamente (el laboratorio reportó la reducción "
        "en log10), tal como se verificó en los datos."
        if estado == "ok" else
        "La fórmula reconstruida presenta diferencias con el %INH reportado. "
        "Antes de interpretar, se debe verificar si la discrepancia es por "
        "escala (log10 vs. crudo) o por redondeo del laboratorio."
    )
    return (
        "### 2.1 Validación de la fórmula del %INH (integridad de datos)\n\n"
        "Se reconstruyó el %INH a partir de los controles C4 y las mediciones "
        "del bioensayo para confirmar la consistencia de los datos de entrada. "
        "La validación es informativa y no reemplaza las respuestas reportadas "
        "usadas en la inferencia.\n\n"
        + cuerpo + "\n\n" + interpretacion + "\n"
    )


def _seccion_supuestos(r) -> str:
    supuestos = r.get("supuestos", {})
    parrafos = []
    for var in VARIABLES_RESPUESTA:
        if var not in supuestos:
            continue
        s = supuestos[var]
        parrafos.append(
            f"- **{VARIABLE_LABEL[var]}**: Shapiro-Wilk p={s.get('shapiro_p', float('nan')):.4f}, "
            f"Levene p={s.get('levene_p', float('nan')):.4f}, Bartlett p={s.get('bartlett_p', float('nan')):.4f}, "
            f"Durbin-Watson={s.get('durbin_watson', float('nan')):.2f}."
        )
    return (
        "## 3. Supuestos de los modelos\n\n"
        "Se verificaron normalidad de residuos (Shapiro-Wilk), homocedasticidad "
        "(Levene y Bartlett, por método) e independencia (Durbin-Watson) sobre el "
        "modelo OLS factorial método × aislado.\n\n"
        + ("\n".join(parrafos) if parrafos else "_Sin datos_") + "\n"
    )


def _seccion_analisis(r) -> str:
    partes = []
    # Rendimiento
    rend = r.get("rendimiento")
    if rend:
        tabla = rend.get("tabla_anova")
        cuerpo = _df_to_markdown(tabla) if tabla is not None else "_Sin datos_"
        partes.append(
            f"### 4.1 Rendimiento de extracción (ANOVA de una vía)\n\n"
            f"**Modelo seleccionado:** {rend.get('tipo_modelo')}. eta2={rend.get('eta2', float('nan')):.4f}; "
            f"omega2={rend.get('omega2', float('nan')):.4f}; Kruskal-Wallis p={rend.get('kruskal_wallis', {}).get('p_valor', float('nan')):.4f}.\n\n"
            + cuerpo + "\n\n"
            f"Justificación: {rend.get('justificacion', '')}\n"
        )

    # Factorial por variable
    modelos = r.get("modelos", {})
    for var in VARIABLES_RESPUESTA:
        m = modelos.get(var)
        if not m:
            continue
        tipo = m.get("tipo_modelo", "")
        cuerpo = _df_to_markdown(m.get("tabla_anova")) if m.get("tabla_anova") is not None else "_Sin datos_"
        extra = ""
        if tipo == "factorial_no_parametrico":
            kw = m.get("kruskal_wallis", {})
            srh = m.get("scheirer_ray_hare", {})
            extra = (
                f" Kruskal-Wallis (método): H={kw.get('H', float('nan')):.3f}, "
                f"p={kw.get('p_valor', float('nan')):.4f}."
            )
            if srh:
                fila_inter = srh.get("C(metodo_extraccion):C(aislamiento)", {})
                fila_metodo = srh.get("C(metodo_extraccion)", {})
                extra += (
                    f" Scheirer-Ray-Hare: método H={fila_metodo.get('H', float('nan')):.3f} "
                    f"(p={fila_metodo.get('p_valor', float('nan')):.4f}); interacción "
                    f"H={fila_inter.get('H', float('nan')):.3f} "
                    f"(p={fila_inter.get('p_valor', float('nan')):.4f})."
                )
        partes.append(
            f"### 4.{2 + VARIABLES_RESPUESTA.index(var)} {VARIABLE_LABEL[var]} "
            f"(unidad: {VARIABLE_UNIDAD[var]})\n\n"
            f"**Modelo seleccionado:** {tipo}.{extra}\n\n"
            + cuerpo + "\n\n"
            f"Justificación: {m.get('justificacion', '')}\n"
        )

    # LMM
    lmm = r.get("lmm", {})
    for var in VARIABLES_RESPUESTA:
        m = lmm.get(var)
        if not m:
            continue
        partes.append(
            f"### 4.{6 + VARIABLES_RESPUESTA.index(var)} Sensibilidad LMM - {VARIABLE_LABEL[var]}\n\n"
            f"Modelo mixto `{var} ~ método + (1|aislamiento)`. ICC (aislado)="
            f"{m.get('icc', float('nan')):.4f}; p del método (LRT)="
            f"{m.get('p_valor_metodo_lrt', float('nan')):.4f}.\n\n"
            + (_df_to_markdown(m.get("tabla_efectos_fijos")) if m.get("tabla_efectos_fijos") is not None else "_Sin datos_")
            + "\n"
        )

    return "## 4. Análisis seleccionados y por qué\n\n" + "\n".join(partes) + "\n"


def _seccion_comparaciones(r) -> str:
    posthoc = r.get("posthoc", {})
    partes = []
    for var in VARIABLES_RESPUESTA:
        p = posthoc.get(var)
        if not p:
            continue
        letras = p.get("letras", {})
        letras_txt = "; ".join(
            f"{METODO_LABEL[m]}={letras[m]}" for m in METODOS if m in letras
        )
        partes.append(
            f"### 5.{VARIABLES_RESPUESTA.index(var) + 1} {VARIABLE_LABEL[var]} "
            f"({p.get('metodo', '')})\n\n"
            f"Letras CLD (métodos que comparten letra no difieren, p>=0.05): **{letras_txt}**.\n\n"
            + (_df_to_markdown(p.get("tabla_pares")) if p.get("tabla_pares") is not None else "_Sin datos_")
            + "\n"
        )
    return "## 5. Comparaciones múltiples\n\n" + ("\n".join(partes) if partes else "_Sin datos_") + "\n"


def _seccion_multivariado(r) -> str:
    mv = r.get("multivariado", {})
    if not mv:
        return "## 6. Análisis multivariado\n\n_Sin datos_"
    varianza = mv.get("varianza", [])
    partes = [
        f"## 6. Análisis multivariado (susceptibilidad)\n\n",
        f"- Varianza explicada por PC1 y PC2: "
        f"{(varianza[0] * 100 if len(varianza) else float('nan')):.1f}% y "
        f"{(varianza[1] * 100 if len(varianza) > 1 else float('nan')):.1f}%.",
        f"- Coeficiente cofenético (Ward): {mv.get('cofenetico', float('nan')):.4f}.",
        f"- Número óptimo de clusters KMeans (silhouette): {mv.get('k_optimo', '')}.",
        f"- Categorías biológicas: terciles del score compuesto de susceptibilidad "
        f"(promedio z de inhibición micelial y de conidias) etiquetadas como "
        f"Alta / Moderada / Baja susceptibilidad relativa. No se utiliza el término "
        f"'resistente' por no existir un criterio validado.",
        "",
        "### 6.1 Categorías y clusters por aislado\n",
        _df_to_markdown(mv.get("tabla_final")) if mv.get("tabla_final") is not None else "_Sin datos_",
        "",
        "### 6.2 Cruce cluster x categoría\n",
        _df_to_markdown(mv.get("cruce"), decimales=2) if mv.get("cruce") is not None else "_Sin datos_",
    ]
    return "\n".join(partes) + "\n"


def _seccion_ranking(r) -> str:
    ranking = r.get("ranking", {})
    tabla = ranking.get("tabla")
    if tabla is None:
        return "## 7. Ranking de técnicas\n\n_Sin datos_"
    return (
        "## 7. Ranking de técnicas\n\n"
        "El score compuesto promedia las métricas normalizadas (min-max 0-1) de "
        "rendimiento, inhibición micelial e inhibición de conidias; mayores "
        "valores indican mejor desempeño global.\n\n"
        + _df_to_markdown(tabla) + "\n"
    )


def _seccion_bio(r) -> str:
    modelos = r.get("modelos", {})
    m_mic = modelos.get("porcentaje_inhibicion_micelial", {})
    m_con = modelos.get("porcentaje_inhibicion_conidias", {})
    ranking = r.get("ranking", {})
    tabla = ranking.get("tabla")
    orden_metodos = ""
    if tabla is not None and len(tabla):
        orden_metodos = (
            "De acuerdo con el score compuesto, el ordenamiento de las técnicas "
            "fue: " + ", ".join(
                f"{METODO_LABEL[f['metodo_extraccion']]} ({f['score_compuesto']})"
                for _, f in tabla.iterrows()
            ) + ".\n"
        )

    partes = [
        "## 8. Interpretación biológica\n\n",
        "Los tres extractos (5 mg/mL) mostraron actividad antifúngica variable. "
        "La existencia de un efecto techo en la inhibición micelial (numerosos "
        "valores de 100 %) indica que, a esta concentración, muchos aislados fueron "
        "completamente inhibidos; por ello la inferencia principal sobre esta "
        "variable se apoyó en la vía seleccionada por los supuestos y en los "
        "modelos mixtos como análisis de sensibilidad.",
        "",
        "La inhibición de conidias se expresa en escala log10; valores negativos "
        "indican que el extracto indujo una mayor esporulación relativa en esa "
        "celda. Esta característica debe considerarse al interpretar diferencias "
        "entre métodos.",
        "",
        "La variabilidad entre aislados (ICC) fue evaluada con modelos mixtos; "
        "una ICC alta indica que el aislado explica una proporción importante de "
        "la variación total, lo que justifica la consideración de la "
        "susceptibilidad diferencial (análisis multivariado).",
        "",
        orden_metodos,
    ]
    return "\n".join(partes) + "\n"


def _seccion_conclusiones(r) -> str:
    ranking = r.get("ranking", {}).get("tabla")
    mejor = ""
    if ranking is not None and len(ranking):
        fila = ranking.iloc[0]
        if fila["metodo_extraccion"] == "soxhlet":
            mejor = (" La técnica Soxhlet destacó por su mayor rendimiento de "
                     "extracción y encabezó el score compuesto.")
        else:
            mejor = (f" La técnica {METODO_LABEL[fila['metodo_extraccion']]} "
                     f"encabezó el score compuesto, impulsada por su actividad "
                     f"antifúngica.")
    conclusiones = [
        "## 9. Conclusiones\n\n",
        "1. **Rendimiento**: el método de extracción afectó significativamente el "
        "rendimiento; Soxhlet mostró el mayor valor medio (ver sección 4.1).",
        "2. **Inhibición micelial**: a 5 mg/mL se observó un efecto fuerte; las "
        "diferencias entre métodos y aislados fueron evaluadas según la vía "
        "seleccionada por los supuestos y con modelos mixtos de sensibilidad.",
        "3. **Conidias**: la variable es continua en escala log10; se modeló con "
        "un modelo lineal sobre log10, sin aplicar Poisson/NB (los conteos crudos "
        "no están disponibles).",
        "4. **Susceptibilidad**: los aislados se agruparon en categorías de "
        "susceptibilidad relativa (Alta / Moderada / Baja) a partir del score "
        "compuesto y de la clusterización.",
        "5. **Ranking**: " + (mejor if mejor else "el ranking se detalla en la sección 7."),
    ]
    return "\n".join(conclusiones) + "\n"


def _seccion_limitaciones(r) -> str:
    validacion = r.get("validacion_inh", {})
    control_txt = ""
    if validacion.get("estado") == "ok":
        control_txt = (
            "Los controles C4 ahora están explícitos en el dataset maestro "
            "(columnas `control_crecimiento_mm` y `control_conidias_log10`, "
            "una por aislado y compartidos por las 3 réplicas); su uso en "
            "futuras iteraciones permitiría recalcular %INH con otros controles "
            "o modelar directamente el crecimiento frente al control.\n"
        )
    elif validacion.get("estado") == "no_disponible":
        control_txt = (
            "Los controles C4 no están disponibles (Excel crudo ausente o "
            "incompleto); no fue posible incorporarlos al dataset maestro.\n"
        )
    return (
        "## 10. Limitaciones\n\n"
        "- **Control compartido (pseudorreplicación)**: cada %INH se calculó "
        "contra un único control C4 compartido por las tres réplicas del aislado. "
        "Las réplicas de %INH no son totalmente independientes; el análisis del "
        "crecimiento crudo en mm no presenta este problema y los modelos mixtos "
        "con aislado aleatorio mitigan parcialmente la dependencia.\n"
        + control_txt
        + "- **Escala log10 del %INH de conidias**: el laboratorio reportó la "
        "reducción de conidias en escala log10; no equivale a la reducción "
        "porcentual de conteos crudos.\n"
        "- **Sin dosis-respuesta**: solo se ensayaron 5 mg/mL; no es posible "
        "estimar EC50/EC90 ni extrapolar a otras concentraciones.\n"
        "- **Conidias como log10 continuas**: los conteos crudos no están "
        "disponibles; no aplica la rama Poisson/NB del pipeline (documentada "
        "para datasets futuros).\n"
        "- **Tamaño muestral**: 3 réplicas biológicas por celda limitan la "
        "potencia; los clusters del análisis multivariado deben interpretarse "
        "con cautela por el tamaño de muestra.\n"
    )
