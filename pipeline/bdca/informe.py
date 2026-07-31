"""Fase 3.7: Generador de informe final MD/HTML para análisis RCBD BDCA.

Combina los resultados de las fases anteriores (resumen del diseño, auditoría,
EDA, supuestos, ANOVA clásico, LMM + ICC, comparaciones post-hoc) en un
informe estructurado y reutilizable. Guarda:
  • md e html (versión limpia sin títulos) en ``bdca/resultados/informe/``
  • tablas en ``bdca/resultados/tablas/``
  • figuras en ``bdca/resultados/figuras/``

Cada variable derivada documenta fuente, fórmula y razón.
"""

from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path

from pipeline.config import (
    guardar_tabla,
    save_figure_pub,
    exportar_excel,
    DIR_FIGURAS,
    DIR_TABLAS,
    DIR_REPORTES,
)


def _generar_resumen_disenio(df: pd.DataFrame) -> str:
    """Devuelve un resumen del diseño RCBD.

    Documenta factores, estructura de observaciones y balance.
    """
    tratamientos = sorted(df["trt"].unique())
    bloques = sorted(df["block"].unique())
    n_total = len(df)
    per_trt = df["trt"].value_counts().to_dict()
    per_block = df["block"].value_counts().to_dict()

    resumen = (
        "# Resumen del diseño RCBD (BDCA)"
        "\n\n"
        f"- **Factor:** Tratamiento (cuatro niveles: {', '.join(tratamientos)})."
        f"- **Factor:** Bloque (nueve niveles: {', '.join(bloques)})."
        f"- **Unidades experimentales:** {n_total} observaciones."
        f"- **Distribución por tratamiento:** {pd.Series(per_trt).to_dict()}."
        f"- **Distribución por bloque:** {pd.Series(per_block).to_dict()}."
        f"- **Balance:** Un cultivo por celda trt × bloque (estándar RCBD)."
        f"- **Variable respuesta:** rendimiento (yield) del cultivo en mg."
        "\n"
    )
    return resumen


def _generar_reporte_auditoria(conclusion: dict) -> str:
    """Genera una sección de texto estructurada de la auditoría de carga."""
    linea = "---"
    secciones = []

    secciones.append("## Auditoría de carga")
    secciones.append(linea)
    secciones.append(f"Total de filas: {conclusion['total_filas']}")
    secciones.append(f"Filas únicas (sin duplicados): {conclusion['filas_unicas']}")
    secciones.append(f"Filas con valores NA: {conclusion['filas_con_na']}")

    bal = conclusion.get("balance_celulas", {})
    if bal:
        secciones.append(f"\nBalance por celda:")
        secciones.append(f"- Total de celdas: {bal.get('total_celdas', 'N/A')}")
        secciones.append(f"- Celdas balanceadas: {bal.get('celdas_balanceadas', 'N/A')}")
        secciones.append(f"- Celdas desbalanceadas: {bal.get('celdas_desbalanceadas', 'N/A')}")
        conteo = bal.get("conteo_por_celda", {})
        if conteo:
            secciones.append(f"\nConteo por celda (trt x block):")
            for (trt, block), cnt in sorted(conteo.items()):
                secciones.append(f"  - {trt} x {block}: {cnt}")

    secciones.append(f"\nConteo por tratamiento: {conclusion.get('conteo_por_trt', {})}")
    secciones.append(f"Conteo por bloque: {conclusion.get('conteo_por_block', {})}")

    extremos = conclusion.get("extremos_yield", {})
    if extremos:
        secciones.append(f"\nExtremos de rendimiento:")
        secciones.append(f"- Mínimo: {extremos.get('minimo', 'N/A')}")
        secciones.append(f"- Máximo: {extremos.get('maximo', 'N/A')}")
        secciones.append(f"- Media: {extremos.get('media', 'N/A')}")
        secciones.append(f"- Std: {extremos.get('std', 'N/A')}")

    secciones.append(f"\nNúmero de cultivos con rendimiento negativo: {conclusion.get('yield_negativo', 0)}")
    secciones.append(f"\n¿Auditoría aprobada? **{'Sí' if conclusion.get('auditado_ok', False) else 'No'}**")
    return "\n".join(secciones)


def _generar_reporte_eda(tabla_descriptivos: pd.DataFrame, fig_paths: list[str]) -> str:
    """Genera una sección MD del informe a partir de la tabla descriptiva y rutas de figuras."""
    secciones = ["## Exploración descriptiva (EDA)", "---"]

    secciones.append("\n### Resumen descriptivo por tratamiento")
    with pd.option_context("display.max_columns", None, "display.width", 180):
        secciones.append(tabla_descriptivos.to_markdown(index=False))

    if fig_paths:
        secciones.append("\n### Figuras exploratorias generadas")
        for fig_path in fig_paths:
            nombre = Path(fig_path).stem
            secciones.append(f'- `{nombre}`: "{nombre}.png"')

    return "\n".join(secciones)


def _generar_reporte_supuestos(resultado_supuestos: dict) -> str:
    """Genera la sección MD del informe para supuestos.

    Documenta la justificación de la ruta de inferencia elegida.
    """
    secciones = ["## Verificación de supuestos", "---"]

    tabla = resultado_supuestos["tabla_supuestos"]
    secciones.append("\n### Tabla de supuestos")
    with pd.option_context("display.max_columns", None, "display.width", 180):
        secciones.append(tabla.to_markdown(index=False))

    secciones.append(f"\n### Decisión: {resultado_supuestos['tipo_modelo'].upper()}")
    secciones.append(f"**Justificación:** {resultado_supuestos['justificacion']}")

    return "\n".join(secciones)


def _generar_reporte_anova(resultado_anova: dict) -> str:
    """Genera la sección MD del informe para ANOVA de bloques.

    Muestra la tabla ANOVA, tamaños de efecto (eta2 parcial) y advertencia sobre
    additividad (no testable por diseño).
    """
    secciones = ["## ANOVA clásico de bloques RCBD", "---"]
    tabla = resultado_anova["tabla_anova"]
    secciones.append("\n### Tabla ANOVA (tipo II)")
    with pd.option_context("display.max_columns", None, "display.width", 180):
        secciones.append(tabla.to_markdown(index=False))

    eta2 = resultado_anova["eta2_parcial"]
    if eta2:
        secciones.append("\n### Tamaños de efecto (eta2 parcial)")
        for fuente, valor in eta2.items():
            secciones.append(f"- **{fuente}**: {valor:.4f}")

    return "\n".join(secciones)


def _generar_reporte_lmm(resultado_lmm: dict) -> str:
    """Genera la sección MD del informe para Modelo Mixto Lineal (LMM)."""
    secciones = ["## Modelo mixto lineal (LMM) de bloques RCBD", "---"]

    tabla = resultado_lmm["tabla_fija"]
    secciones.append("\n### Efectos fijos (tratamiento)")
    with pd.option_context("display.max_columns", None, "display.width", 180):
        secciones.append(tabla.to_markdown(index=False))

    secciones.append(f"\n### Varianzas e ICC")
    secciones.append(f"- Varianza de bloque: {resultado_lmm['var_bloque']:.4f}")
    secciones.append(f"- Varianza residual: {resultado_lmm['var_residual']:.4f}")
    secciones.append(f"- ICC: {resultado_lmm['icc']:.4f}")

    if "limitacion_aditividad" in resultado_lmm:
        secciones.append(f"\n### Limitación de aditividad")
        secciones.append(f"{resultado_lmm['limitacion_aditividad']}")

    return "\n".join(secciones)


def _generar_reporte_posthoc(resultado_posthoc: pd.DataFrame) -> str:
    """Genera la sección MD del informe para comparaciones post-hoc Tukey HSD.

    Resalta las comparaciones vs el control T0 (sin tratar).
    """
    secciones = ["## Comparaciones múltiples post-hoc Tukey HSD", "---"]

    secciones.append("\n### Resultados Tukey HSD (referencia: T0, control sin tratar)")
    with pd.option_context("display.max_columns", None, "display.width", 180):
        secciones.append(resultado_posthoc.to_markdown(index=False))

    # Filtrar para mostrar solo las comparaciones vs T0
    vs_r = resultado_posthoc[resultado_posthoc["vs_referencia_T0"] == True]
    if not vs_r.empty:
        secciones.append("\n### Comparaciones vs Control T0")
        with pd.option_context("display.max_columns", None, "display.width", 180):
            secciones.append(vs_r.to_markdown(index=False))

    return "\n".join(secciones)


def _documentar_derivados() -> str:
    """Documenta todas las variables derivadas, fuentes, fórmulas y razones."""
    secciones = ["## Documentación de variables derivadas", "---"]

    documentos = [
        ("EDA -> medias, IC95%", "pipeline/bdca/eda.py -> resumen_descriptivo", "Cada tratamiento: media de yield, IC95 via t de Student (df=n-1)"),
        ("ANOVA -> eta2 parcial", "pipeline/bdca/modelos.py -> _calcular_eta2_parcial", "Eta2 parcial = SS_efecto / (SS_efecto + SS_residual)"),
        ("LMM -> ICC", "pipeline/bdca/modelos.py -> lmm_bloques", "ICC = var_bloque / (var_bloque + var_residual) para estimar correlación intra-clase de bloque"),
        ("Post-hoc -> Tukey", "pipeline/bdca/comparaciones.py -> posthoc_tukey", "Tukey HSD por pares de tratamientos (statsmodels.pairwise_tukeyhsd) con columna 'vs_referencia_T0'"),
        ("Supuestos -> tabla", "pipeline/bdca/supuestos.py -> analisis_supuestos", "Tabla de supuestos de normalidad (Shapiro-Wilk), homocedasticidad (Levene), independencia serial (Durbin-Watson)"),
    ]

    for nombre, fuente, descripcion in documentos:
        secciones.append(f"\n### {nombre}")
        secciones.append(f"**Fuente:** `{fuente}`")
        secciones.append(f"**Fórmula/Justificación:** {descripcion}")

    return "\n".join(secciones)


def informe_completo(
    df: pd.DataFrame,
    auditoria: dict,
    descriptivos: pd.DataFrame,
    fig_paths_eda: list[str],
    resultado_supuestos: dict,
    resultado_anova: dict,
    resultado_lmm: dict,
    resultado_posthoc: pd.DataFrame,
) -> str:
    """Genera el informe final MD unificado.

    Guarda en ``DIR_REPORTES / "informe_final.md"`` y produce una versión HTML.
    """
    secciones = []

    # 1. Resumen del diseño
    secciones.append(_generar_resumen_disenio(df))

    # 2. Auditoría de carga
    secciones.append(_generar_reporte_auditoria(auditoria))

    # 3. EDA
    secciones.append(_generar_reporte_eda(descriptivos, fig_paths_eda))

    # 4. Supuestos
    secciones.append(_generar_reporte_supuestos(resultado_supuestos))

    # 5. ANOVA clásico
    secciones.append(_generar_reporte_anova(resultado_anova))

    # 6. LMM
    secciones.append(_generar_reporte_lmm(resultado_lmm))

    # 7. Post-hoc
    secciones.append(_generar_reporte_posthoc(resultado_posthoc))

    # 8. Documentación de variables derivadas
    secciones.append(_documentar_derivados())

    # Unir todo
    contenido = "\n\n".join(secciones)

    # Guardar MD
    md_path = DIR_REPORTES / "informe_final.md"
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text(contenido, encoding="utf-8")
    print(f"✓ Informe MD guardado en {md_path}")

    # Convertir a HTML (básico)
    html_path = DIR_REPORTES / "informe_final.html"
    try:
        from markdown import markdown
        html = markdown(contenido)
        html = (
            """<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Informe final - BDCA</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        h1, h2, h3 { color: #1f3864; }
        table { border-collapse: collapse; width: 100%; margin: 10px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
""".join(html.split("\n"))
        )
        html_path.write_text(html, encoding="utf-8")
        print(f"✓ Informe HTML guardado en {html_path}")
    except ImportError:
        print("[ADVERTENCIA] markdown no está instalado; omitiendo HTML.")

    # Guardar la tabla maestra del informe (una fila por sección con metadatos)
    tabla_reporte = pd.DataFrame([
        {"seccion": "diseño", "descripcion": "Resumen del diseño RCBD", "ruta_archivo": str(md_path)},
        {"seccion": "auditoria", "descripcion": "Resultados de auditoría de carga", "ruta_archivo": str(md_path)},
        {"seccion": "eda", "descripcion": "Tabla descriptiva y figuras exploratorias", "figuras_generadas": fig_paths_eda},
        {"seccion": "supuestos", "descripcion": "Tabla de verificación de supuestos", "tipo_inferencia": resultado_supuestos["tipo_modelo"]},
        {"seccion": "anova", "descripcion": "Tabla ANOVA clásico de bloques", "fuente_datos": "modelo OLS yield ~ C(trt) + C(block)"},
        {"seccion": "lmm", "descripcion": "Tabla de efectos fijos LMM y varianzas/ICC", "fuente_datos": "modelo mixto yield ~ C(trt) | block"},
        {"seccion": "posthoc", "descripcion": "Resultados Tukey HSD (vs R)", "n_pares": len(resultado_posthoc)},
    ])
    guardar_tabla(tabla_reporte, "reporte_informe", index=False)

    # Exportar tablas clave a un libro Excel (tarea 3.7: outputs excel)
    try:
        exportar_excel(
            {
                "auditoria": pd.DataFrame([auditoria]),
                "descriptivos": descriptivos,
                "supuestos": resultado_supuestos["tabla_supuestos"],
                "anova": resultado_anova["tabla_anova"],
                "lmm_fijos": resultado_lmm["tabla_fija"],
                "posthoc": resultado_posthoc,
            },
            ruta=DIR_REPORTES.parent / "excel" / "resumen_bdca.xlsx",
        )
        print(f"✓ Libro Excel guardado en {DIR_REPORTES.parent / 'excel' / 'resumen_bdca.xlsx'}")
    except Exception as e:  # noqa: BLE001 — no romper el informe si falla el export
        print(f"[ADVERTENCIA] No se pudo exportar Excel: {e}")

    return contenido


if __name__ == "__main__":
    # Permite ejecución directa para pruebas rápidas
    from pipeline.bdca.cargar import cargar
    from pipeline.bdca.eda import resumen_descriptivo, figuras_eda
    from pipeline.bdca.supuestos import analisis_supuestos
    from pipeline.bdca.modelos import anova_bloques, lmm_bloques
    from pipeline.bdca.comparaciones import posthoc_tukey

    df, auditoria = cargar()
    print("Generando informe completo para BDCA...")
    contenido = informe_completo(
        df=df,
        auditoria=auditoria,
        descriptivos=resumen_descriptivo(df),
        fig_paths_eda=[],
        resultado_supuestos=analisis_supuestos(df),
        resultado_anova=anova_bloques(df),
        resultado_lmm=lmm_bloques(df),
        resultado_posthoc=posthoc_tukey(df),
    )
    # Ejecutar EDA figuras manualmente para prueba rápida (se eligen solo 2 figuras para probar)
    # En realidad usamos un enfoque real: solo las necesarias para prueba rápida.
    print("Ejecución de EDA para obtener rutas de figuras...")
    fig_paths = figuras_eda(df)
    # Reescribir el informe usando las rutas de figuras
    contenido = informe_completo(
        df=df,
        auditoria=auditoria,
        descriptivos=resumen_descriptivo(df),
        fig_paths_eda=fig_paths,
        resultado_supuestos=analisis_supuestos(df),
        resultado_anova=anova_bloques(df),
        resultado_lmm=lmm_bloques(df),
        resultado_posthoc=posthoc_tukey(df),
    )
