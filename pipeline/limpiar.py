"""Fase 3: limpieza y normalizacion del dataset maestro.

Renombra las columnas al esquema canonico en ingles, normaliza los niveles
de los factores, valida la estructura (279 filas, sin nulos, balanceado),
incorpora los controles C4 del Excel crudo (columna control_*) y guarda el
dataset maestro (CSV y XLSX) junto con el diccionario de datos.

También valida de forma informativa la fórmula del %INH reportado por el
laboratorio frente a (1 - C1/C4) * 100.

Provenance de datos crudos:
- DCA: Excel con controles C4 (`datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx`)
- BDCA: CSV `DBCA_Jenkyn_control_mildeo.csv` (`datos_crudos/bdca/`); columnas `plot, trt, block, yield`
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline.config import (
    COLUMNAS_BIOENSAYO,
    DICCIONARIO_MD,
    EXCEL_CRUDO,
    MASTER_CSV,
    MASTER_XLSX,
    METODOS,
    METODO_LABEL,
    guardar_tabla,
    normalizar_texto,
)

RENOMBRES_BIO = {
    "Metodo de extraccion": "metodo_extraccion",
    "Aislamiento": "aislamiento",
    "Replica": "replica",
    "Crecimiento micelial (mm)": "crecimiento_micelial_mm",
    "%INH micelial": "porcentaje_inhibicion_micelial",
    "Conidias (log10/ml)": "conidias_log10_ml",
    "%INH conidias": "porcentaje_inhibicion_conidias",
}

VARIABLES_NUMERICAS_BIO = [
    "crecimiento_micelial_mm",
    "porcentaje_inhibicion_micelial",
    "conidias_log10_ml",
    "porcentaje_inhibicion_conidias",
]

COLUMNAS_CONTROL = ["control_crecimiento_mm", "control_conidias_log10"]

# Disponibilidad de los controles C4 (se actualiza en enriquecer_con_controles).
# La validacion del %INH se saltea si es False.
CONTROLES_DISPONIBLES = False

# Layout del Excel crudo (verificado sobre el archivo; no inventar indices).
# MACERACIÓN y SOXHLET: columna 0 = Aislado, 4 = C4 crecimiento (mm),
# 8 = C4 conidias (log10). ULTRASONIDO (con columna '#'): 1, 5 y 9.
_LAYOUT_HOJAS = {
    "MACERACIÓN": {"col_aislado": 0, "col_c4_mm": 4, "col_c4_log": 8},
    "SOXHLET": {"col_aislado": 0, "col_c4_mm": 4, "col_c4_log": 8},
    "ULTRASONIDO": {"col_aislado": 1, "col_c4_mm": 5, "col_c4_log": 9},
}
_FILAS_ENCABEZADO = 3  # los datos comienzan en la fila 4 (indice 3)


def enriquecer_con_controles(bio: pd.DataFrame) -> pd.DataFrame:
    """Agrega las columnas de control C4 al bioensayo (join por aislado × método).

    Lee el Excel crudo del laboratorio, construye una tabla de una fila por
    aislado × método (93 filas) con ``control_crecimiento_mm`` y
    ``control_conidias_log10``, y la hace join con ``bio``.

    El valor C4 aparece solo en la réplica 1 de cada aislado (o está repetido);
    se propaga a las 3 réplicas con forward-fill dentro de cada hoja. Si el
    archivo o alguna hoja faltan, imprime una advertencia, deja ``bio`` sin
    cambios y fija ``CONTROLES_DISPONIBLES = False``.
    """
    global CONTROLES_DISPONIBLES
    CONTROLES_DISPONIBLES = False

    if not EXCEL_CRUDO.exists():
        print(f"[ADVERTENCIA] No se encontro el Excel crudo {EXCEL_CRUDO}. "
              f"Los controles C4 no se incorporan y la validacion del %INH se saltea.")
        return bio.copy()

    try:
        hojas = pd.ExcelFile(EXCEL_CRUDO).sheet_names
    except Exception as exc:
        print(f"[ADVERTENCIA] No se pudo leer el Excel crudo: {exc}. "
              f"Los controles C4 no se incorporan.")
        return bio.copy()

    faltantes = [sh for sh in _LAYOUT_HOJAS if sh not in hojas]
    if faltantes:
        print(f"[ADVERTENCIA] Faltan hojas en el Excel crudo: {faltantes}. "
              f"Los controles C4 no se incorporan y la validacion se saltea.")
        return bio.copy()

    frames = []
    for hoja, lay in _LAYOUT_HOJAS.items():
        df = pd.read_excel(EXCEL_CRUDO, sheet_name=hoja, header=None)
        df = df.iloc[_FILAS_ENCABEZADO:].copy()
        aisl = df[lay["col_aislado"]].ffill().astype(str).str.strip()
        tmp = pd.DataFrame({
            "metodo_extraccion": normalizar_texto(hoja),
            "aislamiento": aisl,
            "control_crecimiento_mm": pd.to_numeric(df[lay["col_c4_mm"]], errors="coerce").ffill(),
            "control_conidias_log10": pd.to_numeric(df[lay["col_c4_log"]], errors="coerce").ffill(),
        })
        tmp = tmp.dropna(subset=["aislamiento"])
        frames.append(tmp)

    controles = pd.concat(frames, ignore_index=True)
    controles = controles[controles["metodo_extraccion"].isin(METODOS)]

    # El layout produce una fila por réplica (93 filas × 3 réplicas = 279);
    # se deduplica a un control único por aislado × método (93 filas) para
    # que el join con el bioensayo propague el mismo control a las 3 réplicas.
    controles = controles.drop_duplicates(
        subset=["metodo_extraccion", "aislamiento"], keep="first"
    )

    if len(controles) != 93:
        print(f"[ADVERTENCIA] Se esperaban 93 filas únicas de controles, se "
              f"obtuvieron {len(controles)}. Los controles no se incorporan.")
        return bio.copy()

    # Verificación de consistencia: un único control por aislado × método.
    n_unicos = controles.groupby(["metodo_extraccion", "aislamiento"])[
        ["control_crecimiento_mm", "control_conidias_log10"]
    ].nunique().max().max()
    if n_unicos != 1:
        print(f"[ADVERTENCIA] Los controles no son constantes por aislado × "
              f"método (se observaron hasta {n_unicos} valores distintos). "
              f"Los controles no se incorporan.")
        return bio.copy()

    enriquecido = bio.merge(controles, on=["metodo_extraccion", "aislamiento"], how="left")

    nulos_control = int(enriquecido[COLUMNAS_CONTROL].isna().sum().sum())
    if nulos_control:
        print(f"[ADVERTENCIA] Tras el join quedaron {nulos_control} valores nulos "
              f"en las columnas de control; se conservan como NaN.")
    else:
        print(f"  Controles C4 incorporados: {len(controles)} filas, sin nulos.")

    CONTROLES_DISPONIBLES = True
    return enriquecido


def normalizar_master(bio: pd.DataFrame, rend: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normaliza los DataFrames crudos al esquema canonico.

    Aplica renombrado, normalizacion de texto (minusculas sin acentos),
    replica como entero 1-3, tipos numericos, validacion estricta
    (279 filas, sin valores nulos, diseno balanceado) y enriquece el
    bioensayo con los controles C4 (9 columnas en total).
    """
    df_bio = bio.rename(columns=RENOMBRES_BIO)
    df_bio["metodo_extraccion"] = df_bio["metodo_extraccion"].map(normalizar_texto)
    df_bio["aislamiento"] = df_bio["aislamiento"].str.strip()
    df_bio["replica"] = (
        df_bio["replica"].str.strip().str.replace("R", "", regex=False).astype(int)
    )
    for col in VARIABLES_NUMERICAS_BIO:
        df_bio[col] = pd.to_numeric(df_bio[col], errors="coerce")
    df_bio = df_bio[COLUMNAS_BIOENSAYO].reset_index(drop=True)

    df_rend = rend.copy()
    df_rend["metodo_extraccion"] = df_rend["metodo_extraccion"].map(normalizar_texto)
    df_rend["replica_biologica"] = pd.to_numeric(
        df_rend["replica_biologica"], errors="coerce"
    ).astype(int)

    _validar(df_bio, df_rend)
    df_bio = enriquecer_con_controles(df_bio)
    return df_bio, df_rend


def _validar(df_bio: pd.DataFrame, df_rend: pd.DataFrame) -> bool:
    """Valida estructura, completitud y balanceo del dataset canonico."""
    errores = []

    if len(df_bio) != 279:
        errores.append(f"Se esperaban 279 filas de bioensayo, se obtuvieron {len(df_bio)}.")
    if df_bio.isna().any().any():
        errores.append("Existen valores nulos en el dataset de bioensayo.")
    if df_rend.isna().any().any():
        errores.append("Existen valores nulos en el dataset de rendimiento.")

    n_aislados = df_bio["aislamiento"].nunique()
    conteo = df_bio.groupby(["metodo_extraccion", "aislamiento"]).size()
    balanceado = bool((conteo == 3).all())
    if not balanceado:
        errores.append("El diseño no está balanceado: no todas las celdas método × aislado tienen 3 réplicas.")
    if n_aislados != 31:
        errores.append(f"Se esperaban 31 aislados, se obtuvieron {n_aislados}.")
    if set(df_bio["metodo_extraccion"].unique()) != set(METODOS):
        errores.append(
            f"Niveles de método inesperados: {sorted(df_bio['metodo_extraccion'].unique())}."
        )

    if set(df_rend["metodo_extraccion"].unique()) != set(METODOS):
        errores.append(
            f"Niveles de método en rendimiento inesperados: {sorted(df_rend['metodo_extraccion'].unique())}."
        )
    if len(df_rend) != 9:
        errores.append(f"Se esperaban 9 filas de rendimiento, se obtuvieron {len(df_rend)}.")

    if errores:
        raise ValueError("Validacion del dataset maestro fallo:\n  - " + "\n  - ".join(errores))
    return True


def guardar_master(df_bio: pd.DataFrame, df_rend: pd.DataFrame) -> tuple:
    """Guarda el dataset maestro (CSV y XLSX) y el diccionario de datos."""
    df_bio.to_csv(MASTER_CSV, index=False, encoding="utf-8")
    with pd.ExcelWriter(MASTER_XLSX, engine="openpyxl") as writer:
        df_bio.to_excel(writer, sheet_name="Bioensayo", index=False)
        df_rend.to_excel(writer, sheet_name="Rendimiento", index=False)

    lineas = [
        "# Diccionario de datos - Master dataset (pipeline reproducible)",
        "",
        "Fuente: dca/resultados/database/consolidado_tidy.xlsx (hoja 'Consolidado') y",
        "dca/resultados/database/rendimiento_extraccion.csv.",
        "",
        "## Hoja Bioensayo",
        "",
        "| Variable | Unidad | Descripción |",
        "|----------|--------|-------------|",
        "| metodo_extraccion | - | Técnica de extracción (maceración, soxhlet, ultrasonido). Factor fijo. |",
        "| aislamiento | - | Código del aislado de Fusarium spp. (31 aislados). |",
        "| replica | - | Réplica biológica (1-3). |",
        "| crecimiento_micelial_mm | mm | Crecimiento micelial a 5 mg/mL. |",
        "| porcentaje_inhibicion_micelial | % | Inhibición calculada contra el control C4 del propio aislado. |",
        "| conidias_log10_ml | log10(conidias/mL) | Concentración de conidias en escala log10 (transformada por el laboratorio). |",
        "| porcentaje_inhibicion_conidias | % | Reducción de conidias en escala log10; puede ser negativa (mayor esporulación). |",
        "| control_crecimiento_mm | mm | Control C4 del aislado: crecimiento micelial sin extracto (extraído del Excel crudo del laboratorio). Compartido por las 3 réplicas del aislado. |",
        "| control_conidias_log10 | log10(conidias/mL) | Control C4 del aislado: conidias sin extracto en escala log10 (Excel crudo del laboratorio). Compartido por las 3 réplicas del aislado. |",
        "",
        "## Hoja Rendimiento",
        "",
        "| Variable | Unidad | Descripción |",
        "|----------|--------|-------------|",
        "| metodo_extraccion | - | Técnica de extracción. |",
        "| peso_material_seco_g | g | Peso de material vegetal seco. |",
        "| peso_extracto_obtenido_g | g | Peso de extracto obtenido. |",
        "| rendimiento_pct | % | Rendimiento = (peso_extracto / peso_material_seco) × 100. |",
        "| replica_biologica | - | Réplica biológica (1-3). |",
        "",
        "## Variables derivadas y transformaciones",
        "",
        "- Concentración única ensayada: 5 mg/mL (no es factor experimental).",
        "- %INH micelial: 100 indica inhibición completa (efecto techo).",
        "- %INH conidias: calculado por el laboratorio sobre la escala log10.",
        "- Cada %INH se calculó contra un control C4 compartido por las 3 réplicas",
        "  del aislado (pseudorreplicación del control; ver limitaciones).",
        "- Las columnas control_* se extraen del Excel crudo del laboratorio",
        "  (datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx, hojas por técnica) y se propagan",
        "  a las 3 réplicas de cada aislado.",
        "",
    ]
    DICCIONARIO_MD.write_text("\n".join(lineas), encoding="utf-8")

    print("=" * 72)
    print("FASE 3 - Dataset maestro normalizado y validado")
    print("=" * 72)
    print(f"  Bioensayo: {df_bio.shape[0]} filas x {df_bio.shape[1]} columnas")
    print(f"  Aislados unicos: {df_bio['aislamiento'].nunique()}")
    print(f"  Celdas método × aislado con 3 réplicas: {(df_bio.groupby(['metodo_extraccion', 'aislamiento']).size() == 3).all()}")
    print(f"  Guardado: {MASTER_CSV}")
    print(f"  Guardado: {MASTER_XLSX}")
    print(f"  Diccionario: {DICCIONARIO_MD}")
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print("\nVista previa del dataset maestro:")
        print(df_bio.head().to_string())

    return MASTER_CSV, MASTER_XLSX, DICCIONARIO_MD


def _inh_verificado(c1: pd.Series, control: pd.Series) -> pd.Series:
    """Reconstruye %INH = (1 - C1/C4) * 100 con manejo del control nulo.

    Si ``control == 0``: se define 100 cuando C1 == 0 (inhibición total
    consistente) y NaN cuando C1 > 0 (cociente indefinido). En los datos
    reales el control nunca es 0, pero la función no debe romper.
    """
    out = pd.Series(np.nan, index=c1.index, dtype=float)
    c1 = c1.astype(float)
    control = control.astype(float)
    mascara_ok = control != 0
    out[mascara_ok] = (1.0 - c1[mascara_ok] / control[mascara_ok]) * 100.0
    ceros = control == 0
    out[ceros & (c1 == 0)] = 100.0
    return out


def validar_inh(master: pd.DataFrame) -> dict:
    """Valida de forma informativa la fórmula del %INH reportado por el laboratorio.

    Reconstruye ``(1 - C1/C4) * 100`` con los controles incorporados y compara
    contra las columnas ``porcentaje_inhibicion_*`` reportadas. La validación es
    de integridad de datos: **no reemplaza** las respuestas usadas en la
    inferencia (se usa el %INH reportado, instrucción del investigador).

    Para las conidias, la fórmula se aplica sobre la escala log10 directamente
    (el laboratorio reportó %INH de conidias como reducción en log10); esto fue
    verificado contra los valores reportados.

    Devuelve un dict con estado 'ok' / 'discrepancias' / 'no_disponible' y
    guarda ``dca/resultados/tablas/validacion_inh.csv``.
    """
    if not CONTROLES_DISPONIBLES or not all(
        c in master.columns for c in COLUMNAS_CONTROL
    ):
        nota = ("Los controles C4 no están disponibles (Excel crudo ausente o "
                "incompleto); la validación de la fórmula del %INH no se realizó.")
        print(f"[ADVERTENCIA] FASE 7.5 - Validación de %INH: {nota}")
        return {"estado": "no_disponible", "nota": nota}

    pares = [
        ("porcentaje_inhibicion_micelial", "crecimiento_micelial_mm",
         "control_crecimiento_mm", "micelial"),
        ("porcentaje_inhibicion_conidias", "conidias_log10_ml",
         "control_conidias_log10", "conidias"),
    ]
    filas = []
    detalles = {}
    for var_reportada, var_c1, var_control, etiqueta in pares:
        verificado = _inh_verificado(master[var_c1], master[var_control])
        reportado = master[var_reportada].astype(float)
        diff = (reportado - verificado).abs()
        n_validos = int(diff.notna().sum())
        max_diff = float(diff.max()) if n_validos else float("nan")
        n_discrepancias = int((diff > 1e-6).sum()) if n_validos else 0
        estado = "ok" if (n_discrepancias == 0 and n_validos == len(master)) else (
            "discrepancias" if n_discrepancias else "parcial"
        )
        nota = (
            "La fórmula (1 - C1/C4) x 100 coincide con el %INH reportado "
            "(tolerancia 1e-6)."
            if estado == "ok" else
            "Diferencias detectadas entre la fórmula reconstruida y el %INH reportado."
        )
        filas.append({
            "variable": var_reportada,
            "n_verificadas": n_validos,
            "max_diff_abs": round(max_diff, 10),
            "n_discrepancias": n_discrepancias,
            "estado": estado,
            "nota": nota,
        })
        detalles[etiqueta] = {
            "n_verificadas": n_validos,
            "max_diff_abs": max_diff,
            "n_discrepancias": n_discrepancias,
            "estado": estado,
            "ejemplos": _ejemplos_discrepancias(
                master, var_reportada, var_c1, var_control, verificado
            ),
        }

    tabla = pd.DataFrame(filas)
    guardar_tabla(tabla, "validacion_inh", index=False)

    print("=" * 72)
    print("FASE 7.5 - Validación de la fórmula del %INH (informativa)")
    print("=" * 72)
    print("  La inferencia principal usa el %INH reportado por el laboratorio;")
    print("  esta validación solo verifica la integridad de los datos.")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(tabla.to_string(index=False))

    return {
        "estado": "ok" if all(d["estado"] == "ok" for d in detalles.values()) else
                 ("discrepancias" if any(d["estado"] == "discrepancias" for d in detalles.values()) else "parcial"),
        "tabla": tabla,
        "detalles": detalles,
    }


def _ejemplos_discrepancias(master, var_reportada, var_c1, var_control, verificado) -> list:
    """Devuelve ejemplos (hasta 5) de filas con discrepancia > 1e-6."""
    reportado = master[var_reportada].astype(float)
    diff = (reportado - verificado).abs()
    mascara = diff > 1e-6
    if not mascara.any():
        return []
    cols = ["metodo_extraccion", "aislamiento", "replica",
            var_c1, var_control, var_reportada]
    sub = master.loc[mascara, cols].head(5).copy()
    return sub.to_dict(orient="records")
