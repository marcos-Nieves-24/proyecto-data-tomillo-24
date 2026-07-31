"""Fases 1 y 2: carga de datos y auditoria de calidad.

Fase 1 (carga): lectura del dataset consolidado tidy (Excel, hoja
'Consolidado') y del rendimiento de extraccion (CSV), con verificacion de
dimensiones, tipos y vista previa.

Fase 2 (auditoria): identificacion de valores faltantes, duplicados,
columnas constantes, valores inconsistentes y atipicos por rango
intercuartilico (IQR). Los atipicos se marcan pero NUNCA se eliminan
(regla de integridad de datos del proyecto).
"""

from __future__ import annotations

import pandas as pd

from pipeline.config import (
    CSV_RENDIMIENTO,
    EXCEL_TIDY,
    METODO_LABEL,
    guardar_tabla,
)

# Restricciones de valores esperados por columna (formato crudo del Excel).
COLUMNAS_VALIDACION = {
    "Metodo de extraccion": set(METODO_LABEL.values()),
    "Replica": {"R1", "R2", "R3"},
}
LIMITES_RANGOS = {
    "%INH micelial": (0.0, 100.0),
    "Crecimiento micelial (mm)": (0.0, 100.0),
    "Conidias (log10/ml)": (0.0, None),  # solo limite inferior
}


def cargar_datos() -> dict:
    """Carga los datos del bioensayo y del rendimiento de extraccion.

    Devuelve un diccionario con 'bio' (DataFrame) y 'rend' (DataFrame) e
    imprime un resumen de dimensiones, columnas, tipos y vista previa.
    """
    bio = pd.read_excel(EXCEL_TIDY, sheet_name="Consolidado")
    rend = pd.read_csv(CSV_RENDIMIENTO)

    print("=" * 72)
    print("FASE 1 - Carga de datos")
    print("=" * 72)
    for nombre, df in (("Bioensayo", bio), ("Rendimiento", rend)):
        print(f"\n[{nombre}]")
        print(f"  Dimensiones: {df.shape[0]} filas x {df.shape[1]} columnas")
        print(f"  Columnas:    {list(df.columns)}")
        print(f"  Tipos:       {dict(df.dtypes.astype(str))}")
        print("  Vista previa (5 filas):")
        with pd.option_context("display.max_columns", None, "display.width", 180):
            print(df.head(5).to_string())

    return {"bio": bio, "rend": rend}


def auditoria_calidad(df, nombre="bioensayo") -> tuple[pd.DataFrame, dict]:
    """Audita la calidad de un DataFrame columna a columna.

    Reporta n de no nulos, porcentaje de faltantes, filas duplicadas exactas,
    columnas constantes, tipos, valores inconsistentes (chequeos por columna)
    y atipicos por rango intercuartilico (1.5xIQR).

    Los atipicos se flaguean pero no se eliminan. Devuelve (tabla, resumen).
    """
    filas = df.shape[0]
    resumen = {
        "n_filas": int(filas),
        "n_columnas": int(df.shape[1]),
        "n_filas_duplicadas": int(df.duplicated().sum()),
    }

    registros = []
    for columna in df.columns:
        serie = df[columna]
        n_no_nulos = int(serie.notna().sum())
        pct_faltantes = float((1 - n_no_nulos / filas) * 100) if filas else 0.0
        es_constante = bool(serie.nunique(dropna=False) <= 1)
        inconsistentes = _contar_inconsistentes(df, columna)
        n_atipicos = _contar_atipicos_iqr(serie)
        if n_atipicos:
            resumen.setdefault("atipicos_por_iqr", {})[columna] = n_atipicos

        registros.append({
            "variable": columna,
            "tipo_dato": str(serie.dtype),
            "n_no_nulos": n_no_nulos,
            "pct_faltantes": round(pct_faltantes, 2),
            "columna_constante": es_constante,
            "n_filas_duplicadas": resumen["n_filas_duplicadas"],
            "n_valores_inconsistentes": inconsistentes,
            "n_atipicos_iqr": n_atipicos,
        })

    tabla = pd.DataFrame(registros)
    guardar_tabla(tabla, "auditoria_calidad", index=False)

    print("=" * 72)
    print("FASE 2 - Auditoría de calidad de datos")
    print("=" * 72)
    print(
        f"Filas: {resumen['n_filas']} | Columnas: {resumen['n_columnas']} | "
        f"Filas duplicadas exactas: {resumen['n_filas_duplicadas']}"
    )
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print(tabla.to_string(index=False))
    for col, n in resumen.get("atipicos_por_iqr", {}).items():
        print(
            f"  [IQR] Columna '{col}': {n} valor(es) fuera de 1.5xIQR "
            "(flaggeados, no eliminados)."
        )

    return tabla, resumen


def _contar_inconsistentes(df: pd.DataFrame, columna: str) -> int:
    """Cuenta valores que violan las restricciones de rango/set por columna."""
    serie = df[columna]
    nulos = serie.isna()
    if columna in COLUMNAS_VALIDACION:
        permitidos = COLUMNAS_VALIDACION[columna]
        return int((~nulos & ~serie.isin(permitidos)).sum())
    if columna in LIMITES_RANGOS:
        inf, sup = LIMITES_RANGOS[columna]
        viola_inf = serie < inf
        viola_sup = serie > sup if sup is not None else pd.Series(False, index=serie.index)
        return int((~nulos & (viola_inf | viola_sup)).sum())
    return 0


def _contar_atipicos_iqr(serie) -> int:
    """Cuenta valores fuera de 1.5xIQR para columnas numericas."""
    if not pd.api.types.is_numeric_dtype(serie):
        return 0
    s = serie.dropna()
    if s.empty:
        return 0
    q1, q3 = s.quantile([0.25, 0.75])
    iqr = q3 - q1
    if iqr == 0:
        return 0
    lim_inf = q1 - 1.5 * iqr
    lim_sup = q3 + 1.5 * iqr
    return int(((s < lim_inf) | (s > lim_sup)).sum())
