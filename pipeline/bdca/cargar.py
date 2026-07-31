"""Fase 3.1: Cargar y auditar CSV de rendimiento BDCA.

Carga ``datos_crudos/bdca/DBCA_Jenkyn_control_mildeo.csv`` (36 filas: plot, trt, block, yield),
audita el diseño RCBD balanceado y documenta cualquier anomalía. No se imputan
valores.
"""

from __future__ import annotations

import pandas as pd
import numpy as np
import os

from pipeline.config import guardar_tabla, DIR_CRUDOS


CSV_BDCA = DIR_CRUDOS / "DBCA_Jenkyn_control_mildeo.csv"


def cargar() -> tuple[pd.DataFrame, dict]:
    """Carga el CSV BDCA y ejecuta la auditoría.

    Retorna:
        df: DataFrame limpio con columnas 'trt', 'block', 'yield' (float).
        auditoria: registro estructurado con contadores, anomalías y resúmenes.
    """
    # ---------- Carga ----------
    if not CSV_BDCA.exists():
        raise FileNotFoundError(f"Archivo CSV BDCA no encontrado: {CSV_BDCA}")

    df = pd.read_csv(CSV_BDCA, dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # ---------- Validación de esquema ----------
    esquemas_esperados = {"plot", "trt", "block", "yield"}
    esquemas = set(df.columns)
    faltantes = esquemas_esperados - esquemas
    extras = esquemas - esquemas_esperados

    # ---------- Conversión de tipos ----------
    df["yield"] = pd.to_numeric(df["yield"], errors="coerce")
    df["plot"] = pd.to_numeric(df["plot"], errors="coerce")
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()
    df["block"] = df["block"].astype(str).str.strip().str.upper()

    # ---------- Auditoría ----------
    total_filas = len(df)
    filas_duplicadas = df.duplicated().sum()
    filas_unicas = total_filas - filas_duplicadas

    valores_na = df.isna().sum().to_dict()
    filas_con_na = int(df.isna().any(axis=1).sum())

    # Balance: 1 obs por celda trt × block
    balance = {}
    para_balance = df[["trt", "block"]].dropna()
    if not para_balance.empty:
        conteo = para_balance.groupby(["trt", "block"]).size().reset_index(name="conteo")
        conteo_unicos = conteo[conteo["conteo"] == 1]
        balance["total_celdas"] = conteo["trt"].nunique() * conteo["block"].nunique()
        balance["celdas_balanceadas"] = len(conteo_unicos)
        balance["celdas_desbalanceadas"] = len(conteo) - balance["celdas_balanceadas"]
        balance["conteo_por_celda"] = conteo.set_index(["trt", "block"])["conteo"].to_dict()

    # Conteo por trt (R, T0, T1, T2)
    conteo_por_trt = df["trt"].value_counts().to_dict()
    # Conteo por block (B1-B9)
    conteo_por_block = df["block"].value_counts().to_dict()

    # Valores extremos (yield)
    yield_series = df["yield"].dropna()
    extremos = {
        "minimo": float(yield_series.min()) if not yield_series.empty else None,
        "maximo": float(yield_series.max()) if not yield_series.empty else None,
        "media": float(yield_series.mean()) if not yield_series.empty else None,
        "std": float(yield_series.std(ddof=1)) if len(yield_series) > 1 else None,
    }

    # Valores prohibidos (yield negativo)
    yield_negativo = int((yield_series < 0).sum())

    # Requisitos de conclusión
    conclusion = {
        "total_filas": total_filas,
        "filas_unicas": filas_unicas,
        "filas_duplicadas": int(filas_duplicadas),
        "filas_con_na": filas_con_na,
        "valores_na": valores_na,
        "balance_celulas": balance,
        "conteo_por_trt": conteo_por_trt,
        "conteo_por_block": conteo_por_block,
        "extremos_yield": extremos,
        "yield_negativo": yield_negativo,
        "auditado_ok": (
            filas_unicas == total_filas
            and filas_con_na == 0
            and (yield_negativo == 0)
            and balance.get("celdas_desbalanceadas", 0) == 0
        ),
    }

    # Guardar auditoría
    guardar_tabla(pd.DataFrame([conclusion]), "auditoria_bdca", index=False)

    # Devolver solo las columnas utilizadas en el análisis
    df_final = df[["trt", "block", "yield"]].copy()
    return df_final, conclusion


def _ejemplos_anomalias(df: pd.DataFrame, col: str, max_rows: int = 5) -> list:
    """Devuelve ejemplos (hasta max_rows) de filas con anomalías en la columna dada."""
    anomalies = []
    if col not in df.columns:
        return anomalies

    # Duplicados (mismas trt+block)
    if col == "duplicados":
        dup_mask = df.duplicated(subset=["trt", "block"], keep=False)
        for _, row in df[dup_mask].head(max_rows).iterrows():
            anomalies.append(row.to_dict())
        return anomalies

    # Valores NA
    if col == "na":
        na_mask = df["yield"].isna()
        for _, row in df[na_mask].head(max_rows).iterrows():
            anomalies.append(row.to_dict())
        return anomalies

    # Yield negativo
    if col == "negativo_yield":
        neg_mask = df["yield"] < 0
        for _, row in df[neg_mask].head(max_rows).iterrows():
            anomalies.append(row.to_dict())
        return anomalies

    # Bloques desbalanceados
    if col == "desbalanceado_celda":
        conteo = df.groupby(["trt", "block"]).size()
        desbalance = conteo[conteo != 1]
        for (trt, block), cnt in desbalance.head(max_rows).items():
            anomalies.append({"trt": trt, "block": block, "conteo": int(cnt)})
        return anomalies

    return anomalies


# Script principal para ejecución independiente
if __name__ == "__main__":
    df, audit = cargar()
    print(f"BDCA cargado: {len(df)} filas únicas.")
    print(f"Tratamientos: {sorted(df['trt'].unique())}")
    print(f"Bloques: {sorted(df['block'].unique())}")
    print(f"Balance OK: {audit['auditado_ok']}")
    if not audit["auditado_ok"]:
        print("Advertencias/anomalías detectadas (ver auditoria_bdca.csv)")
