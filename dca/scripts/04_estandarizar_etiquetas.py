"""
04_estandarizar_etiquetas.py — Normalización de etiquetas en los datasets.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import numpy as np
from config import DIR_DATOS, NORMALIZAR_AISLADOS

ARCHIVOS = [
    ("crecimiento_micelial_raw.csv", "crecimiento_micelial_estandarizado.csv"),
    ("conidias_raw.csv", "conidias_estandarizado.csv"),
]

NORMALIZAR_METODO = {
    "maceracion": "maceracion", "soxhlet": "soxhlet", "ultrasonido": "ultrasonido",
}


def estandarizar(df):
    if "metodo_extraccion" in df.columns:
        df["metodo_extraccion"] = df["metodo_extraccion"].str.strip().str.lower()
    if "aislado_id" in df.columns:
        df["aislado_id"] = df["aislado_id"].apply(
            lambda x: " ".join(str(x).strip().split())
        )
        df["aislado_id"] = df["aislado_id"].replace(NORMALIZAR_AISLADOS)
    if "grupo_experimental" in df.columns:
        df["grupo_experimental"] = df["grupo_experimental"].apply(
            lambda x: " ".join(str(x).strip().split()) if pd.notna(x) else None
        )
    sort_cols = [c for c in ["metodo_extraccion", "aislado_id", "concentracion_mg_ml", "replica_biologica"]
                 if c in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols).reset_index(drop=True)
    return df


print("📄 Estandarizando etiquetas...")

for entrada, salida in ARCHIVOS:
    ruta_entrada = DIR_DATOS / entrada
    ruta_salida = DIR_DATOS / salida
    if not ruta_entrada.exists():
        print(f"  ⚠ No encontrado: {ruta_entrada}")
        continue
    df = pd.read_csv(ruta_entrada)
    print(f"  {entrada}: {len(df)} filas → ", end="")
    df = estandarizar(df)
    df.to_csv(ruta_salida, index=False)
    print(f"{len(df)} filas → {salida}")

print("✅ Estandarización completa.")
