"""
01_extraer_rendimiento.py — Extracción de datos de rendimiento.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import numpy as np
from config import EXCEL_ORIGINAL, DIR_DATOS

HOJA = "RENDIMIENTOS"
SALIDA = DIR_DATOS / "rendimiento_extraccion.csv"

print(f"📄 Extrayendo datos de rendimiento desde hoja '{HOJA}'...")

df_raw = pd.read_excel(EXCEL_ORIGINAL, sheet_name=HOJA, engine="openpyxl")

print(f"  Filas leídas: {len(df_raw)}")

rendimiento = pd.DataFrame({
    "metodo_extraccion": df_raw["TOMILLO - METANOL"].str.strip().str.lower(),
    "peso_material_seco_g": df_raw["PESO MATERIAL SECO (g)"],
    "peso_extracto_obtenido_g": df_raw["PESO EXTRACTO OBTENIDO (g)"],
    "rendimiento_pct": df_raw["RENDIMIENTO"],
})
rendimiento["replica_biologica"] = rendimiento.groupby("metodo_extraccion").cumcount() + 1

print(f"\n📊 Resumen por método:")
print(rendimiento.groupby("metodo_extraccion")["rendimiento_pct"].describe().to_string())

assert (rendimiento["rendimiento_pct"] >= 0).all(), "¡Rendimiento negativo!"
assert (rendimiento["peso_material_seco_g"] > 0).all(), "¡Peso <= 0!"

rendimiento.to_csv(SALIDA, index=False)
print(f"\n✅ Guardado: {SALIDA} ({len(rendimiento)} filas)")
print(rendimiento.to_string())
