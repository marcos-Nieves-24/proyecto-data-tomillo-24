"""
03_extraer_conidias.py — Extracción de datos de producción de conidias.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import numpy as np
from config import METODOS_EXTRACCION, DIR_DATOS, parsear_bloque_metodo


def parsear_conidias(hoja):
    df = parsear_bloque_metodo(hoja, bloque="conidias")
    df = df.rename(columns={"crecimiento_mm": "conidias_log10"})
    df["conidias_crudo"] = 10 ** df["conidias_log10"]
    return df


SALIDA = DIR_DATOS / "conidias_raw.csv"
print("📄 Extrayendo datos de producción de conidias...")

todas = []
for hoja in METODOS_EXTRACCION:
    print(f"  Procesando hoja: {hoja}...")
    df_hoja = parsear_conidias(hoja)
    print(f"    → {len(df_hoja)} filas, conidias no nulas: {df_hoja['conidias_log10'].notna().sum()}")
    todas.append(df_hoja)

conidias = pd.concat(todas, ignore_index=True)

print(f"\n📊 Total: {len(conidias)} filas")
print(f"  Rango log10: {conidias['conidias_log10'].min():.2f} – {conidias['conidias_log10'].max():.2f}")

resumen = conidias.groupby(["metodo_extraccion", "concentracion_mg_ml"])[
    "conidias_log10"
].apply(lambda x: x.notna().sum())
print(resumen.to_string())

conidias.to_csv(SALIDA, index=False)
print(f"\n✅ Guardado: {SALIDA} ({len(conidias)} filas)")
