"""
02_extraer_crecimiento.py — Extracción de datos de crecimiento micelial.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
from config import METODOS_EXTRACCION, DIR_DATOS, parsear_crecimiento

SALIDA = DIR_DATOS / "crecimiento_micelial_raw.csv"
print("📄 Extrayendo datos de crecimiento micelial...")

todas = []
for hoja in METODOS_EXTRACCION:
    print(f"  Procesando hoja: {hoja}...")
    df_hoja = parsear_crecimiento(hoja)
    print(f"    → {len(df_hoja)} filas extraídas")
    todas.append(df_hoja)

crecimiento = pd.concat(todas, ignore_index=True)

print(f"\n📊 Total: {len(crecimiento)} filas")
print(f"  Aislados: {crecimiento['aislado_id'].nunique()}")
print(f"  Métodos: {crecimiento['metodo_extraccion'].unique()}")
print(f"  Concentraciones: {sorted(crecimiento['concentracion_mg_ml'].unique())}")

resumen = crecimiento.groupby(["metodo_extraccion", "concentracion_mg_ml"])[
    "crecimiento_mm"
].apply(lambda x: x.notna().sum())
print(resumen.to_string())

crecimiento.to_csv(SALIDA, index=False)
print(f"\n✅ Guardado: {SALIDA} ({len(crecimiento)} filas)")
