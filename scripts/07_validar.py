"""
07_validar.py — Validación de calidad de los datos procesados.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import numpy as np
from config import DIR_DATOS, DIR_REPORTES

ARCHIVOS = [
    {
        "ruta": DIR_DATOS / "rendimiento_extraccion.csv",
        "nombre": "Rendimiento de extracción",
        "rangos": {"rendimiento_pct": (0, 100), "peso_material_seco_g": (0, 100)},
    },
    {
        "ruta": DIR_DATOS / "crecimiento_micelial_con_inhibicion.csv",
        "nombre": "Crecimiento micelial",
        "rangos": {"crecimiento_mm": (0, 90), "control_crecimiento_mm": (0, 90),
                    "porcentaje_inhibicion": (-100, 100)},
    },
    {
        "ruta": DIR_DATOS / "conidias_con_inhibicion.csv",
        "nombre": "Conidias",
        "rangos": {"conidias_log10": (0, 10), "control_conidias_log10": (0, 10),
                    "porcentaje_inhibicion": (-500, 100),   # crudo permite >100% de aumento
                    "porcentaje_inhibicion_log10": (-100, 100)},
    },
]


def validar(cfg):
    ruta, nombre = cfg["ruta"], cfg["nombre"]
    if not ruta.exists():
        print(f"  ❌ {nombre}: ARCHIVO NO ENCONTRADO")
        return

    df = pd.read_csv(ruta)
    print(f"\n  ── {nombre} ({len(df)} filas × {len(df.columns)} cols) ──")

    # Faltantes
    nulos = {c: int(df[c].isna().sum()) for c in df.columns if df[c].isna().sum() > 0}
    if nulos:
        print(f"    ⚠ Valores faltantes:")
        for c, n in sorted(nulos.items(), key=lambda x: -x[1]):
            print(f"      {c}: {n} ({100*n/len(df):.0f}%)")
    else:
        print(f"    ✅ Sin valores faltantes")

    # Rangos
    for col, (lo, hi) in cfg.get("rangos", {}).items():
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce")
            fuera = ((vals < lo) | (vals > hi)).sum()
            if fuera:
                print(f"    ❌ {col}: {fuera} valores fuera de [{lo}, {hi}]")

    # Duplicados
    dups = df.duplicated().sum()
    if dups:
        print(f"    ⚠ {dups} filas duplicadas")

    # Específicos
    if "inhibicion_completa" in df.columns:
        print(f"    🔬 Inhibición completa: {df['inhibicion_completa'].sum()}")
    if "inhibicion_negativa" in df.columns:
        n_neg = (~df["es_control"] & df["inhibicion_negativa"]).sum()
        if n_neg:
            print(f"    ⚠ Inhibición negativa: {n_neg}")


print("📄 Validando datasets...")
for cfg in ARCHIVOS:
    validar(cfg)
print("\n✅ Validación completa.")
