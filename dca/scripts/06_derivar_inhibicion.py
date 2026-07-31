"""
06_derivar_inhibicion.py — Cálculo de % de inhibición y validación.

Para crecimiento micelial:
  %INH = (1 − crecimiento_tratamiento / crecimiento_control) × 100

Para conidias:
  %INH se calcula sobre la ESCALA CRUDA (conteos), NO sobre log₁₀.
  Esto da una reducción porcentual biológicamente interpretable.
  %INH_crudo = (1 − 10^log10_tratamiento / 10^log10_control) × 100

  Se conserva %INH_log10 como referencia (el cálculo de la hoja original).
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import numpy as np
from config import DIR_DATOS

print("📄 Derivando porcentajes de inhibición...")

# ─── 1. CRECIMIENTO MICELIAL ────────────────────────────────────────
print("\n  Crecimiento micelial (%INH en escala original)...")
ruta_crec = DIR_DATOS / "crecimiento_micelial_con_controles.csv"
salida_crec = DIR_DATOS / "crecimiento_micelial_con_inhibicion.csv"

if not ruta_crec.exists():
    print(f"    ⚠ No encontrado: {ruta_crec}")
else:
    df = pd.read_csv(ruta_crec)
    df["porcentaje_inhibicion"] = np.nan
    df["inhibicion_completa"] = False
    df["inhibicion_negativa"] = False
    df["diferencia_con_hoja"] = np.nan
    df["logit_inhibicion"] = np.nan

    mask = (
        ~df["es_control"]
        & df["crecimiento_mm"].notna()
        & df["control_crecimiento_mm"].notna()
        & (df["control_crecimiento_mm"] > 0)
    )
    df.loc[mask, "porcentaje_inhibicion"] = (
        100 * (1 - df.loc[mask, "crecimiento_mm"] / df.loc[mask, "control_crecimiento_mm"])
    )
    df["inhibicion_completa"] = ~df["es_control"] & (df["crecimiento_mm"] == 0) & df["crecimiento_mm"].notna()
    # Marcar negativos ANTES de clipear (trazabilidad)
    df["inhibicion_negativa"] = ~df["es_control"] & df["porcentaje_inhibicion"].notna() & (df["porcentaje_inhibicion"] < 0)
    # Clipear a 0: biológicamente, un tratamiento no puede tener inhibición negativa
    df.loc[mask, "porcentaje_inhibicion"] = df.loc[mask, "porcentaje_inhibicion"].clip(lower=0)

    if "porcentaje_inhibicion_hoja" in df.columns:
        mask_h = df["porcentaje_inhibicion"].notna() & df["porcentaje_inhibicion_hoja"].notna()
        df["diferencia_con_hoja"] = np.where(
            mask_h, df["porcentaje_inhibicion"] - df["porcentaje_inhibicion_hoja"], np.nan
        )

    mask_l = df["porcentaje_inhibicion"].notna() & (df["porcentaje_inhibicion"] > 0) & (df["porcentaje_inhibicion"] < 100)
    p = df.loc[mask_l, "porcentaje_inhibicion"] / 100
    df.loc[mask_l, "logit_inhibicion"] = np.log(p / (1 - p))

    df.to_csv(salida_crec, index=False)
    print(f"    completa={df['inhibicion_completa'].sum()}, "
          f"negativa={df['inhibicion_negativa'].sum()}")
    print(f"    → {salida_crec.name}")

# ─── 2. CONIDIAS — sobre ESCALA CRUDA ───────────────────────────────
print("\n  Conidias (%INH sobre escala CRUDA de conteos)...")
ruta_coni = DIR_DATOS / "conidias_con_controles.csv"
salida_coni = DIR_DATOS / "conidias_con_inhibicion.csv"

if not ruta_coni.exists():
    print(f"    ⚠ No encontrado: {ruta_coni}")
else:
    df = pd.read_csv(ruta_coni)
    df["porcentaje_inhibicion"] = np.nan           # %INH sobre CRUDO (primario)
    df["porcentaje_inhibicion_log10"] = np.nan     # %INH sobre LOG10 (referencia)
    df["inhibicion_completa"] = False
    df["inhibicion_negativa"] = False
    df["diferencia_con_hoja"] = np.nan
    df["logit_inhibicion"] = np.nan

    # Asegurar que columnas crudas existen
    if "conidias_crudo" not in df.columns:
        df["conidias_crudo"] = 10 ** df["conidias_log10"]
    if "control_conidias_crudo" not in df.columns:
        df["control_conidias_crudo"] = 10 ** df["control_conidias_log10"]

    mask = (
        ~df["es_control"]
        & df["conidias_crudo"].notna()
        & df["control_conidias_crudo"].notna()
        & (df["control_conidias_crudo"] > 0)
    )
    # %INH sobre escala cruda (biológicamente interpretable)
    df.loc[mask, "porcentaje_inhibicion"] = (
        100 * (1 - df.loc[mask, "conidias_crudo"] / df.loc[mask, "control_conidias_crudo"])
    )
    # %INH sobre escala log10 (referencia / hoja original)
    mask_log = (
        ~df["es_control"]
        & df["conidias_log10"].notna()
        & df["control_conidias_log10"].notna()
        & (df["control_conidias_log10"] > 0)
    )
    df.loc[mask_log, "porcentaje_inhibicion_log10"] = (
        100 * (1 - df.loc[mask_log, "conidias_log10"] / df.loc[mask_log, "control_conidias_log10"])
    )

    # Flags (ANTES de clipear — trazabilidad)
    df["inhibicion_completa"] = mask & (df["conidias_crudo"] == 0)
    df["inhibicion_negativa"] = mask & (df["porcentaje_inhibicion"] < 0)

    # Clipear a 0: un tratamiento antifúngico no puede tener inhibición negativa
    df.loc[mask, "porcentaje_inhibicion"] = df.loc[mask, "porcentaje_inhibicion"].clip(lower=0)
    df.loc[mask_log, "porcentaje_inhibicion_log10"] = df.loc[mask_log, "porcentaje_inhibicion_log10"].clip(lower=0)

    # Diferencia con hoja (solo para C1 = 5 mg/mL, igual que en crecimiento)
    if "porcentaje_inhibicion_hoja" in df.columns:
        # La hoja calculó %INH sobre log10
        mask_h = df["porcentaje_inhibicion"].notna() & df["porcentaje_inhibicion_hoja"].notna()
        df["diferencia_con_hoja"] = np.where(
            mask_h, df["porcentaje_inhibicion"] - df["porcentaje_inhibicion_hoja"], np.nan
        )

    # Logit sobre %INH crudo
    mask_l = df["porcentaje_inhibicion"].notna() & (df["porcentaje_inhibicion"] > 0) & (df["porcentaje_inhibicion"] < 100)
    p = df.loc[mask_l, "porcentaje_inhibicion"] / 100
    df.loc[mask_l, "logit_inhibicion"] = np.log(p / (1 - p))

    df.to_csv(salida_coni, index=False)
    print(f"    completa={df['inhibicion_completa'].sum()}, "
          f"negativa={df['inhibicion_negativa'].sum()}")
    print(f"    → {salida_coni.name}")

    # Mostrar diferencias entre escalas
    comp = df[mask & df["porcentaje_inhibicion_log10"].notna()][
        ["aislado_id", "metodo_extraccion", "concentracion_mg_ml", "replica_biologica",
         "conidias_crudo", "control_conidias_crudo",
         "porcentaje_inhibicion", "porcentaje_inhibicion_log10"]
    ].copy()
    comp["diferencia"] = comp["porcentaje_inhibicion"] - comp["porcentaje_inhibicion_log10"]
    print(f"\n    Comparación %INH crudo vs log10 (primeros 5 casos):")
    print(f"    {'Aislado':15s} {'Conc':6s} {'Crudo':>7s} {'Control':>8s} "
          f"{'%INH_crudo':>10s} {'%INH_log10':>10s} {'Diff':>7s}")
    for _, row in comp.head(5).iterrows():
        print(f"    {row['aislado_id']:15s} {row['concentracion_mg_ml']:5.1f} "
              f"{row['conidias_crudo']:>7.0f} {row['control_conidias_crudo']:>8.0f} "
              f"{row['porcentaje_inhibicion']:>10.2f} {row['porcentaje_inhibicion_log10']:>10.2f} "
              f"{row['diferencia']:>7.2f}")
    print(f"\n    Media diff crudo − log10: {comp['diferencia'].mean():.2f}")

print("\n✅ Inhibición derivada.")
