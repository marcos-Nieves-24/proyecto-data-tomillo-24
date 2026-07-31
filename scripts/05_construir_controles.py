"""
05_construir_controles.py — Extracción y desnormalización de controles.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import numpy as np
from config import DIR_DATOS

ENTRADAS = [
    ("crecimiento_micelial_estandarizado.csv", "crecimiento_micelial_con_controles.csv", "crecimiento_mm"),
    ("conidias_estandarizado.csv", "conidias_con_controles.csv", "conidias_log10"),
]


def construir_controles(df, col_respuesta):
    controles = df[df["es_control"] & df[col_respuesta].notna()][
        ["aislado_id", "metodo_extraccion", col_respuesta]
    ].drop_duplicates(subset=["aislado_id", "metodo_extraccion"])
    controles = controles.rename(columns={col_respuesta: "control_valor"})

    df_con = df.merge(controles, on=["aislado_id", "metodo_extraccion"], how="left")

    if "crecimiento_mm" in df.columns:
        df_con = df_con.rename(columns={"control_valor": "control_crecimiento_mm"})
    elif "conidias_log10" in df.columns:
        df_con = df_con.rename(columns={"control_valor": "control_conidias_log10"})
        df_con["control_conidias_crudo"] = 10 ** df_con["control_conidias_log10"]
    return df_con


print("📄 Construyendo controles desnormalizados...")

for entrada, salida, col_resp in ENTRADAS:
    ruta_entrada = DIR_DATOS / entrada
    ruta_salida = DIR_DATOS / salida
    if not ruta_entrada.exists():
        print(f"  ⚠ No encontrado: {ruta_entrada}")
        continue

    df = pd.read_csv(ruta_entrada)
    n_trat = (~df["es_control"]).sum()
    df_res = construir_controles(df, col_resp)
    col_ctrl = "control_crecimiento_mm" if col_resp == "crecimiento_mm" else "control_conidias_log10"
    n_con_ctrl = df_res[~df_res["es_control"]][col_ctrl].notna().sum()
    df_res.to_csv(ruta_salida, index=False)
    print(f"  {entrada}: {n_con_ctrl}/{n_trat} tratamientos con control → {salida}")

print("✅ Controles construidos.")
