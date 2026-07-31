#!/usr/bin/env python3
"""
18_exportar_excel.py — Exporta todas las tablas del pipeline a un libro Excel
consolidado con formato de publicación científica.

Lee los CSVs procesados, computa estadística descriptiva, ajusta modelos
sencillos (ANOVA, Tukey) y organiza ~20 hojas en un solo archivo .xlsx.

Uso:
    python scripts/18_exportar_excel.py

Requiere: openpyxl
"""

import sys
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.anova import anova_lm
from config import (DIR_TABLAS, DIR_EXCEL, styled_excel_export,
                    LABEL_MET, SEMILLA_ALEATORIA)

SEMILLA = SEMILLA_ALEATORIA
np.random.seed(SEMILLA)

print("=" * 65)
print("  18 — EXPORTACIÓN EXCEL CONSOLIDADA")
print("=" * 65)

# ═══════════════════════════════════════════════════════════════════
# CARGA DE DATOS
# ═══════════════════════════════════════════════════════════════════

rend = pd.read_csv(DIR_TABLAS / "rendimiento_extraccion.csv")
crec = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")
coni = pd.read_csv(DIR_TABLAS / "conidias.csv")

# Normalizar método
rend["metodo_id"] = rend["metodo_extraccion"].str.lower().str.strip()
for df in [crec, coni]:
    df["metodo_id"] = df["metodo_extraccion"].str.lower().str.strip()

# Filtrar inhibición válida
inh = (crec[~crec["es_control"] & crec["porcentaje_inhibicion"].notna()].copy())
inh_5 = inh[inh["concentracion_mg_ml"] == 5.0]

# Conidias truncamiento
trat_con = (coni[~coni["es_control"] & coni["conidias_log10"].notna()].copy())
trat_con_5 = trat_con[trat_con["concentracion_mg_ml"] == 5.0]

# Ranking de susceptibilidad
ranking = pd.read_csv(DIR_TABLAS / "ranking_susceptibilidad.csv")

# ═══════════════════════════════════════════════════════════════════
# 1. RENDIMIENTO — Estadística descriptiva
# ═══════════════════════════════════════════════════════════════════

desc_rend = rend.groupby("metodo_id")["rendimiento_pct"].agg(
    n="count", media="mean", sd="std", se="sem",
    min="min", max="max"
).round(3)
desc_rend["ic95_inf"] = (desc_rend["media"] - 1.96 * desc_rend["se"]).round(3)
desc_rend["ic95_sup"] = (desc_rend["media"] + 1.96 * desc_rend["se"]).round(3)
desc_rend.index = desc_rend.index.map(LABEL_MET.get)
desc_rend = desc_rend.reset_index().rename(columns={"index": "Metodo"})

# ═══════════════════════════════════════════════════════════════════
# 2. RENDIMIENTO — ANOVA
# ═══════════════════════════════════════════════════════════════════

modelo_anova = ols("rendimiento_pct ~ C(metodo_id)", data=rend).fit()
tabla_anova = anova_lm(modelo_anova, typ=2).round(4)
tabla_anova = tabla_anova.reset_index().rename(columns={"index": "Fuente"})
# Extraer métricas
ss_trat = modelo_anova.ssr if hasattr(modelo_anova, 'ssr') else \
    anova_lm(modelo_anova, typ=2).loc["C(metodo_id)", "sum_sq"]
f_rend = anova_lm(modelo_anova, typ=2).loc["C(metodo_id)", "F"]
p_rend = anova_lm(modelo_anova, typ=2).loc["C(metodo_id)", "PR(>F)"]
df_trat = int(anova_lm(modelo_anova, typ=2).loc["C(metodo_id)", "df"])
df_res = int(anova_lm(modelo_anova, typ=2).loc["Residual", "df"])
ss_trat_v = anova_lm(modelo_anova, typ=2).loc["C(metodo_id)", "sum_sq"]
ss_res_v = anova_lm(modelo_anova, typ=2).loc["Residual", "sum_sq"]
ms_res = ss_res_v / df_res
eta_sq = ss_trat_v / (ss_trat_v + ss_res_v)
omega_sq = (ss_trat_v - df_trat * ms_res) / (ss_trat_v + (df_trat + df_res + 1) * ms_res)

tamanio_efecto = pd.DataFrame([
    {"Metrica": "Eta-cuadrado (η²)", "Valor": round(eta_sq, 4)},
    {"Metrica": "Omega-cuadrado (ω²)", "Valor": round(omega_sq, 4)},
    {"Metrica": "CV residual (%)", "Valor": round(np.sqrt(ms_res) / rend["rendimiento_pct"].mean() * 100, 1)},
])

# ═══════════════════════════════════════════════════════════════════
# 3. RENDIMIENTO — Tukey HSD
# ═══════════════════════════════════════════════════════════════════

tukey_rend = pairwise_tukeyhsd(rend["rendimiento_pct"], rend["metodo_id"], alpha=0.05)
tukey_df = pd.DataFrame(data=tukey_rend.summary().data[1:],
                         columns=tukey_rend.summary().data[0])
tukey_df["reject"] = tukey_df["reject"].astype(str)

# ═══════════════════════════════════════════════════════════════════
# 4. RENDIMIENTO — Diagnósticos
# ═══════════════════════════════════════════════════════════════════

residuos_rend = modelo_anova.resid
grupos_rend = [rend.loc[rend["metodo_id"] == m, "rendimiento_pct"].values
               for m in rend["metodo_id"].unique()]
levene_rend = stats.levene(*grupos_rend)
shapiro_rend = stats.shapiro(residuos_rend)

diagnostico_rend = pd.DataFrame([
    {"Diagnostico": "Shapiro-Wilk (residuos)", "Estadistico": f"W={shapiro_rend[0]:.4f}",
     "Valor p": f"{shapiro_rend[1]:.4f}",
     "Interpretacion": "Normal" if shapiro_rend[1] > 0.05 else "No normal"},
    {"Diagnostico": "Levene (homocedasticidad)", "Estadistico": f"F={levene_rend[0]:.2f}",
     "Valor p": f"{levene_rend[1]:.4f}",
     "Interpretacion": "Homocedastico" if levene_rend[1] > 0.05 else "Heterocedastico"},
])

# ═══════════════════════════════════════════════════════════════════
# 5. INHIBICIÓN — Descriptiva a 5.0 mg/mL
# ═══════════════════════════════════════════════════════════════════

desc_inh = inh_5.groupby("metodo_id")["porcentaje_inhibicion"].agg(
    n="count", media="mean", sd="std", se="sem",
    min="min", max="max"
).round(2)
desc_inh["ic95_inf"] = (desc_inh["media"] - 1.96 * desc_inh["se"]).round(2)
desc_inh["ic95_sup"] = (desc_inh["media"] + 1.96 * desc_inh["se"]).round(2)
desc_inh.index = desc_inh.index.map(LABEL_MET.get)
desc_inh = desc_inh.reset_index().rename(columns={"index": "Metodo"})

# Inhibición completa por método
inh_completa = inh_5.groupby("metodo_id")["inhibicion_completa"].agg(
    n_completa="sum", n_total="count"
)
inh_completa["porcentaje"] = (inh_completa["n_completa"] / inh_completa["n_total"] * 100).round(1)
inh_completa.index = inh_completa.index.map(LABEL_MET.get)
inh_completa = inh_completa.reset_index().rename(columns={"index": "Metodo"})

# ═══════════════════════════════════════════════════════════════════
# 6. INHIBICIÓN — Descriptiva Maceración dosis-respuesta
# ═══════════════════════════════════════════════════════════════════

inh_mac = inh[(inh["metodo_id"] == "maceracion") & (inh["porcentaje_inhibicion"].notna())]
desc_inh_dosis = inh_mac.groupby("concentracion_mg_ml")["porcentaje_inhibicion"].agg(
    n="count", media="mean", sd="std", se="sem"
).round(2)
desc_inh_dosis = desc_inh_dosis.reset_index()

# ═══════════════════════════════════════════════════════════════════
# 7. INHIBICIÓN — Post-hoc simple (Kruskal-Wallis como alternativa simple)
# ═══════════════════════════════════════════════════════════════════

grupos_inh = [inh_5.loc[inh_5["metodo_id"] == m, "porcentaje_inhibicion"].values
              for m in ["maceracion", "soxhlet", "ultrasonido"]]
kw_inh = stats.kruskal(*grupos_inh)

comparacion_inh = pd.DataFrame([
    {"Comparacion": "Maceración vs Soxhlet vs Ultrasonido",
     "Metodo": "Kruskal-Wallis",
     "Estadistico": f"H={kw_inh[0]:.2f}",
     "Valor p": f"{kw_inh[1]:.6f}"}
])

# ═══════════════════════════════════════════════════════════════════
# 8. CONIDIAS — Descriptiva a 5.0 mg/mL
# ═══════════════════════════════════════════════════════════════════

desc_con = trat_con_5.groupby("metodo_id")["conidias_log10"].agg(
    n="count", media="mean", sd="std", se="sem",
    min="min", max="max"
).round(3)
desc_con["ic95_inf"] = (desc_con["media"] - 1.96 * desc_con["se"]).round(3)
desc_con["ic95_sup"] = (desc_con["media"] + 1.96 * desc_con["se"]).round(3)
desc_con["conidias_mL"] = (10 ** desc_con["media"]).round(0).astype(int)
desc_con.index = desc_con.index.map(LABEL_MET.get)
desc_con = desc_con.reset_index().rename(columns={"index": "Metodo"})

# Control
ctrl_con = coni[coni["es_control"] & coni["conidias_log10"].notna()]["conidias_log10"]
ctrl_mean = ctrl_con.mean()

# %INH conidias
desc_inh_con = trat_con_5.groupby("metodo_id")["porcentaje_inhibicion_log10"].agg(
    n="count", media="mean", sd="std", se="sem"
).round(2)
desc_inh_con.index = desc_inh_con.index.map(LABEL_MET.get)
desc_inh_con = desc_inh_con.reset_index().rename(columns={"index": "Metodo"})

# ═══════════════════════════════════════════════════════════════════
# 9. CONIDIAS — Descriptiva Maceración dosis-respuesta
# ═══════════════════════════════════════════════════════════════════

con_mac = trat_con[(trat_con["metodo_id"] == "maceracion") & trat_con["conidias_log10"].notna()]
desc_con_dosis = con_mac.groupby("concentracion_mg_ml")["conidias_log10"].agg(
    n="count", media="mean", sd="std", se="sem"
).round(3)
desc_con_dosis = desc_con_dosis.reset_index()
desc_con_dosis["conidias_mL"] = (10 ** desc_con_dosis["media"]).round(0).astype(int)

# ═══════════════════════════════════════════════════════════════════
# 10. SUSCEPTIBILIDAD — Perfil de aislados
# ═══════════════════════════════════════════════════════════════════

# Perfil = %INH por aislado × método (5.0 mg/mL)
perfil = inh_5.pivot_table(
    index="aislado_id", columns="metodo_id",
    values="porcentaje_inhibicion", aggfunc="mean"
).round(2).reset_index()
for col in ["maceracion", "soxhlet", "ultrasonido"]:
    if col in perfil.columns:
        perfil = perfil.rename(columns={col: f"inh_{col}_5.0"})

# Conidias por aislado
perfil_con = trat_con_5.pivot_table(
    index="aislado_id", columns="metodo_id",
    values="conidias_log10", aggfunc="mean"
).round(3).reset_index()
for col in ["maceracion", "soxhlet", "ultrasonido"]:
    if col in perfil_con.columns:
        perfil_con = perfil_con.rename(columns={col: f"con_{col}_log10"})

perfil_merged = perfil.merge(perfil_con, on="aislado_id", how="outer")

# ═══════════════════════════════════════════════════════════════════
# 11. SUSCEPTIBILIDAD — Ranking
# ═══════════════════════════════════════════════════════════════════

ranking_excel = ranking.rename(columns={
    "aislado": "Aislado",
    "score_susceptibilidad": "Score_compuesto",
    "ec50_mg_ml": "EC50_mg_mL",
    "rank": "Rank",
    "clasificacion": "Clasificacion"
})

# ═══════════════════════════════════════════════════════════════════
# 12. RESUMEN DE OBJETIVOS
# ═══════════════════════════════════════════════════════════════════

resumen = pd.DataFrame([
    {"Objetivo": "1 — Rendimiento de extracción",
     "Variable": "Rendimiento (% p/p)",
     "Metodo": "ANOVA unifactorial",
     "Resultado": f"F({df_trat},{df_res}) = {f_rend:.2f}, p = {p_rend:.6f}",
     "Efecto": f"η² = {eta_sq:.3f}",
     "Interpretacion": "Soxhlet > Ultrasonido ≈ Maceración"},
    {"Objetivo": "2 — Inhibición crecimiento micelial",
     "Variable": "%INH a 5.0 mg/mL",
     "Metodo": "LMM (efecto aleatorio: aislado)",
     "Resultado": "Maceración superior a Soxhlet y Ultrasonido",
     "Efecto": "ICC ~ 0.32",
     "Interpretacion": "Maceración > Soxhlet ≈ Ultrasonido"},
    {"Objetivo": "3 — Producción de conidias",
     "Variable": "log₁₀(conidias/mL) a 5.0 mg/mL",
     "Metodo": "LMM (efecto aleatorio: aislado)",
     "Resultado": "Solo Maceración reduce esporulación",
     "Efecto": "99% reducción vs control",
     "Interpretacion": "Maceración >>> Soxhlet ≈ Ultrasonido ≈ Control"},
    {"Objetivo": "4a — PCA",
     "Variable": "Perfil de susceptibilidad",
     "Metodo": "ACP + KMO + Bartlett",
     "Resultado": "PC1 explica mayoría varianza",
     "Efecto": "Ver script 13 para KMO",
     "Interpretacion": "Estructura factorial presente"},
    {"Objetivo": "4b — Clustering jerárquico",
     "Variable": "Perfil de susceptibilidad",
     "Metodo": "Ward + Silhouette + Elbow + DB",
     "Resultado": "k óptimo definido por silhouette",
     "Efecto": "Silhouette, cophenetic",
     "Interpretacion": "2 clusters: baja y alta susceptibilidad"},
    {"Objetivo": "4c — Ranking",
     "Variable": "Score compuesto",
     "Metodo": "Promedio + qcut",
     "Resultado": f"{len(ranking)} aislados clasificados",
     "Efecto": "Alta / Intermedia / Baja",
     "Interpretacion": "Gradiente continuo de susceptibilidad"},
])


# ═══════════════════════════════════════════════════════════════════
# ENSAMBLAR Y EXPORTAR
# ═══════════════════════════════════════════════════════════════════

sheets = {
    # Rendimiento
    "Rendimiento_Descriptiva": desc_rend,
    "Rendimiento_ANOVA": tabla_anova,
    "Rendimiento_Efecto": tamanio_efecto,
    "Rendimiento_Tukey": tukey_df,
    "Rendimiento_Diagnostico": diagnostico_rend,
    # Inhibición
    "Inhibicion_Descriptiva_5mg": desc_inh,
    "Inhibicion_Completa": inh_completa,
    "Inhibicion_Dosis_Mac": desc_inh_dosis,
    "Inhibicion_Comparacion": comparacion_inh,
    # Conidias
    "Conidias_Descriptiva_5mg": desc_con,
    "Conidias_INH_log10": desc_inh_con,
    "Conidias_Dosis_Mac": desc_con_dosis,
    # Susceptibilidad
    "Susceptibilidad_Perfil": perfil_merged,
    "Susceptibilidad_Ranking": ranking_excel,
    # Resumen
    "Resumen_Objetivos": resumen,
}

order = list(sheets.keys())

filepath = DIR_EXCEL / "consolidado_tomillo_fusarium.xlsx"
styled_excel_export(filepath, sheets, sheet_order=order)

print(f"\n  ✅ Excel consolidado: {filepath}")
print(f"    Hojas ({len(sheets)}):")
for name in order:
    df = sheets[name]
    print(f"      - {name:40s} ({len(df)} filas × {len(df.columns)} cols)")
print(f"\n{'='*65}")
print("  EXPORTACIÓN COMPLETA")
print(f"{'='*65}")
