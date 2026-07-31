#!/usr/bin/env python3
"""
12_modelo_conidias.py — Objetivo 3: Efecto del tratamiento sobre la
producción de conidias.

Estrategia:
  A) LMM sobre log10(conidias/mL) a 5.0 mg/mL — compara los 3 métodos
  B) LMM sobre log10(conidias) — Maceración dosis-respuesta
  C) LMM sobre %INH crudo de conidias a 5.0 mg/mL — comparable a Objetivo 2

Los datos de conidias fueron reportados por el laboratorio en escala
log10. Los conteos crudos se derivan por back-transformación.
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
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools import add_constant
import statsmodels.api as sm
from config import (DIR_TABLAS, DIR_REPORTES, DIR_FIGURAS, COLOR_MET,
                    LABEL_MET, setup_figure_style, save_figure_pub,
                    diagnostic_durbin_watson, diagnostic_breusch_pagan,
                    interpretar_dw)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

SEMILLA = 42
np.random.seed(SEMILLA)
setup_figure_style()

# ═══════════════════════════════════════════════════════════════════
print("=" * 65)
print("  OBJETIVO 3 — PRODUCCIÓN DE CONIDIAS")
print("=" * 65)

# ─── 1. Cargar datos ─────────────────────────────────────────────
coni = pd.read_csv(DIR_TABLAS / "conidias.csv")
trat = (coni[~coni["es_control"] & coni["conidias_log10"].notna()]
        .copy())
trat["metodo_id"] = trat["metodo_extraccion"].map(
    {"maceracion": "maceracion", "maceración": "maceracion",
     "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"}
)
trat["conc_cat"] = trat["concentracion_mg_ml"].astype("category")
trat["conc_log"] = np.log(trat["concentracion_mg_ml"] + 0.01)

print(f"\n  Total obs con conidias: {len(trat)}")
print(f"  Aislados: {trat['aislado_id'].nunique()}")
print(f"  Rango log10(conidias): [{trat['conidias_log10'].min():.2f}, {trat['conidias_log10'].max():.2f}]")

# Control pool para referencia
ctrl = coni[coni["es_control"] & coni["conidias_log10"].notna()]
print(f"  Controles: {len(ctrl)} obs, media log10={ctrl['conidias_log10'].mean():.2f}")

# ─── 2. MODELO A: log10 conidias a 5.0 mg/mL ─────────────────────
print("\n" + "─" * 65)
print("  MODELO A: log10(conidias) a 5.0 mg/mL — comparación de métodos")
print("─" * 65)

df_a = trat[trat["concentracion_mg_ml"] == 5.0].copy()
print(f"  Datos: {len(df_a)} obs, {df_a['aislado_id'].nunique()} aislados")
print(f"  Media control (global): {ctrl['conidias_log10'].mean():.2f}")

try:
    df_a_ml = pd.get_dummies(df_a, columns=["metodo_id"], drop_first=True, dtype=float)
    groups_a = df_a_ml["aislado_id"]
    exog_a = add_constant(df_a_ml[["metodo_id_soxhlet", "metodo_id_ultrasonido"]])
    endog_a = df_a_ml["conidias_log10"]

    modelo_a = MixedLM(endog_a, exog_a, groups=groups_a).fit(reml=True, maxiter=200)
    print(f"\n  Convergió: {modelo_a.converged}  |  Log-Lik: {modelo_a.llf:.1f}")
    print(f"\n  Efectos fijos:\n{modelo_a.fe_params.to_string()}")

    coefs_a = pd.DataFrame({
        "Coef": modelo_a.fe_params, "EE": modelo_a.bse_fe,
        "z": modelo_a.tvalues, "p_valor": modelo_a.pvalues,
    })
    coefs_a["IC95_inf"] = coefs_a["Coef"] - 1.96 * coefs_a["EE"]
    coefs_a["IC95_sup"] = coefs_a["Coef"] + 1.96 * coefs_a["EE"]
    print(f"\n{coefs_a.round(3).to_string()}")

    var_iso = modelo_a.cov_re.iloc[0, 0]
    var_res = modelo_a.scale
    icc_a = var_iso / (var_iso + var_res)
    print(f"\n  Var(aislado) = {var_iso:.3f}")
    print(f"  Var(residual) = {var_res:.3f}")
    print(f"  ICC (aislado) = {icc_a:.3f}")

    # Medias marginales (en escala log10)
    pred_mac_a = modelo_a.fe_params["const"]
    pred_sox_a = modelo_a.fe_params["const"] + modelo_a.fe_params["metodo_id_soxhlet"]
    pred_ult_a = modelo_a.fe_params["const"] + modelo_a.fe_params["metodo_id_ultrasonido"]
    se_mac_a = modelo_a.bse_fe["const"]
    for lbl, est, se in [("Maceración", pred_mac_a, se_mac_a),
                          ("Soxhlet", pred_sox_a, modelo_a.bse_fe["metodo_id_soxhlet"]),
                          ("Ultrasonido", pred_ult_a, modelo_a.bse_fe["metodo_id_ultrasonido"])]:
        # SE around estimate (for methods vs ref, use difference SE)
        print(f"    {lbl:15s} log10 = {est:.2f}  ({10**est:.0f} conidias/mL)")

    # Back-transform to % of control
    ctrl_mean = ctrl["conidias_log10"].mean()
    print(f"\n  Reducción vs control (control log10 = {ctrl_mean:.2f}):")
    for lbl, est in [("Maceración", pred_mac_a), ("Soxhlet", pred_sox_a),
                      ("Ultrasonido", pred_ult_a)]:
        reduccion_log = ctrl_mean - est
        reduccion_pct = (1 - 10**(-reduccion_log)) * 100
        print(f"    {lbl:15s}: -{reduccion_log:.2f} log10 = {reduccion_pct:.1f}% reducción")

    modelo_a_ok = True

    # Diagnóstico de residuos — Modelo A
    dw_a = diagnostic_durbin_watson(modelo_a.resid)
    print(f"  Durbin-Watson: {dw_a:.3f} — {interpretar_dw(dw_a)}")
    bp_a = diagnostic_breusch_pagan(modelo_a, exog_a)
    print(f"  Breusch-Pagan: LM={bp_a['lm_stat']:.2f}, p={bp_a['p_val']:.4f} "
          f"{'⚠ Heterocedástico' if bp_a['p_val'] < 0.05 else '✅ Homocedástico'}")
except Exception as e:
    print(f"  ⚠ Modelo A falló: {e}")
    import traceback; traceback.print_exc()
    modelo_a_ok = False


# ─── 3. MODELO B: Maceración dosis-respuesta (log10 conidias) ────
print("\n" + "─" * 65)
print("  MODELO B: Maceración — dosis-respuesta en log10(conidias)")
print("─" * 65)

df_b = trat[trat["metodo_id"] == "maceracion"].copy()
print(f"  Datos: {len(df_b)} obs")
print(f"  Concentraciones: {sorted(df_b['concentracion_mg_ml'].unique())}")

try:
    dummies_conc = pd.get_dummies(df_b["conc_cat"], drop_first=True, dtype=float, prefix="conc")
    df_b_ml = pd.concat([df_b.reset_index(drop=True), dummies_conc.reset_index(drop=True)], axis=1)
    groups_b = df_b_ml["aislado_id"]
    exog_b = add_constant(df_b_ml[[c for c in dummies_conc.columns]])
    endog_b = df_b_ml["conidias_log10"]

    modelo_b = MixedLM(endog_b, exog_b, groups=groups_b).fit(reml=True, maxiter=200)
    print(f"\n  Convergió: {modelo_b.converged}  |  Log-Lik: {modelo_b.llf:.1f}")

    coefs_b = pd.DataFrame({
        "Coef": modelo_b.fe_params, "EE": modelo_b.bse_fe,
        "z": modelo_b.tvalues, "p_valor": modelo_b.pvalues,
    })
    coefs_b["IC95_inf"] = coefs_b["Coef"] - 1.96 * coefs_b["EE"]
    coefs_b["IC95_sup"] = coefs_b["Coef"] + 1.96 * coefs_b["EE"]
    print(f"\n{coefs_b.round(3).to_string()}")

    var_iso_b = modelo_b.cov_re.iloc[0, 0]
    var_res_b = modelo_b.scale
    icc_b = var_iso_b / (var_iso_b + var_res_b)
    print(f"\n  Var(aislado) = {var_iso_b:.3f}")
    print(f"  Var(residual) = {var_res_b:.3f}")
    print(f"  ICC (aislado) = {icc_b:.3f}")

    conc_ref = sorted(df_b['concentracion_mg_ml'].unique())[0]
    const_b = modelo_b.fe_params["const"]
    print(f"\n  Medias marginales (ref = {conc_ref} mg/mL):")
    for conc in sorted(df_b['concentracion_mg_ml'].unique()):
        if conc == conc_ref:
            est = const_b
        else:
            col = f"conc_{conc}"
            est = const_b + modelo_b.fe_params[col]
        print(f"    {conc:>5} mg/mL: log10 = {est:.2f}  ({10**est:.0f} conidias/mL)")
    modelo_b_ok = True

    # Diagnóstico de residuos — Modelo B
    dw_b = diagnostic_durbin_watson(modelo_b.resid)
    print(f"  Durbin-Watson: {dw_b:.3f} — {interpretar_dw(dw_b)}")
except Exception as e:
    print(f"  ⚠ Modelo B falló: {e}")
    import traceback; traceback.print_exc()
    modelo_b_ok = False


# ─── 4. MODELO C: %INH conidias (escala log₁₀ de la hoja) ─────────
print("\n" + "─" * 65)
print("  MODELO C: %INH conidias (escala log₁₀, reportada por el laboratorio)")
print("─" * 65)

df_c = trat[(trat["concentracion_mg_ml"] == 5.0)
            & trat["porcentaje_inhibicion_log10"].notna()].copy()
print(f"  Datos: {len(df_c)} obs")
print("  NOTA: Se usa %INH_log10 en lugar de %INH crudo porque la escala")
print("  cruda tiene skew=-12.2 y kurtosis=176 (no transformable).")
print("  La escala log10 (de la hoja Excel original) es más estable.")

try:
    df_c_ml = pd.get_dummies(df_c, columns=["metodo_id"], drop_first=True, dtype=float)
    groups_c = df_c_ml["aislado_id"]
    exog_c = add_constant(df_c_ml[["metodo_id_soxhlet", "metodo_id_ultrasonido"]])
    endog_c = df_c_ml["porcentaje_inhibicion_log10"]

    modelo_c = MixedLM(endog_c, exog_c, groups=groups_c).fit(reml=True, maxiter=200)
    print(f"\n  Convergió: {modelo_c.converged}  |  Log-Lik: {modelo_c.llf:.1f}")

    coefs_c = pd.DataFrame({
        "Coef": modelo_c.fe_params, "EE": modelo_c.bse_fe,
        "z": modelo_c.tvalues, "p_valor": modelo_c.pvalues,
    })
    coefs_c["IC95_inf"] = coefs_c["Coef"] - 1.96 * coefs_c["EE"]
    coefs_c["IC95_sup"] = coefs_c["Coef"] + 1.96 * coefs_c["EE"]
    print(f"\n{coefs_c.round(3).to_string()}")

    pred_mac_c = modelo_c.fe_params["const"]
    pred_sox_c = modelo_c.fe_params["const"] + modelo_c.fe_params["metodo_id_soxhlet"]
    pred_ult_c = modelo_c.fe_params["const"] + modelo_c.fe_params["metodo_id_ultrasonido"]
    print(f"\n  %INH estimado (escala log₁₀):")
    for lbl, est in [("Maceración", pred_mac_c), ("Soxhlet", pred_sox_c),
                      ("Ultrasonido", pred_ult_c)]:
        print(f"    {lbl:15s}: {est:.1f}%")
    modelo_c_ok = True
except Exception as e:
    print(f"  ⚠ Modelo C falló: {e}")
    import traceback; traceback.print_exc()
    modelo_c_ok = False


# ═══════════════════════════════════════════════════════════════════
# 5. FIGURAS
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  FIGURAS")
print("─" * 65)

# 5a. Boxplot log10(conidias) a 5.0 mg/mL
fig, ax = plt.subplots(figsize=(8, 5))
order = ["maceracion", "soxhlet", "ultrasonido"]
sns.boxplot(data=df_a, x="metodo_id", y="conidias_log10", hue="metodo_id",
            palette=COLOR_MET, ax=ax, order=order, legend=False)
sns.stripplot(data=df_a, x="metodo_id", y="conidias_log10",
              color="black", alpha=0.3, size=4, ax=ax, jitter=True, order=order)
# Línea del control
ax.axhline(ctrl["conidias_log10"].mean(), color="red", ls="--", alpha=0.6, label="Control")
ax.set_xlabel("Método de extracción")
ax.set_ylabel("log₁₀(conidias/mL)")
ax.set_title("Conidias a 5.0 mg/mL por método")
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([LABEL_MET[l] for l in order])
ax.legend(loc="lower right")
save_figure_pub(fig, "obj3_conidias_log10_5mg_ml.png")
print("  ✅ obj3_conidias_log10_5mg_ml.png")

# 5b. Dosis-respuesta Maceración (log10 conidias)
fig, ax = plt.subplots(figsize=(8, 5))
order_b = sorted(df_b["conc_cat"].unique().astype(float))
sns.boxplot(data=df_b, x="conc_cat", y="conidias_log10", hue="conc_cat",
            palette="Blues", ax=ax, legend=False)
sns.stripplot(data=df_b, x="conc_cat", y="conidias_log10",
              color="black", alpha=0.3, size=4, ax=ax, jitter=True)
ax.axhline(ctrl[ctrl["metodo_extraccion"]=="maceracion"]["conidias_log10"].mean(),
           color="red", ls="--", alpha=0.6, label="Control")
ax.set_xlabel("Concentración (mg/mL)")
ax.set_ylabel("log₁₀(conidias/mL)")
ax.set_title("Maceración — dosis-respuesta en conidias")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
fig.subplots_adjust(right=0.82)
save_figure_pub(fig, "obj3_dosis_respuesta_mac_conidias.png")
print("  ✅ obj3_dosis_respuesta_mac_conidias.png")

# 5c. %INH conidias (escala log₁₀) a 5.0 mg/mL
fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=df_c, x="metodo_id", y="porcentaje_inhibicion_log10", hue="metodo_id",
            palette=COLOR_MET, ax=ax, order=order, legend=False)
sns.stripplot(data=df_c, x="metodo_id", y="porcentaje_inhibicion_log10",
              color="black", alpha=0.3, size=4, ax=ax, jitter=True, order=order)
ax.axhline(0, color="gray", ls=":", alpha=0.5)
ax.set_xlabel("Método de extracción")
ax.set_ylabel("Reducción de conidias (%) — escala log₁₀")
ax.set_title("Inhibición de conidias (log₁₀, escala de la hoja)")
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([LABEL_MET[l] for l in order])
save_figure_pub(fig, "obj3_inh_conidias_5mg_ml.png")
print("  ✅ obj3_inh_conidias_5mg_ml.png")

# 5d. Correlación %INH crecimiento vs %INH conidias
fig, ax = plt.subplots(figsize=(8, 6))
crec = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")
crec_inh = crec[~crec["es_control"] & crec["porcentaje_inhibicion"].notna()][
    ["aislado_id", "metodo_extraccion", "concentracion_mg_ml", "replica_biologica",
     "porcentaje_inhibicion"]
].rename(columns={"porcentaje_inhibicion": "inh_crecimiento"})
coni_inh = trat[trat["porcentaje_inhibicion"].notna()][
    ["aislado_id", "metodo_extraccion", "concentracion_mg_ml", "replica_biologica",
     "porcentaje_inhibicion"]
].rename(columns={"porcentaje_inhibicion": "inh_conidias"})
merged = crec_inh.merge(coni_inh, on=["aislado_id", "metodo_extraccion",
                                       "concentracion_mg_ml", "replica_biologica"])
if len(merged) > 10:
    ax.scatter(merged["inh_crecimiento"], merged["inh_conidias"],
               alpha=0.4, c="#2e86ab", edgecolors="none")
    r, p = stats.pearsonr(merged["inh_crecimiento"], merged["inh_conidias"])
    ax.set_xlabel("Inhibición crecimiento (%)")
    ax.set_ylabel("Inhibición conidias (%) — crudo")
    ax.set_title(f"Correlación %INH crecimiento vs. conidias (r={r:.3f}, p={p:.4f})")
    ax.axhline(0, color="gray", ls=":", alpha=0.5)
    ax.axvline(0, color="gray", ls=":", alpha=0.5)
    save_figure_pub(fig, "obj3_correlacion_inh_crecimiento_conidias.png")
    print("  ✅ obj3_correlacion_inh_crecimiento_conidias.png")
else:
    plt.close(fig)
    print("  ⚠ Pocos datos para correlación (%d)" % len(merged))

# 5e. Residuos Modelo A
if modelo_a_ok:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(wspace=0.4)
    resid_a = modelo_a.resid
    stats.probplot(resid_a, dist="norm", plot=axes[0])
    axes[0].set_title("Q-Q residuos (Modelo A — log10 conidias)")
    axes[1].scatter(modelo_a.fittedvalues, resid_a, alpha=0.5, c="#2e86ab", edgecolors="none")
    axes[1].axhline(0, color="gray", ls="--", alpha=0.5)
    axes[1].set_xlabel("Ajustados"); axes[1].set_ylabel("Residuos")
    axes[1].set_title("Residuos vs. ajustados")
    if len(resid_a) >= 3 and len(resid_a) <= 5000:
        _, shap_p = stats.shapiro(resid_a)
        axes[1].text(0.05, 0.95, f"Shapiro-Wilk p={shap_p:.4f}",
                     transform=axes[1].transAxes, va="top", fontsize=9,
                     bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    save_figure_pub(fig, "obj3_diagnostico_modelo_a.png")
    print("  ✅ obj3_diagnostico_modelo_a.png")


# ═══════════════════════════════════════════════════════════════════
# 6. REPORTE
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  GENERANDO REPORTE")
print("─" * 65)

with open(DIR_REPORTES / "05_objetivo3_conidias.md", "w", encoding="utf-8") as f:
    f.write("# Objetivo 3: Producción de conidias\n\n")
    f.write(f"**Fecha:** 2026-07-29\n\n")
    f.write("## Datos\n\n")
    f.write(f"- {len(trat)} observaciones con recuento de conidias\n")
    f.write(f"- Variable respuesta: log₁₀(conidias/mL) — reportada por el laboratorio\n")
    f.write(f"- Control (sin tratamiento): media log₁₀ = {ctrl['conidias_log10'].mean():.2f}\n")
    f.write(f"- {trat['aislado_id'].nunique()} aislados de Fusarium\n\n")

    f.write("## Modelo A: log₁₀(conidias) a 5.0 mg/mL\n\n")
    if modelo_a_ok:
        f.write(f"**Modelo:** Linear Mixed Model (REML)\n")
        f.write(f"- Efecto fijo: método de extracción\n")
        f.write(f"- Efecto aleatorio: aislado (intercepto)\n")
        f.write(f"- Log-Lik: {modelo_a.llf:.1f}\n\n")
        f.write("### Coeficientes\n\n")
        f.write("| Parámetro | Coeficiente | EE | z | p | IC95% |\n")
        f.write("|-----------|-------------|-----|----|----|-------|\n")
        for idx, row in coefs_a.iterrows():
            f.write(f"| {idx} | {row['Coef']:.2f} | {row['EE']:.2f} | "
                    f"{row['z']:.2f} | {row['p_valor']:.4f} | "
                    f"[{row['IC95_inf']:.1f}, {row['IC95_sup']:.1f}] |\n")
        f.write("\n")
        f.write("### Medias marginales y reducción vs control\n\n")
        f.write("| Método | log₁₀(conidias) | Conidias/mL | Reducción vs control |\n")
        f.write("|--------|:-:|:-:|:-:|\n")
        for lbl, est in [("Maceración", pred_mac_a), ("Soxhlet", pred_sox_a),
                         ("Ultrasonido", pred_ult_a)]:
            red_pct = (1 - 10**(-(ctrl_mean - est))) * 100
            f.write(f"| {lbl} | {est:.2f} | {10**est:.0f} | {red_pct:.1f}% |\n")
        f.write(f"| Control | {ctrl_mean:.2f} | {10**ctrl_mean:.0f} | — |\n")
        f.write("\n")
        f.write(f"### Varianza\n\n")
        f.write(f"- Var(aislado) = {var_iso:.3f}\n")
        f.write(f"- Var(residual) = {var_res:.3f}\n")
        f.write(f"- ICC = {icc_a:.3f}\n\n")
        f.write("### Diagnóstico de residuos\n\n")
        f.write(f"- Durbin-Watson: {dw_a:.3f} — {interpretar_dw(dw_a)}\n")
        f.write(f"- Breusch-Pagan: LM = {bp_a['lm_stat']:.2f}, ")
        f.write(f"p = {bp_a['p_val']:.4f} ")
        f.write(f"({'⚠ Heterocedástico' if bp_a['p_val'] < 0.05 else '✅ Homocedástico'})\n\n")
        f.write("### Interpretación\n\n")
        f.write("- **Maceración reduce drásticamente** la producción de conidias ")
        f.write(f"(log₁₀ {ctrl_mean:.2f} → {pred_mac_a:.2f}, ")
        f.write(f"reducción del {(1-10**(-(ctrl_mean-pred_mac_a)))*100:.0f}%)\n")
        f.write("- Soxhlet y Ultrasonido tienen poca o ninguna reducción real ")
        f.write("en la esporulación a 5.0 mg/mL\n")
        f.write("- La Maceración afecta tanto el crecimiento micelial como la ")
        f.write("esporulación, mientras que Soxhlet/Ultrasonido solo afectan el crecimiento\n\n")
    else:
        f.write("*El modelo no convergió.*\n\n")

    f.write("## Modelo B: Dosis-respuesta en Maceración\n\n")
    if modelo_b_ok:
        ctrl_mac = ctrl[ctrl["metodo_extraccion"]=="maceracion"]["conidias_log10"].mean()
        f.write(f"**Modelo:** Linear Mixed Model (REML)\n")
        f.write(f"**Control Maceración:** log₁₀ = {ctrl_mac:.2f}\n\n")
        f.write("| Concentración | log₁₀(conidias) | Conidias/mL | Reducción vs control |\n")
        f.write("|:-:|:-:|:-:|:-:|\n")
        for conc in sorted(df_b['concentracion_mg_ml'].unique()):
            if conc == conc_ref:
                est = const_b
            else:
                col = f"conc_{conc}"
                est = const_b + modelo_b.fe_params[col]
            red_pct = (1 - 10**(-(ctrl_mac - est))) * 100 if est < ctrl_mac else 0
            if red_pct < 0:
                red_pct = 0
            f.write(f"| {conc} mg/mL | {est:.2f} | {10**est:.0f} | {red_pct:.1f}% |\n")
        f.write(f"| Control | {ctrl_mac:.2f} | {10**ctrl_mac:.0f} | — |\n")
        f.write("\n### Interpretación\n\n")
        f.write("- La Maceración muestra un claro efecto dosis-dependiente ")
        f.write("sobre la producción de conidias\n")
        f.write("- A 0.2 mg/mL el efecto es mínimo o nulo\n")
        f.write("- A 5.0 mg/mL la reducción es sustancial\n\n")
    else:
        f.write("*El modelo no convergió.*\n\n")

    f.write("## Modelo C: %INH conidias (escala log₁₀ de la hoja) a 5.0 mg/mL\n\n")
    f.write("**Nota:** No se modela %INH en escala cruda porque su distribución ")
    f.write("(skew=−12.2, kurtosis=176) no es transformable. ")
    f.write("Se usa %INH_log10 (escala original del laboratorio) como variable de modelado ")
    f.write("y %INH crudo solo para visualización/interpretación.\n\n")
    if modelo_c_ok:
        f.write("**Modelo:** Linear Mixed Model (REML)\n\n")
        f.write("| Método | %INH conidias (log₁₀) |\n")
        f.write("|--------|:-:|\n")
        for lbl, est in [("Maceración", pred_mac_c), ("Soxhlet", pred_sox_c),
                         ("Ultrasonido", pred_ult_c)]:
            f.write(f"| {lbl} | {est:.1f}% |\n")
        f.write("\n### Interpretación\n\n")
        f.write("- Maceración mantiene la mayor reducción de conidias ")
        f.write("en esta escala\n")
        f.write("- La escala log₁₀ (de la hoja original) es más estable y ")
        f.write("reproducible que la escala cruda\n\n")

    f.write("## Correlación con Objetivo 2\n\n")
    if len(merged) > 10:
        f.write(f"- Correlación %INH crecimiento vs %INH conidias: r={r:.3f} (p={p:.4f})\n")
        f.write("- La baja correlación sugiere que los mecanismos de inhibición ")
        f.write("del crecimiento micelial y la esporulación son independientes ")
        f.write("o responden diferente a los extractos\n\n")

    f.write("## Figuras\n\n")
    f.write("- `obj3_conidias_log10_5mg_ml.png` — conidias por método a 5.0 mg/mL\n")
    f.write("- `obj3_dosis_respuesta_mac_conidias.png` — dosis-respuesta Maceración\n")
    f.write("- `obj3_inh_conidias_5mg_ml.png` — %INH conidias (escala log₁₀)\n")
    f.write("- `obj3_correlacion_inh_crecimiento_conidias.png` — correlación con Objetivo 2\n")
    f.write("- `obj3_diagnostico_modelo_a.png` — diagnóstico de residuos\n")

print(f"\n  ✅ Reporte guardado: {DIR_REPORTES / '05_objetivo3_conidias.md'}")
print(f"\n{'='*65}")
print("  OBJETIVO 3 — COMPLETO")
print(f"{'='*65}")
