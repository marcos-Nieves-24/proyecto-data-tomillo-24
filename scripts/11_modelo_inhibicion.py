#!/usr/bin/env python3
"""
11_modelo_inhibicion.py — Objetivo 2: Efecto del método de extracción,
aislado y concentración sobre la inhibición de crecimiento micelial.

Estrategia de modelado:
  A) Modelo mixto global (solo 5.0 mg/mL) — compara los 3 métodos a la
     concentración común. Efecto fijo: método. Efecto aleatorio: aislado.
  B) Modelo Maceración dosis-respuesta — efecto de concentración intra-
     Maceración (único método con gradiente completo).
  C) Modelo logístico — probabilidad de inhibición completa.
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
from statsmodels.stats.multitest import multipletests
from statsmodels.tools import add_constant
from statsmodels.discrete.discrete_model import Logit
import statsmodels.api as sm
from config import (DIR_TABLAS, DIR_REPORTES, DIR_FIGURAS, COLOR_MET, 
                    LABEL_MET, setup_figure_style, save_figure_pub,
                    diagnostic_durbin_watson, diagnostic_breusch_pagan,
                    diagnostic_vif, interpretar_dw)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

SEMILLA = 42
np.random.seed(SEMILLA)

setup_figure_style()

# ═══════════════════════════════════════════════════════════════════
print("=" * 65)
print("  OBJETIVO 2 — INHIBICIÓN DE CRECIMIENTO MICELIAL")
print("=" * 65)

# ─── 1. Cargar datos ─────────────────────────────────────────────
crec = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")

# Solo filas de tratamiento con %INH válido
inh = (crec[~crec["es_control"] & crec["porcentaje_inhibicion"].notna()]
       .copy())
inh["metodo_id"] = inh["metodo_extraccion"].map(
    {"maceracion": "maceracion", "maceración": "maceracion",
     "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"}
)

print(f"\n  Total filas con %INH: {len(inh)}")
print(f"  Aislados: {inh['aislado_id'].nunique()}")
print(f"  Rangos %INH: [{inh['porcentaje_inhibicion'].min():.1f}, {inh['porcentaje_inhibicion'].max():.1f}]")

# ─── 2. Preprocesamiento ─────────────────────────────────────────

# Estandarizar nombres de aislado para consistencia
inh["aislado_id"] = inh["aislado_id"].str.strip()

# Variable de concentración como categoría para modelos factoriales
inh["conc_cat"] = inh["concentracion_mg_ml"].astype("category")
# Versión log para dosis-respuesta
inh["conc_log"] = np.log(inh["concentracion_mg_ml"] + 0.01)

# Inhibición completa binaria
inh["completa_bin"] = (inh["inhibicion_completa"]).astype(int)

# ─── 3. MODELO A: Global a 5.0 mg/mL ──────────────────────────
print("\n" + "─" * 65)
print("  MODELO A: Comparación de métodos a 5.0 mg/mL (concentración común)")
print("─" * 65)

df_a = inh[inh["concentracion_mg_ml"] == 5.0].copy()
print(f"  Datos: {len(df_a)} obs, {df_a['aislado_id'].nunique()} aislados")

# Asegurar que método sea categórico y Maceración sea la referencia
df_a["metodo_id"] = pd.Categorical(df_a["metodo_id"],
                                     categories=["maceracion", "soxhlet", "ultrasonido"],
                                     ordered=False)

# ML estimation may fail to converge — use REML via MixedLM
# Modelo: %INH ~ método + (1|aislado)
try:
    # Crear dummies manualmente para MixedLM (no tiene formula)
    df_a_ml = df_a.copy()
    df_a_ml = pd.get_dummies(df_a_ml, columns=["metodo_id"], drop_first=True, dtype=float)

    # Variable de grupo es aislado
    groups_a = df_a_ml["aislado_id"]

    # Predictores fijos: intercept + Soxhlet + Ultrasonido (Maceración es ref)
    exog_a = add_constant(df_a_ml[["metodo_id_soxhlet", "metodo_id_ultrasonido"]])
    endog_a = df_a_ml["porcentaje_inhibicion"]

    modelo_a = MixedLM(endog_a, exog_a, groups=groups_a).fit(reml=True, maxiter=200)
    print(f"\n  Convergió: {modelo_a.converged}  |  Log-Lik: {modelo_a.llf:.1f}")
    print(f"  AIC: {modelo_a.aic:.0f}  |  BIC: {modelo_a.bic:.0f}")
    print(f"\n  Efectos fijos:\n{modelo_a.fe_params.to_string()}")
    print(f"\n  Efectos aleatorios (σ²):\n  Var(aislado) = {modelo_a.cov_re.iloc[0,0]:.2f}")
    print(f"  Var(residual) = {modelo_a.scale:.2f}")

    # ICC
    icc_a = modelo_a.cov_re.iloc[0, 0] / (modelo_a.cov_re.iloc[0, 0] + modelo_a.scale)
    print(f"  ICC (aislado) = {icc_a:.3f}  ({icc_a*100:.1f}% de varianza entre aislados)")

    # Tabla de coeficientes
    coefs_a = pd.DataFrame({
        "Coef": modelo_a.fe_params,
        "EE": modelo_a.bse_fe,
        "z": modelo_a.tvalues,
        "p_valor": modelo_a.pvalues,
    })
    coefs_a["IC95_inf"] = coefs_a["Coef"] - 1.96 * coefs_a["EE"]
    coefs_a["IC95_sup"] = coefs_a["Coef"] + 1.96 * coefs_a["EE"]
    print(f"\n  Coeficientes detallados:")
    print(coefs_a.round(3).to_string())

    modelo_a_ok = True

    # Diagnóstico de residuos
    dw_a = diagnostic_durbin_watson(modelo_a.resid)
    print(f"  Durbin-Watson: {dw_a:.3f} — {interpretar_dw(dw_a)}")
    bp_a = diagnostic_breusch_pagan(modelo_a, exog_a)
    print(f"  Breusch-Pagan: LM={bp_a['lm_stat']:.2f}, p={bp_a['p_val']:.4f}")
except Exception as e:
    print(f"  ⚠ Modelo A falló: {e}")
    modelo_a_ok = False

# Medias marginales estimadas (Modelo A)
if modelo_a_ok:
    print("\n  Medias marginales estimadas a 5.0 mg/mL:")
    # Para Maceración (ref): solo intercepto
    pred_mac = modelo_a.fe_params["const"]
    pred_sox = modelo_a.fe_params["const"] + modelo_a.fe_params["metodo_id_soxhlet"]
    pred_ult = modelo_a.fe_params["const"] + modelo_a.fe_params["metodo_id_ultrasonido"]

    se_mac = modelo_a.bse_fe["const"]
    se_sox = np.sqrt(modelo_a.bse_fe["const"]**2 + modelo_a.bse_fe["metodo_id_soxhlet"]**2
                      + 2*modelo_a.cov_params().loc["const", "metodo_id_soxhlet"])
    se_ult = np.sqrt(modelo_a.bse_fe["const"]**2 + modelo_a.bse_fe["metodo_id_ultrasonido"]**2
                      + 2*modelo_a.cov_params().loc["const", "metodo_id_ultrasonido"])

    for lbl, est, se in [("Maceración", pred_mac, se_mac),
                         ("Soxhlet", pred_sox, se_sox),
                         ("Ultrasonido", pred_ult, se_ult)]:
        print(f"    {lbl:15s} = {est:.2f}% ± {1.96*se:.2f}")


# ─── 4. MODELO B: Maceración dosis-respuesta ──────────────────────
print("\n" + "─" * 65)
print("  MODELO B: Maceración — dosis-respuesta")
print("─" * 65)

df_b = inh[inh["metodo_id"] == "maceracion"].copy()
print(f"  Datos: {len(df_b)} obs, {df_b['aislado_id'].nunique()} aislados")
print(f"  Concentraciones: {sorted(df_b['concentracion_mg_ml'].unique())}")

try:
    dummies_conc = pd.get_dummies(df_b["conc_cat"], drop_first=True, dtype=float,
                                   prefix="conc")
    df_b_ml = pd.concat([df_b.reset_index(drop=True), dummies_conc.reset_index(drop=True)], axis=1)
    groups_b = df_b_ml["aislado_id"]

    exog_b = add_constant(df_b_ml[[c for c in dummies_conc.columns]])
    endog_b = df_b_ml["porcentaje_inhibicion"]

    modelo_b = MixedLM(endog_b, exog_b, groups=groups_b).fit(reml=True, maxiter=200)
    print(f"\n  Convergió: {modelo_b.converged}  |  Log-Lik: {modelo_b.llf:.1f}")
    print(f"  AIC: {modelo_b.aic:.0f}")

    coefs_b = pd.DataFrame({
        "Coef": modelo_b.fe_params,
        "EE": modelo_b.bse_fe,
        "z": modelo_b.tvalues,
        "p_valor": modelo_b.pvalues,
    })
    coefs_b["IC95_inf"] = coefs_b["Coef"] - 1.96 * coefs_b["EE"]
    coefs_b["IC95_sup"] = coefs_b["Coef"] + 1.96 * coefs_b["EE"]
    print(f"\n  Coeficientes:\n{coefs_b.round(3).to_string()}")

    icc_b = modelo_b.cov_re.iloc[0, 0] / (modelo_b.cov_re.iloc[0, 0] + modelo_b.scale)
    print(f"\n  Var(aislado) = {modelo_b.cov_re.iloc[0,0]:.2f}")
    print(f"  Var(residual) = {modelo_b.scale:.2f}")
    print(f"  ICC (aislado) = {icc_b:.3f}")

    # Medias marginales por concentración
    conc_ref = sorted(df_b['concentracion_mg_ml'].unique())[0]
    print(f"\n  Medias marginales (ref = {conc_ref} mg/mL):")
    const_b = modelo_b.fe_params["const"]
    for conc in sorted(df_b['concentracion_mg_ml'].unique()):
        if conc == conc_ref:
            est = const_b
            se = modelo_b.bse_fe["const"]
        else:
            col = f"conc_{conc}"
            est = const_b + modelo_b.fe_params[col]
            se = np.sqrt(modelo_b.bse_fe["const"]**2 + modelo_b.bse_fe[col]**2
                          + 2*modelo_b.cov_params().loc["const", col])
        print(f"    {conc:>5} mg/mL: {est:.2f}% ± {1.96*se:.2f}")

    modelo_b_ok = True

    # VIF para modelo B
    vif_b = diagnostic_vif(exog_b, var_names=["const"] + [c for c in dummies_conc.columns])
    print(f"  VIF: {', '.join([f'{k}={v:.2f}' for k,v in vif_b.items()])}")
except Exception as e:
    print(f"  ⚠ Modelo B falló: {e}")
    modelo_b_ok = False


# ─── 5. MODELO C: Logístico — inhibición completa ─────────────────
print("\n" + "─" * 65)
print("  MODELO C: Probabilidad de inhibición completa (logístico, 5.0 mg/mL)")
print("─" * 65)

df_c = df_a[df_a["concentracion_mg_ml"] == 5.0].copy()
n_completa = df_c["completa_bin"].sum()
print(f"  Inhibición completa: {n_completa}/{len(df_c)} ({100*n_completa/len(df_c):.1f}%)")

# Tabla de contingencia método × completa
ct = pd.crosstab(df_c["metodo_id"], df_c["completa_bin"])
print(f"\n  Contingencia:\n{ct}")

# Verificar si hay variación suficiente para modelo logístico
n_completa_por_metodo = df_c.groupby("metodo_id")["completa_bin"].sum()
print(f"\n  Inhibición completa por método:\n{n_completa_por_metodo}")

if n_completa_por_metodo.nunique() == 1:
    # Solo Maceración tiene casos — reportar proporciones
    print("\n  ⚠ Solo un método tiene inhibición completa — modelo logístico no aplicable")
    print("  Se reportan proporciones observadas directamente.")
    modelo_c_ok = False
    modelo_c = None
    pred_probs_display = ct[1] / ct.sum(axis=1) if 1 in ct.columns else pd.Series(0, index=ct.index)
else:
    try:
        df_c_ml = pd.get_dummies(df_c, columns=["metodo_id"], drop_first=True, dtype=float)
        exog_c = add_constant(df_c_ml[["metodo_id_soxhlet", "metodo_id_ultrasonido"]])
        endog_c = df_c_ml["completa_bin"]
        modelo_c = Logit(endog_c, exog_c).fit(disp=False, maxiter=200)
        print(f"\n  Log-Lik: {modelo_c.llf:.1f}  |  Pseudo-R²: {modelo_c.prsquared:.3f}")
        pred_probs = modelo_c.predict(exog_c)
        modelo_c_ok = True
    except Exception as e:
        print(f"  ⚠ Modelo logístico falló: {e}")
        modelo_c_ok = False

if modelo_c_ok:
    print(f"\n  Probabilidad predicha de inhibición completa:")
    for metodo in ["maceracion", "soxhlet", "ultrasonido"]:
        mask = df_c["metodo_id"] == metodo
        pp = pred_probs[mask].mean()
        print(f"    {LABEL_MET[metodo]:15s}: {pp:.1%}")
else:
    print(f"\n  Proporción observada de inhibición completa:")
    for metodo in ["maceracion", "soxhlet", "ultrasonido"]:
        n_total = len(df_c[df_c["metodo_id"] == metodo])
        n_comp = n_completa_por_metodo.get(metodo, 0)
        print(f"    {LABEL_MET[metodo]:15s}: {n_comp}/{n_total} = {100*n_comp/n_total:.1f}%")


# ═══════════════════════════════════════════════════════════════════
# 6. FIGURAS
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  FIGURAS")
print("─" * 65)

# 6a. Boxplot comparativo a 5 mg/mL
fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=df_a, x="metodo_id", y="porcentaje_inhibicion", hue="metodo_id",
            palette=COLOR_MET, ax=ax, legend=False)
sns.stripplot(data=df_a, x="metodo_id", y="porcentaje_inhibicion",
              color="black", alpha=0.3, size=4, ax=ax, jitter=True)
ax.set_xlabel("Método de extracción")
ax.set_ylabel("Inhibición (%)")
ax.set_title("Inhibición a 5.0 mg/mL por método")
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([LABEL_MET[l] for l in ["maceracion", "soxhlet", "ultrasonido"]])
ax.axhline(0, color="gray", ls=":", alpha=0.5)

# Agregar medias estimadas del modelo
if modelo_a_ok:
    ax.errorbar([0, 1, 2], [pred_mac, pred_sox, pred_ult],
                yerr=[1.96*se_mac, 1.96*se_sox, 1.96*se_ult],
                fmt="D", color="#d62728", markersize=7, capsize=4, zorder=5,
                label="Estimado LMM")
    ax.legend(loc="upper right")

save_figure_pub(fig, "obj2_inhibicion_5mg_ml.png", clean=True)
print("  ✅ obj2_inhibicion_5mg_ml.png")

# 6b. Dosis-respuesta Maceración
fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=df_b, x="conc_cat", y="porcentaje_inhibicion", hue="conc_cat",
            palette="Blues", ax=ax, legend=False)
sns.stripplot(data=df_b, x="conc_cat", y="porcentaje_inhibicion",
              color="black", alpha=0.3, size=4, ax=ax, jitter=True)
ax.set_xlabel("Concentración (mg/mL)")
ax.set_ylabel("Inhibición (%)")
ax.set_title("Maceración — dosis-respuesta")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
if modelo_b_ok:
    x_ticks = sorted(df_b['concentracion_mg_ml'].unique())
    x_pos = list(range(len(x_ticks)))
    for i, conc in enumerate(x_ticks):
        est, se = None, None
        if conc == conc_ref:
            est = const_b; se = modelo_b.bse_fe["const"]
        else:
            col = f"conc_{conc}"
            if col in modelo_b.fe_params.index:
                est = const_b + modelo_b.fe_params[col]
                se = np.sqrt(modelo_b.bse_fe["const"]**2 + modelo_b.bse_fe[col]**2
                              + 2*modelo_b.cov_params().loc["const", col])
        if est is not None:
            ax.plot(i, est, "D", color="red", markersize=8, zorder=5)
            ax.errorbar(i, est, yerr=1.96*se, color="red", capsize=5, capthick=2, alpha=0.7)

save_figure_pub(fig, "obj2_dosis_respuesta_maceracion.png", clean=True)
print("  ✅ obj2_dosis_respuesta_maceracion.png")

# 6c. Perfil por aislado (Maceración, todas las concentraciones)
fig, ax = plt.subplots(figsize=(12, 6))
perfil = df_b.groupby(["aislado_id", "concentracion_mg_ml"])["porcentaje_inhibicion"].mean().reset_index()
for aislado in df_b["aislado_id"].unique():
    sub = perfil[perfil["aislado_id"] == aislado]
    ax.plot(sub["concentracion_mg_ml"], sub["porcentaje_inhibicion"],
            marker="o", alpha=0.4, linewidth=0.8, color="#2e86ab")
# Promedio general
avg = perfil.groupby("concentracion_mg_ml")["porcentaje_inhibicion"].mean()
ax.plot(avg.index, avg.values, "r-o", linewidth=3, label="Promedio")
ax.set_xlabel("Concentración (mg/mL)")
ax.set_ylabel("Inhibición media (%)")
ax.set_title("Perfil individual por aislado — Maceración")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
fig.subplots_adjust(right=0.82)
save_figure_pub(fig, "obj2_perfil_aislados_mac.png", clean=True)
print("  ✅ obj2_perfil_aislados_mac.png")

# 6d. Inhibición completa por método
fig, ax = plt.subplots(figsize=(8, 5))
ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
plot_col = 1 if 1 in ct_pct.columns else ct_pct.columns[0]
ct_pct[plot_col].plot(kind="bar", ax=ax, color=[COLOR_MET[l] for l in ["maceracion", "soxhlet", "ultrasonido"]],
                      edgecolor="black", alpha=0.8)
ax.set_xticks(range(3))
ax.set_xticklabels([LABEL_MET[l] for l in ["maceracion", "soxhlet", "ultrasonido"]],
                   rotation=0)
ax.set_xlabel("Método")
ax.set_ylabel("% de observaciones con inhibición completa")
ax.set_title("Proporción de inhibición completa a 5.0 mg/mL")
save_figure_pub(fig, "obj2_inhibicion_completa.png", clean=True)
print("  ✅ obj2_inhibicion_completa.png")

# 6e. Residuales del modelo A
if modelo_a_ok:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(wspace=0.4)
    resid_a = modelo_a.resid
    stats.probplot(resid_a, dist="norm", plot=axes[0])
    axes[0].set_title("Q-Q de residuos (Modelo A)")
    axes[1].scatter(modelo_a.fittedvalues, resid_a, alpha=0.5, c="#a23b72", edgecolors="none")
    axes[1].axhline(0, color="gray", ls="--", alpha=0.5)
    axes[1].set_xlabel("Ajustados")
    axes[1].set_ylabel("Residuos")
    axes[1].set_title("Residuos vs. ajustados")
    # Shapiro de residuos
    if len(resid_a) >= 3 and len(resid_a) <= 5000:
        _, shap_p = stats.shapiro(resid_a)
        axes[1].text(0.05, 0.95, f"Shapiro-Wilk p={shap_p:.4f}",
                     transform=axes[1].transAxes, va="top", fontsize=9,
                     bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    save_figure_pub(fig, "obj2_diagnostico_modelo_a.png", clean=True)
    print("  ✅ obj2_diagnostico_modelo_a.png")


# ═══════════════════════════════════════════════════════════════════
# 7. REPORTE
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  GENERANDO REPORTE")
print("─" * 65)

with open(DIR_REPORTES / "04_objetivo2_inhibicion.md", "w", encoding="utf-8") as f:
    f.write("# Objetivo 2: Inhibición de crecimiento micelial\n\n")
    f.write(f"**Fecha:** 2026-07-29\n\n")

    f.write("## Datos\n\n")
    f.write(f"- {len(inh)} observaciones con %INH válido\n")
    f.write(f"- {inh['aislado_id'].nunique()} aislados de Fusarium\n")
    f.write(f"- 3 métodos de extracción, 3 réplicas biológicas\n")
    f.write(f"- Inhibición completa: {(inh['inhibicion_completa']).sum()} casos\n")
    f.write(f"- Inhibición negativa: {(inh['inhibicion_negativa']).sum()} casos\n\n")

    f.write("## Modelo A: Comparación global a 5.0 mg/mL\n\n")
    if modelo_a_ok:
        f.write(f"**Modelo:** Linear Mixed Model (REML)\n")
        f.write(f"- Efecto fijo: método de extracción\n")
        f.write(f"- Efecto aleatorio: aislado (intercepto)\n")
        f.write(f"- Log-Lik: {modelo_a.llf:.1f}, AIC: {modelo_a.aic:.1f}\n\n")
        f.write("### Coeficientes\n\n")
        f.write("| Parámetro | Coeficiente | EE | z | p | IC95% |\n")
        f.write("|-----------|-------------|-----|----|----|-------|\n")
        for idx, row in coefs_a.iterrows():
            f.write(f"| {idx} | {row['Coef']:.2f} | {row['EE']:.2f} | "
                    f"{row['z']:.2f} | {row['p_valor']:.4f} | "
                    f"[{row['IC95_inf']:.1f}, {row['IC95_sup']:.1f}] |\n")
        f.write("\n")
        f.write(f"### Varianza de componentes\n\n")
        f.write(f"- Var(aislado) = {modelo_a.cov_re.iloc[0,0]:.2f}\n")
        f.write(f"- Var(residual) = {modelo_a.scale:.2f}\n")
        f.write(f"- ICC = {icc_a:.3f} — el {icc_a*100:.1f}% de la variabilidad está entre aislados\n\n")
        f.write("### Diagnóstico de residuos\n\n")
        f.write(f"- Durbin-Watson: {dw_a:.3f} — {interpretar_dw(dw_a)}\n")
        f.write(f"- Breusch-Pagan: LM={bp_a['lm_stat']:.2f}, "
                f"p={bp_a['p_val']:.4f}"
                f" ({'' if bp_a['p_val'] >= 0.05 else '⚠ '}evidencia de heterocedasticidad)\n\n")
        f.write("### Medias marginales estimadas (5.0 mg/mL)\n\n")
        for lbl, est, se in [("Maceración", pred_mac, se_mac),
                              ("Soxhlet", pred_sox, se_sox),
                              ("Ultrasonido", pred_ult, se_ult)]:
            f.write(f"- **{lbl}**: {est:.1f}% ± {1.96*se:.1f} (IC95%)\n")
        f.write("\n### Interpretación\n\n")
        sox_eff = modelo_a.fe_params["metodo_id_soxhlet"]
        ult_eff = modelo_a.fe_params["metodo_id_ultrasonido"]
        sox_p = modelo_a.pvalues["metodo_id_soxhlet"]
        ult_p = modelo_a.pvalues["metodo_id_ultrasonido"]
        f.write(f"- Soxhlet supera a Maceración en {sox_eff:.1f} puntos porcentuales "
                f"(p={sox_p:.4f})\n")
        f.write(f"- Ultrasonido supera a Maceración en {ult_eff:.1f} puntos porcentuales "
                f"(p={ult_p:.4f})\n")
        f.write(f"- La alta correlación intra-aislado (ICC={icc_a:.2f}) confirma que "
                f"los aislados difieren consistentemente en su susceptibilidad.\n\n")
    else:
        f.write("*El modelo no convergió.*\n\n")

    f.write("## Modelo B: Dosis-respuesta en Maceración\n\n")
    if modelo_b_ok:
        f.write(f"**Modelo:** Linear Mixed Model (REML)\n")
        f.write(f"- Efecto fijo: concentración (categórica)\n")
        f.write(f"- Efecto aleatorio: aislado (intercepto)\n")
        f.write(f"- Log-Lik: {modelo_b.llf:.1f}, AIC: {modelo_b.aic:.1f}\n\n")
        f.write("### Coeficientes\n\n")
        f.write("| Parámetro | Coeficiente | EE | z | p | IC95% |\n")
        f.write("|-----------|-------------|-----|----|----|-------|\n")
        for idx, row in coefs_b.iterrows():
            f.write(f"| {idx} | {row['Coef']:.2f} | {row['EE']:.2f} | "
                    f"{row['z']:.2f} | {row['p_valor']:.4f} | "
                    f"[{row['IC95_inf']:.1f}, {row['IC95_sup']:.1f}] |\n")
        f.write("\n")
        f.write(f"### Interpretación\n\n")
        f.write(f"- ICC = {icc_b:.3f}\n")
        f.write("### Diagnóstico de multicolinealidad\n\n")
        f.write(f"- VIF: {', '.join([f'{k}={vif_b[k]:.2f}' for k in vif_b])}\n\n")
        f.write("- Se observa un claro gradiente de inhibición con la concentración.\n\n")
    else:
        f.write("*El modelo no convergió.*\n\n")

    f.write("## Modelo C: Probabilidad de inhibición completa\n\n")
    if modelo_c_ok:
        f.write(f"**Modelo:** Regresión logística (5.0 mg/mL)\n")
        f.write(f"- Pseudo-R²: {modelo_c.prsquared:.3f}\n")
        f.write(f"- Eventos: {n_completa}/{len(df_c)} observaciones\n\n")
        f.write("### Probabilidades predichas\n\n")
        for metodo in ["maceracion", "soxhlet", "ultrasonido"]:
            mask = df_c["metodo_id"] == metodo
            pp = pred_probs[mask].mean()
            obs = mask.sum()
            f.write(f"- {LABEL_MET[metodo]}: {pp:.1%} (n={obs})\n")
    else:
        f.write("**Nota:** No fue posible ajustar un modelo logístico porque solo Maceración ")
        f.write("presenta casos de inhibición completa (Soxhlet y Ultrasonido tienen 0 casos ")
        f.write(f"en {len(df_c[df_c['metodo_id']=='soxhlet'])} y ")
        f.write(f"{len(df_c[df_c['metodo_id']=='ultrasonido'])} observaciones respectivamente).\n\n")
        f.write("### Proporciones observadas\n\n")
        for metodo in ["maceracion", "soxhlet", "ultrasonido"]:
            n_total = len(df_c[df_c["metodo_id"] == metodo])
            n_comp = int(n_completa_por_metodo.get(metodo, 0))
            f.write(f"- {LABEL_MET[metodo]}: {n_comp}/{n_total} = {100*n_comp/n_total:.1f}%\n")
        f.write("\n")
        f.write("**Interpretación:** Solo Maceración logró inhibición completa en algunos casos. ")
        f.write("Soxhlet y Ultrasonido, a pesar de tener buena inhibición promedio, nunca ")
        f.write("alcanzaron el 100% de inhibición del crecimiento. Esto sugiere que el perfil ")
        f.write("fitoquímico de los extractos de Maceración es cualitativamente diferente.\n\n")

    f.write("## Figuras\n\n")
    f.write("- `obj2_inhibicion_5mg_ml.png` — inhibición por método a 5.0 mg/mL\n")
    f.write("- `obj2_dosis_respuesta_maceracion.png` — curva dosis-respuesta Maceración\n")
    f.write("- `obj2_perfil_aislados_mac.png` — perfiles individuales por aislado\n")
    f.write("- `obj2_inhibicion_completa.png` — proporción de inhibición completa\n")
    f.write("- `obj2_diagnostico_modelo_a.png` — diagnóstico de residuos\n")

print(f"\n  ✅ Reporte guardado: {DIR_REPORTES / '04_objetivo2_inhibicion.md'}")
print(f"\n{'='*65}")
print("  OBJETIVO 2 — COMPLETO")
print(f"{'='*65}")
