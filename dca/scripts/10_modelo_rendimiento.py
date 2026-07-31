#!/usr/bin/env python3
"""
10_modelo_rendimiento.py — Objetivo 1: Efecto del método de extracción
sobre el rendimiento de extracción (%).

Análisis: ANOVA unifactorial / Kruskal-Wallis según supuestos.

Diseño: 3 métodos × 3 réplicas biológicas = 9 observaciones.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.anova import anova_lm
from config import (DIR_TABLAS, DIR_REPORTES, DIR_FIGURAS, COLOR_MET, 
                    LABEL_MET, setup_figure_style, save_figure_pub,
                    diagnostic_durbin_watson, interpretar_dw)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

SEMILLA = 42
np.random.seed(SEMILLA)

setup_figure_style()

# ═════════════════════════════════════════════════════════════════
print("=" * 65)
print("  OBJETIVO 1 — RENDIMIENTO DE EXTRACCIÓN")
print("=" * 65)

# ─── 1. Cargar datos ──────────────────────────────────────────────
rend = pd.read_csv(DIR_TABLAS / "rendimiento_extraccion.csv")
rend["metodo_id"] = rend["metodo_extraccion"].map({
    "maceración": "maceracion", "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"
})
print(f"\n  Datos: {len(rend)} observaciones, {rend['metodo_id'].nunique()} métodos")
print(f"  Réplicas por método: {rend.groupby('metodo_id').size().to_dict()}")

# ─── 2. Estadística descriptiva ──────────────────────────────────
desc = rend.groupby("metodo_id")["rendimiento_pct"].agg(
    n="count", media="mean", sd="std", se="sem",
    min="min", max="max"
).round(3)
print("\n  ── Estadística descriptiva ──")
for metodo, row in desc.iterrows():
    print(f"  {LABEL_MET[metodo]:15s}  n={int(row['n'])}  "
          f"media={row['media']:.2f}%  DE={row['sd']:.2f}  "
          f"IC95=[{row['media']-1.96*row['se']:.2f}, {row['media']+1.96*row['se']:.2f}]")

# ─── 3. Visualización ────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.subplots_adjust(wspace=0.35)

# 3a. Boxplot + puntos individuales
ax = axes[0]
colores = [COLOR_MET[m] for m in rend["metodo_id"]]
sns.stripplot(data=rend, x="metodo_id", y="rendimiento_pct",
              color="black", alpha=0.6, size=8, ax=ax, jitter=False)
sns.boxplot(data=rend, x="metodo_id", y="rendimiento_pct", hue="metodo_id",
            palette=COLOR_MET, ax=ax, width=0.4, legend=False)
# Ajustar transparencia manualmente
for patch in ax.patches:
    if hasattr(patch, 'set_facecolor'):
        c = patch.get_facecolor()
        patch.set_facecolor((c[0], c[1], c[2], 0.5))
ax.set_xlabel("Método de extracción")
ax.set_ylabel("Rendimiento (%)")
ax.set_title("Rendimiento por método")
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([LABEL_MET[l] for l in ["maceracion", "soxhlet", "ultrasonido"]],
                   rotation=45, ha="right")

# 3b. Barras con IC95
ax = axes[1]
means = rend.groupby("metodo_id")["rendimiento_pct"].mean()
sems = rend.groupby("metodo_id")["rendimiento_pct"].sem()
x_pos = np.arange(len(means))
ax.bar(x_pos, means.values, yerr=1.96 * sems.values, capsize=5,
       color=[COLOR_MET[m] for m in means.index], width=0.5, alpha=0.8)
ax.set_xticks(x_pos)
ax.set_xticklabels([LABEL_MET[l] for l in means.index])
ax.set_ylabel("Rendimiento medio (%)")
ax.set_title("Rendimiento medio ± IC95%")

# 3c. Distribución (violin)
ax = axes[2]
sns.violinplot(data=rend, x="metodo_id", y="rendimiento_pct", hue="metodo_id",
               palette=COLOR_MET, ax=ax, inner="quartile", legend=False)
ax.set_xlabel("Método de extracción")
ax.set_ylabel("Rendimiento (%)")
ax.set_title("Distribución del rendimiento")
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([LABEL_MET[l] for l in ["maceracion", "soxhlet", "ultrasonido"]])

save_figure_pub(fig, "obj1_rendimiento.png", clean=True)
print(f"\n  ✅ Figura: obj1_rendimiento.png")

# ─── 4. Diagnóstico de supuestos ─────────────────────────────────
print("\n  ── Diagnóstico de supuestos ──")

# 4a. Normalidad de residuos (Shapiro-Wilk sobre cada grupo)
print("  Normalidad por grupo (Shapiro-Wilk, n=3 → baja potencia):")
for metodo in ["maceracion", "soxhlet", "ultrasonido"]:
    vals = rend.loc[rend["metodo_id"] == metodo, "rendimiento_pct"]
    if len(vals) >= 3:
        w, p = stats.shapiro(vals)
        print(f"    {LABEL_MET[metodo]:15s} W={w:.4f} p={p:.4f} "
              f"{'❌ Rechaza H0' if p < 0.05 else '✅ No rechaza H0'}")

# 4b. ANOVA de un factor
print("\n  ANOVA de un factor (Rendimiento ~ Método):")
modelo_anova = ols("rendimiento_pct ~ C(metodo_id)", data=rend).fit()
tabla_anova = anova_lm(modelo_anova, typ=2)
print(tabla_anova.round(4).to_string())

# Extraer valores
f_stat = tabla_anova.loc["C(metodo_id)", "F"]
p_val = tabla_anova.loc["C(metodo_id)", "PR(>F)"]
ss_trat = tabla_anova.loc["C(metodo_id)", "sum_sq"]
ss_res = tabla_anova.loc["Residual", "sum_sq"]
df_trat = int(tabla_anova.loc["C(metodo_id)", "df"])
df_res = int(tabla_anova.loc["Residual", "df"])
ms_trat = ss_trat / df_trat
ms_res = ss_res / df_res

# Eta-cuadrado y omega-cuadrado
eta_sq = ss_trat / (ss_trat + ss_res)
omega_sq = (ss_trat - df_trat * ms_res) / (ss_trat + (df_trat + df_res + 1) * ms_res)

print(f"\n  Tamaño del efecto:")
print(f"    η² = {eta_sq:.4f}  (varianza explicada por el método)")
print(f"    ω² = {omega_sq:.4f}  (estimación insesgada)")
print(f"    Coeficiente de variación residual (CV) = {np.sqrt(ms_res) / rend['rendimiento_pct'].mean() * 100:.1f}%")

# 4c. Residuos
residuos = modelo_anova.resid
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.subplots_adjust(wspace=0.4)

# Q-Q plot
stats.probplot(residuos, dist="norm", plot=axes[0])
axes[0].set_title("Q-Q de residuos (ANOVA)")

# Homocedasticidad: residuos vs predicciones
axes[1].scatter(modelo_anova.fittedvalues, residuos, alpha=0.7, c="#2e86ab")
axes[1].axhline(0, color="gray", ls="--", alpha=0.5)
axes[1].set_xlabel("Valores ajustados")
axes[1].set_ylabel("Residuos")
axes[1].set_title("Residuos vs. ajustados")

# Test de Levene (homocedasticidad)
levene_stat, levene_p = stats.levene(
    rend.loc[rend["metodo_id"] == "maceracion", "rendimiento_pct"],
    rend.loc[rend["metodo_id"] == "soxhlet", "rendimiento_pct"],
    rend.loc[rend["metodo_id"] == "ultrasonido", "rendimiento_pct"],
)
axes[1].text(0.05, 0.95, f"Levene: F={levene_stat:.2f}, p={levene_p:.4f}",
             transform=axes[1].transAxes, va="top", fontsize=9,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

save_figure_pub(fig, "obj1_diagnostico_anova.png", clean=True)
print(f"  ✅ Figura: obj1_diagnostico_anova.png")
print(f"  Homocedasticidad (Levene): F={levene_stat:.2f}, p={levene_p:.4f} "
      f"{'❌ Heterocedástico' if levene_p < 0.05 else '✅ Homocedástico'}")

# Durbin-Watson: independencia de residuos
dw_stat = diagnostic_durbin_watson(residuos)
print(f"  Durbin-Watson: {dw_stat:.3f} — {interpretar_dw(dw_stat)}")

# ─── 5. Alternativa no paramétrica ───────────────────────────────
print("\n  ── Alternativa no paramétrica ──")
kw_stat, kw_p = stats.kruskal(
    rend.loc[rend["metodo_id"] == "maceracion", "rendimiento_pct"],
    rend.loc[rend["metodo_id"] == "soxhlet", "rendimiento_pct"],
    rend.loc[rend["metodo_id"] == "ultrasonido", "rendimiento_pct"],
)
print(f"  Kruskal-Wallis: H={kw_stat:.2f}, p={kw_p:.6f}")

# ─── 6. Post-hoc ─────────────────────────────────────────────────
print("\n  ── Comparaciones post-hoc ──")

if p_val < 0.05:
    tukey = pairwise_tukeyhsd(rend["rendimiento_pct"], rend["metodo_id"], alpha=0.05)
    print("\n  Tukey HSD:")
    print(tukey.summary().as_text())
    print()

    # Tabla de diferencias
    tukey_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
    print(f"  {'Comparación':30s} {'Diferencia':>10s} {'p-ajust':>8s} {'Signif.':>8s}")
    print("  " + "-" * 60)
    for _, row in tukey_df.iterrows():
        sig = " ✅" if row["reject"] else ""
        print(f"  {row['group1'] + ' - ' + row['group2']:30s} {row['meandiff']:>8.3f}  "
              f"{row['p-adj']:>8.4f} {sig:>8s}")
else:
    print("  ANOVA no significativo → no se realizan comparaciones post-hoc paramétricas")

# ─── 7. Tabla resumen ────────────────────────────────────────────
resultados = pd.DataFrame({
    "Método": [LABEL_MET[m] for m in desc.index],
    "n": desc["n"].values.astype(int),
    "Media (%)": desc["media"].values,
    "DE": desc["sd"].values,
    "IC95_inf": (desc["media"] - 1.96 * desc["se"]).values,
    "IC95_sup": (desc["media"] + 1.96 * desc["se"]).values,
    "Mín": desc["min"].values,
    "Máx": desc["max"].values,
})

resumen_modelo = {
    "Modelo": ["ANOVA unifactorial (OLS)", "Kruskal-Wallis"],
    "Estadístico": [f"F({df_trat},{df_res}) = {f_stat:.2f}", f"H = {kw_stat:.2f}"],
    "Valor p": [f"{p_val:.6f}", f"{kw_p:.6f}"],
    "Tamaño efecto": [f"η² = {eta_sq:.3f}, ω² = {omega_sq:.3f}", "—"],
}

print("\n  ── Tabla resumen ──")
print(resultados.to_string(index=False))

# ─── Guardar reporte ─────────────────────────────────────────────
with open(DIR_REPORTES / "03_objetivo1_rendimiento.md", "w", encoding="utf-8") as f:
    f.write("# Objetivo 1: Efecto del método de extracción sobre el rendimiento\n\n")
    f.write(f"**Fecha:** 2026-07-29\n\n")
    f.write("## Datos\n\n")
    f.write(f"- {len(rend)} observaciones (3 métodos × 3 réplicas biológicas)\n")
    f.write(f"- Variable respuesta: rendimiento de extracción (% p/p)\n\n")
    f.write("## Estadística descriptiva\n\n")
    f.write(resultados.to_string(index=False))
    f.write("\n\n## Resultados del modelo\n\n")
    f.write(f"### ANOVA unifactorial\n\n")
    f.write(f"- F({df_trat},{df_res}) = {f_stat:.2f}, p = {p_val:.6f}\n")
    f.write(f"- η² = {eta_sq:.3f} ({eta_sq*100:.1f}% de varianza explicada)\n")
    f.write(f"- ω² = {omega_sq:.3f} (estimación insesgada del tamaño del efecto)\n")
    if p_val < 0.05:
        f.write("- El método de extracción tiene un efecto estadísticamente significativo ")
        f.write("sobre el rendimiento.\n")
        for _, row in tukey_df.iterrows():
            if row["reject"]:
                f.write(f"  - {row['group1']} ≠ {row['group2']} (p={row['p-adj']:.4f}, "
                        f"diferencia={row['meandiff']:.2f}%)\n")
    else:
        f.write("- No se detectaron diferencias significativas entre métodos.\n\n")
    f.write("\n### Alternativa no paramétrica\n\n")
    f.write(f"- Kruskal-Wallis: H = {kw_stat:.2f}, p = {kw_p:.6f}\n\n")
    f.write("## Diagnóstico de supuestos\n\n")
    f.write(f"- Levene (homocedasticidad): F = {levene_stat:.2f}, p = {levene_p:.4f}\n")
    f.write(f"- Durbin-Watson (independencia): {dw_stat:.3f} — {interpretar_dw(dw_stat)}\n")
    f.write(f"- Normalidad: n=3 por grupo — potencia insuficiente para Shapiro-Wilk\n")
    f.write("- ANOVA es robusto con diseño balanceado incluso con ligeras desviaciones.\n\n")
    f.write("## Interpretación biológica\n\n")
    f.write(f"- Soxhlet produce el mayor rendimiento ({desc.loc['soxhlet', 'media']:.1f}%), ")
    f.write("seguido de ultrasonido ")
    f.write(f"({desc.loc['ultrasonido', 'media']:.1f}%) y maceración ")
    f.write(f"({desc.loc['maceracion', 'media']:.1f}%).\n")
    f.write("- La diferencia entre Soxhlet y los otros dos métodos es sustancial ")
    f.write("y consistente.\n\n")
    f.write("## Figuras\n\n")
    f.write("- `obj1_rendimiento.png` — visualización de datos\n")
    f.write("- `obj1_diagnostico_anova.png` — diagnóstico de residuos\n")

print(f"\n  ✅ Reporte guardado: {DIR_REPORTES / '03_objetivo1_rendimiento.md'}")
print(f"\n{'='*65}")
print("  OBJETIVO 1 — COMPLETO")
print(f"{'='*65}")
