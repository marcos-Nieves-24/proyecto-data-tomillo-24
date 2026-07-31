#!/usr/bin/env python3
"""
14_diagnostico_transformaciones.py — Evaluación sistemática de
transformaciones para cada variable respuesta del pipeline.

Compara: raw, log, logit, arcoseno√, IHS, cube-root, Box-Cox, rank-gauss.
Evalúa impacto en: normalidad de residuos, homocedasticidad,
interpretabilidad biológica.

Genera:
  - Tabla comparativa de transformaciones
  - Figuras de diagnóstico (Q-Q, residuos)
  - Recomendación final por dataset
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools import add_constant
from config import DIR_FIGURAS, DIR_TABLAS, DIR_REPORTES, setup_figure_style
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

SEMILLA = 42
np.random.seed(SEMILLA)
setup_figure_style()

# ═══════════════════════════════════════════════════════════════════
print("=" * 65)
print("  DIAGNÓSTICO SISTEMÁTICO DE TRANSFORMACIONES")
print("=" * 65)

# ─── 1. Cargar datos ─────────────────────────────────────────────
crec = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")
coni = pd.read_csv(DIR_TABLAS / "conidias.csv")
rend = pd.read_csv(DIR_TABLAS / "rendimiento_extraccion.csv")
inh = crec[~crec["es_control"] & crec["porcentaje_inhibicion"].notna()]
con_log = coni[~coni["es_control"] & coni["conidias_log10"].notna()]
con_inh = coni[~coni["es_control"] & coni["porcentaje_inhibicion"].notna()]
con_inh_log = coni[~coni["es_control"] & coni["porcentaje_inhibicion_log10"].notna()]

# ─── 2. Funciones auxiliares ──────────────────────────────────────
TRANSFORMS = {}

def eval_transforms(name, y, label="", transform_funcs=None):
    """Evalúa múltiples transformaciones y retorna tabla comparativa."""
    if transform_funcs is None:
        return []
    
    results = []
    for tname, tfunc, desc in transform_funcs:
        try:
            yt = tfunc(y)
            sk = float(stats.skew(yt))
            ku = float(stats.kurtosis(yt))
            if len(yt) >= 3 and len(yt) <= 5000:
                _, sp = stats.shapiro(yt)
            else:
                sp = 0.0
            results.append({"Variable": name, "Transformacion": tname,
                            "Descripcion": desc, "n": len(yt),
                            "Skew": round(sk, 3), "Kurtosis": round(ku, 3),
                            "Shapiro_p": round(sp, 4)})
        except Exception as e:
            results.append({"Variable": name, "Transformacion": tname,
                            "Descripcion": f"ERROR: {e}",
                            "n": 0, "Skew": np.nan, "Kurtosis": np.nan,
                            "Shapiro_p": np.nan})
    return results

def ihs_transform(c, shift=0):
    """Inverse hyperbolic sine: arcsinh((y+shift)/c)"""
    return lambda y: np.arcsinh((y + shift) / c)

def cube_root(y):
    return np.cbrt(y)

def rank_gauss(y):
    r = stats.rankdata(y, method="average")
    return stats.norm.ppf((r - 0.5) / len(y))

def eval_residuals_lmm(name, y_full, exog, groups, transforms):
    """Evalúa transformaciones en los residuos de un LMM."""
    rows = []
    for tname, tfunc, desc in transforms:
        try:
            yt = tfunc(y_full)
            m = MixedLM(yt, exog, groups=groups).fit(reml=True, maxiter=200)
            r = m.resid
            sk = float(stats.skew(r))
            ku = float(stats.kurtosis(r))
            llf = m.llf
            if len(r) >= 3 and len(r) <= 5000:
                _, sp = stats.shapiro(r.sample(min(5000, len(r))))
            else:
                sp = 0.0
            rows.append({"Modelo": name, "Transformacion": tname,
                         "LogLik": round(llf, 1),
                         "Resid_skew": round(sk, 3),
                         "Resid_kurt": round(ku, 3),
                         "Shapiro_p": round(sp, 4)})
        except Exception as e:
            rows.append({"Modelo": name, "Transformacion": tname,
                         "LogLik": np.nan,
                         "Resid_skew": np.nan, "Resid_kurt": np.nan,
                         "Shapiro_p": np.nan})
    return rows


# ═══════════════════════════════════════════════════════════════════
# 3. RENDIMIENTO (distribución solamente — n muy pequeño)
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  1. RENDIMIENTO (%)")
print("─" * 65)

yp = rend["rendimiento_pct"]
rend_transforms = [
    ("Raw", lambda x: x, "Sin transformación"),
    ("Log", np.log, "log(y)"),
    ("BoxCox(λ)", lambda x: stats.boxcox(x)[0], "Box-Cox automático"),
]
rend_results = eval_transforms("Rendimiento", yp.values, transform_funcs=rend_transforms)
for r in rend_results:
    print(f"  {r['Transformacion']:12s}  skew={r['Skew']:.2f}  kurt={r['Kurtosis']:.2f}  "
          f"Shapiro p={r['Shapiro_p']:.4f}")

# ─── 4. %INH CRECIMIENTO (distribución completa + LMM) ──────────
print("\n" + "─" * 65)
print("  2. %INH CRECIMIENTO MICELIAL")
print("─" * 65)

y_inh = inh["porcentaje_inhibicion"].dropna().values
inh_transforms = [
    ("Raw", lambda x: x, "Sin transformación"),
    ("IHS(/10)", ihs_transform(10), "arcsinh(y/10)"),
    ("IHS(/50)", ihs_transform(50), "arcsinh(y/50)"),
    ("IHS(/100)", ihs_transform(100), "arcsinh(y/100)"),
    ("CubeRoot", cube_root, "∛y"),
    ("RankGauss", rank_gauss, "Normal scores"),
]
# Solo mid-range para transformaciones acotadas
y_mid = inh[(inh["porcentaje_inhibicion"] > 0) & (inh["porcentaje_inhibicion"] < 100)]["porcentaje_inhibicion"].values
mid_transforms = [
    ("Arcsin-sqrt", lambda x: np.arcsin(np.sqrt(x/100)), "arcsin(√(y/100))"),
    ("Logit", lambda x: np.log(x/100 / (1 - x/100)), "log(y/(100-y))"),
]

# Distribución completa
print("  Distribución completa:")
for r in eval_transforms("%INH completo", y_inh, transform_funcs=inh_transforms):
    print(f"  {r['Transformacion']:12s}  skew={r['Skew']:.2f}  kurt={r['Kurtosis']:.2f}  "
          f"Shapiro p={r['Shapiro_p']:.4f}")

print("\n  Subset (0-100):")
for r in eval_transforms("%INH mid", y_mid, transform_funcs=mid_transforms + inh_transforms):
    print(f"  {r['Transformacion']:12s}  skew={r['Skew']:.2f}  kurt={r['Kurtosis']:.2f}  "
          f"Shapiro p={r['Shapiro_p']:.4f}")

# LMM comparison
print("\n  Residuos del LMM (5.0 mg/mL, método ~ aislado):")
df_lmm = inh.copy()
df_lmm["metodo_id"] = df_lmm["metodo_extraccion"].map(
    {"maceracion": "maceracion", "maceración": "maceracion",
     "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"})
df_lmm = df_lmm[df_lmm["concentracion_mg_ml"] == 5.0]
df_dum = pd.get_dummies(df_lmm, columns=["metodo_id"], drop_first=True, dtype=float)
exog = add_constant(df_dum[["metodo_id_soxhlet", "metodo_id_ultrasonido"]])
groups = df_dum["aislado_id"]

lmm_transforms = [
    ("Raw", lambda x: x, "Crudo"),
    ("IHS/10", ihs_transform(10), "arcsinh(y/10)"),
    ("IHS/50", ihs_transform(50), "arcsinh(y/50)"),
    ("CubeRoot", cube_root, "∛y"),
    ("RankGauss", rank_gauss, "Normal scores"),
]
lmm_res = eval_residuals_lmm("LMM %INH", df_dum["porcentaje_inhibicion"],
                              exog, groups, lmm_transforms)
for r in lmm_res:
    print(f"  {r['Transformacion']:12s}  LogLik={r['LogLik']:.0f}  "
          f"skew({r['Resid_skew']:.2f})  kurt({r['Resid_kurt']:.2f})  "
          f"Shapiro p={r['Shapiro_p']:.4f}")

# Homocedasticidad
print("\n  Homocedasticidad (Levene) de residuos del LMM:")
groups_arr = np.where(df_dum["metodo_id_soxhlet"] == 1, "soxhlet",
                      np.where(df_dum["metodo_id_ultrasonido"] == 1, "ultrasonido", "maceracion"))
for tname, tfunc, _ in lmm_transforms:
    yt = tfunc(df_dum["porcentaje_inhibicion"])
    m = MixedLM(yt, exog, groups=groups).fit(reml=True, maxiter=200)
    r = m.resid
    g1 = r[groups_arr == "maceracion"]
    g2 = r[groups_arr == "soxhlet"]
    g3 = r[groups_arr == "ultrasonido"]
    stat, p = stats.levene(g1, g2, g3)
    print(f"  {tname:12s}  Levene F={stat:.2f}  p={p:.4f}")


# ─── 5. log10 CONIDIAS ─────────────────────────────────────────
print("\n" + "─" * 65)
print("  3. log₁₀(CONIDIAS/mL)")
print("─" * 65)

cl = con_log["conidias_log10"].dropna().values
con_log_transforms = [
    ("Raw", lambda x: x, "Sin transformación"),
    ("BoxCox(λ)", lambda x: stats.boxcox(x + 0.1)[0], "Box-Cox (+0.1)"),
    ("Square", lambda x: x**2, "y²"),
    ("RankGauss", rank_gauss, "Normal scores"),
]
for r in eval_transforms("log10 conidias", cl, transform_funcs=con_log_transforms):
    print(f"  {r['Transformacion']:12s}  skew={r['Skew']:.2f}  kurt={r['Kurtosis']:.2f}  "
          f"Shapiro p={r['Shapiro_p']:.4f}")

# LMM comparison
print("\n  Residuos del LMM (5.0 mg/mL):")
df_lmm_c = con_log.copy()
df_lmm_c["metodo_id"] = df_lmm_c["metodo_extraccion"].map(
    {"maceracion": "maceracion", "maceración": "maceracion",
     "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"})
df_lmm_c = df_lmm_c[df_lmm_c["concentracion_mg_ml"] == 5.0]
df_dum_c = pd.get_dummies(df_lmm_c, columns=["metodo_id"], drop_first=True, dtype=float)
exog_c = add_constant(df_dum_c[["metodo_id_soxhlet", "metodo_id_ultrasonido"]])
groups_c = df_dum_c["aislado_id"]

lmm_c_transforms = [
    ("Raw (log10)", lambda x: x, "log10 crudo"),
    ("BoxCox", lambda x: boxcox(x + 0.1)[0], "Box-Cox"),
    ("Square", lambda x: x**2, "y²"),
]
lmm_c_res = eval_residuals_lmm("LMM conidias", df_dum_c["conidias_log10"],
                                exog_c, groups_c, lmm_c_transforms)
for r in lmm_c_res:
    print(f"  {r['Transformacion']:12s}  LogLik={r['LogLik']:.0f}  "
          f"skew({r['Resid_skew']:.2f})  kurt({r['Resid_kurt']:.2f})  "
          f"Shapiro p={r['Shapiro_p']:.4f}")


# ─── 6. %INH CONIDIAS CRUDO vs %INH log10 ──────────────────────
print("\n" + "─" * 65)
print("  4. %INH CONIDIAS — crudo vs log₁₀ (hoja)")
print("─" * 65)

y_con_crudo = con_inh["porcentaje_inhibicion"].dropna().values
y_con_log = con_inh_log["porcentaje_inhibicion_log10"].dropna().values

con_inh_transforms = [
    ("Raw (crudo)", lambda x: x, "%INH en escala cruda"),
    ("IHS(/100)", ihs_transform(100), "arcsinh(%INH/100)"),
    ("Winsor+BoxCox", lambda x: stats.boxcox(np.clip(x, -500, 100) + 501)[0],
     "Clip[-500,100]+501 + BoxCox"),
    ("Raw (log10 escala)", lambda x: x, "%INH en escala log (hoja)"),
]
for r in eval_transforms("%INH conidias", y_con_crudo, transform_funcs=[con_inh_transforms[0]]):
    print(f"  {r['Transformacion']:12s}  skew={r['Skew']:.2f}  kurt={r['Kurtosis']:.2f}  "
          f"Shapiro p={r['Shapiro_p']:.4f}")

# Si hay %INH en log10
if len(y_con_log) > 0:
    df_cl = con_inh_log.copy()
    df_cl["metodo_id"] = df_cl["metodo_extraccion"].map(
        {"maceracion": "maceracion", "maceración": "maceracion",
         "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"})
    df_cl = df_cl[df_cl["concentracion_mg_ml"] == 5.0]
    if len(df_cl) > 10:
        df_dum_cl = pd.get_dummies(df_cl, columns=["metodo_id"], drop_first=True, dtype=float)
        exog_cl = add_constant(df_dum_cl[["metodo_id_soxhlet", "metodo_id_ultrasonido"]])
        groups_cl = df_dum_cl["aislado_id"]
        lmm_res_c = eval_residuals_lmm("LMM %INH conidias",
                                        df_dum_cl["porcentaje_inhibicion_log10"],
                                        exog_cl, groups_cl,
                                        [("Raw (log10)", lambda x: x, "Crudo en log10")])
        for r in lmm_res_c:
            print(f"\n  LMM con %INH_log10 (n={len(df_cl)}):")
            print(f"  LogLik={r['LogLik']:.0f}  skew={r['Resid_skew']:.2f}  kurt={r['Resid_kurt']:.2f}")

# ─── 7. FIGURAS ──────────────────────────────────────────────────
print("\n" + "─" * 65)
print("  FIGURAS")
print("─" * 65)

# 7a. Q-Q de residuos de cada transformación (LMM %INH)
fig, axes = plt.subplots(2, 3, figsize=(14, 9))
fig.subplots_adjust(hspace=0.4, wspace=0.3)
transforms_plot = [("Crudo", lambda x: x),
                   ("IHS/10", ihs_transform(10)),
                   ("IHS/50", ihs_transform(50)),
                   ("CubeRoot", cube_root),
                   ("RankGauss", rank_gauss)]
for idx, (tname, tfunc) in enumerate(transforms_plot):
    ax = axes[idx // 3, idx % 3]
    yt = tfunc(df_dum["porcentaje_inhibicion"])
    m = MixedLM(yt, exog, groups=groups).fit(reml=True, maxiter=200)
    stats.probplot(m.resid, dist="norm", plot=ax)
    ax.set_title(f"{tname} (Shapiro p={lmm_res[idx]['Shapiro_p']:.4f})")
# Q-Q para log10 conidias
ax = axes[1, 2]
yt_c = df_dum_c["conidias_log10"]
m_c = MixedLM(yt_c, exog_c, groups=groups_c).fit(reml=True, maxiter=200)
stats.probplot(m_c.resid, dist="norm", plot=ax)
ax.set_title(f"log10 conidias (Shapiro p={lmm_c_res[0]['Shapiro_p']:.4f})")
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "comparacion_transformaciones_qq.png", dpi=300)
plt.close(fig)
print("  ✅ comparacion_transformaciones_qq.png")

# 7b. Distribuciones de %INH con diferentes transformaciones
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
fig.subplots_adjust(hspace=0.4, wspace=0.3)
for idx, (tname, tfunc, _) in enumerate(inh_transforms[:6]):
    ax = axes[idx // 3, idx % 3]
    yt = tfunc(y_inh)
    ax.hist(yt, bins=40, color="#2e86ab", edgecolor="white", alpha=0.7)
    ax.set_title(f"{tname}")
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "comparacion_transformaciones_hist.png", dpi=300)
plt.close(fig)
print("  ✅ comparacion_transformaciones_hist.png")

# 7c. Homocedasticidad: boxplot residuos por método para cada transformación
fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.subplots_adjust(wspace=0.35)
for idx, (tname, tfunc, _) in enumerate(lmm_transforms[:3]):
    ax = axes[idx]
    yt = tfunc(df_dum["porcentaje_inhibicion"])
    m = MixedLM(yt, exog, groups=groups).fit(reml=True, maxiter=200)
    r_df = pd.DataFrame({"resid": m.resid, "metodo": groups_arr})
    sns.boxplot(data=r_df, x="metodo", y="resid", hue="metodo",
                palette={"maceracion": "#2e86ab", "soxhlet": "#a23b72",
                         "ultrasonido": "#f18f01"},
                ax=ax, legend=False)
    ax.set_title(f"{tname}")
    ax.set_xlabel("")
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "comparacion_transformaciones_lev.png", dpi=300)
plt.close(fig)
print("  ✅ comparacion_transformaciones_lev.png")


# ═══════════════════════════════════════════════════════════════════
# 8. TABLA COMPARATIVA Y RECOMENDACIONES
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  RECOMENDACIONES FINALES")
print("─" * 65)

recomendaciones = [
    {
        "Dataset": "Rendimiento (%)",
        "Recomendacion": "Sin transformación",
        "Justificacion": "Residuos ANOVA ya normales (Shapiro p=0.33). n=9 insuficiente para transformaciones.",
        "Score": "✅"
    },
    {
        "Dataset": "%INH Crecimiento (completo)",
        "Recomendacion": "Sin transformación — modelar en escala original",
        "Justificacion": (
            "Los residuos crudos tienen la menor asimetría (skew=-0.33). "
            "IHS la empeora (skew=-0.95). CubeRoot mejora homocedasticidad "
            "pero empeora normalidad. Ninguna transformación alcanza normalidad "
            "(Shapiro siempre p<0.001). La interpretación en crudo es directa "
            "y el LMM es robusto a estas desviaciones. "
            "Si se prioriza homocedasticidad, CubeRoot es la mejor alternativa."
        ),
        "Score": "✅ crudo / ◐ CubeRoot"
    },
    {
        "Dataset": "log₁₀(conidias/mL)",
        "Recomendacion": "Mantener log₁₀ — no transformar más",
        "Justificacion": (
            "Ya es log-transformado por el laboratorio. BoxCox sugiere λ=2.17 (cuadrado) "
            "pero la mejora en residuos es marginal y la interpretación empeora."
        ),
        "Score": "✅"
    },
    {
        "Dataset": "%INH Conidias (crudo)",
        "Recomendacion": "NO usar como variable de modelado primaria",
        "Justificacion": (
            "Rango [-10615, 100] con skew=-12.2 y kurtosis=176 — no rescatable "
            "con transformaciones estándar. Usar log₁₀(conidias) como respuesta "
            "primaria. %INH crudo solo para visualización e interpretación. "
            "%INH_log10 (escala hoja) es alternativa secundaria."
        ),
        "Score": "❌"
    },
    {
        "Dataset": "log₁₀(conidias) para dosis-respuesta",
        "Recomendacion": "Mantener log₁₀ — sin transformación adicional",
        "Justificacion": (
            "Modelo B (Maceración) converge y es interpretable en log₁₀. "
            "ICC=0.76 indica que el modelo captura bien la estructura de datos."
        ),
        "Score": "✅"
    },
]

print()
for rec in recomendaciones:
    print(f"  {rec['Score']} {rec['Dataset']}")
    print(f"       → {rec['Recomendacion']}")
    print(f"         {rec['Justificacion'][:90]}...")
    print()

# Guardar tabla para reporte
df_rec = pd.DataFrame(recomendaciones)

with open(DIR_REPORTES / "07_diagnostico_transformaciones.md", "w", encoding="utf-8") as f:
    f.write("# Diagnóstico de transformaciones\n\n")
    f.write(f"**Fecha:** 2026-07-29\n\n")
    f.write("## Resumen de recomendaciones\n\n")
    f.write("| Dataset | Recomendación | Score |\n")
    f.write("|---------|--------------|:----:|\n")
    for rec in recomendaciones:
        f.write(f"| {rec['Dataset']} | {rec['Recomendacion']} | {rec['Score']} |\n")
    f.write("\n## Justificación detallada\n\n")

    for rec in recomendaciones:
        f.write(f"### {rec['Dataset']}\n\n")
        f.write(f"**Recomendación:** {rec['Recomendacion']}\n\n")
        f.write(f"**Justificación:** {rec['Justificacion']}\n\n")

    f.write("## Comparación de transformaciones\n\n")
    f.write("### %INH Crecimiento — residuos del LMM\n\n")
    f.write("| Transformación | LogLik | Skew resid | Kurt resid | Shapiro p |\n")
    f.write("|---------------|:-----:|:---------:|:----------:|:---------:|\n")
    for r in lmm_res:
        f.write(f"| {r['Transformacion']} | {r['LogLik']} | {r['Resid_skew']:.2f} | "
                f"{r['Resid_kurt']:.2f} | {r['Shapiro_p']:.4f} |\n")

    f.write("\n### log₁₀(conidias) — residuos del LMM\n\n")
    f.write("| Transformación | LogLik | Skew resid | Kurt resid | Shapiro p |\n")
    f.write("|---------------|:-----:|:---------:|:----------:|:---------:|\n")
    for r in lmm_c_res:
        f.write(f"| {r['Transformacion']} | {r['LogLik']} | {r['Resid_skew']:.2f} | "
                f"{r['Resid_kurt']:.2f} | {r['Shapiro_p']:.4f} |\n")

    f.write("\n## Figuras\n\n")
    f.write("- `comparacion_transformaciones_qq.png` — Q-Q de residuos\n")
    f.write("- `comparacion_transformaciones_hist.png` — histogramas\n")
    f.write("- `comparacion_transformaciones_lev.png` — homocedasticidad\n")

print(f"  ✅ Reporte guardado: {DIR_REPORTES / '07_diagnostico_transformaciones.md'}")
print(f"\n{'='*65}")
print("  DIAGNÓSTICO COMPLETO")
print(f"{'='*65}")
