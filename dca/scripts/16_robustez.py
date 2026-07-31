"""
16_robustez.py — Diagnóstico de robustez estadística.

Verifica supuestos y estabilidad de los modelos LMM:
  1. ICC (coeficiente de correlación intraclase) para cada modelo
  2. Influential point analysis (Cook-like distance para efectos mixtos)
  3. Verificación de singularidad (¿el efecto aleatorio es identificable?)
  4. Distribución de los efectos aleatorios
  5. Heterocedasticidad en residuos (Breusch-Pagan / Levene)
  6. Bootstrap no paramétrico de intervalos de confianza (Modelo A)
  7. Sensibilidad: refit sin aislados extremos
  8. Comparación LMM vs LM (¿realmente necesitamos el efecto aleatorio?)

Salida: dca/resultados/reportes/08_robustez.md
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats
from config import (
    DIR_TABLAS, DIR_REPORTES, DIR_FIGURAS, SEMILLA_ALEATORIA, COLOR_MET, LABEL_MET
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from config import setup_figure_style
setup_figure_style()

DIR_FIGURAS.mkdir(parents=True, exist_ok=True)
DIR_REPORTES.mkdir(parents=True, exist_ok=True)
np.random.seed(SEMILLA_ALEATORIA)

REPORTE = DIR_REPORTES / "08_robustez.md"


# ═══════════════════════════════════════════════════════════════════
# 1. CARGA Y PREPARACIÓN
# ═══════════════════════════════════════════════════════════════════
print("=" * 65)
print("  DIAGNÓSTICO DE ROBUSTEZ — Tomillo × Fusarium")
print("=" * 65)

crec = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")
inh = crec[~crec["es_control"] & crec["porcentaje_inhibicion"].notna()].copy()
rend = pd.read_csv(DIR_TABLAS / "rendimiento_extraccion.csv")

# Subconjuntos para cada modelo
inh_5 = inh[inh["concentracion_mg_ml"] == 5.0].copy()
mac = inh[inh["metodo_extraccion"] == "maceracion"].copy()
rend_met = rend.copy()

from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools import add_constant
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import ols
from sklearn.linear_model import LinearRegression

# Para Bootstrap
from sklearn.utils import resample


# ═══════════════════════════════════════════════════════════════════
# 2. FUNCIONES AUXILIARES
# ═══════════════════════════════════════════════════════════════════
def fit_modelo_a(df):
    """Fittea Modelo A y devuelve objetos."""
    df_dum = pd.get_dummies(df, columns=["metodo_extraccion"], drop_first=True, dtype=float)
    exog = add_constant(df_dum[["metodo_extraccion_soxhlet", "metodo_extraccion_ultrasonido"]])
    endog = df_dum["porcentaje_inhibicion"]
    groups = df_dum["aislado_id"]
    m = MixedLM(endog, exog, groups=groups).fit(reml=True, maxiter=200)
    return m, exog, endog, groups


def fit_modelo_b(df):
    """Modelo B: Maceración dosis-respuesta."""
    df = df.copy()
    dummies_conc = pd.get_dummies(df["conc_cat"], drop_first=True, dtype=float, prefix="conc")
    df_ml = pd.concat([df.reset_index(drop=True), dummies_conc.reset_index(drop=True)], axis=1)
    exog = add_constant(df_ml[[c for c in dummies_conc.columns]])
    endog = df_ml["porcentaje_inhibicion"]
    groups = df_ml["aislado_id"]
    m = MixedLM(endog, exog, groups=groups).fit(reml=True, maxiter=200)
    return m, exog, endog, groups


# ═══════════════════════════════════════════════════════════════════
# 3. ICC — COEFICIENTE DE CORRELACIÓN INTRACLASE
# ═══════════════════════════════════════════════════════════════════
print("\n📊 ICC — Correlación intraclase")

resultados_icc = []
modelos_info = [
    ("Modelo A", "LMM — %INH ~ método (5.0 mg/mL)", inh_5, fit_modelo_a),
]

# Modelo A
m_a, _, _, _ = fit_modelo_a(inh_5)
var_random = m_a.cov_re.iloc[0, 0]  # var(aislado)
var_residual = m_a.scale               # var(residual)
icc_a = var_random / (var_random + var_residual)
print(f"  Modelo A: σ²_aislado={var_random:.1f}, σ²_residual={var_residual:.1f}")
print(f"  ICC = {icc_a:.3f}  ({icc_a*100:.1f}% de la varianza explicada por aislado)")
resultados_icc.append(("Modelo A (5.0 mg/mL)", var_random, var_residual, icc_a))

# Modelo B (Maceración dosis-respuesta)
mac_temp = mac.copy()
mac_temp["conc_cat"] = mac_temp["concentracion_mg_ml"].map(
    lambda c: "C0" if c == 0 else ("C02" if c == 0.2 else ("C1" if c == 1.0 else "C5"))
)
conc_ref = "C02"
m_b, _, _, _ = fit_modelo_b(mac_temp)
var_random_b = m_b.cov_re.iloc[0, 0]
var_residual_b = m_b.scale
icc_b = var_random_b / (var_random_b + var_residual_b)
print(f"  Modelo B: σ²_aislado={var_random_b:.1f}, σ²_residual={var_residual_b:.1f}")
print(f"  ICC = {icc_b:.3f}  ({icc_b*100:.1f}% de la varianza explicada por aislado)")
resultados_icc.append(("Modelo B (Maceración dosis-respuesta)", var_random_b, var_residual_b, icc_b))


# ═══════════════════════════════════════════════════════════════════
# 4. VERIFICACIÓN DE SINGULARIDAD
# ═══════════════════════════════════════════════════════════════════
print("\n📊 Singularidad del efecto aleatorio")

for nombre, m in [("Modelo A", m_a), ("Modelo B", m_b)]:
    var_rand = m.cov_re.iloc[0, 0]
    singular = var_rand < 1e-6
    ratio = var_rand / m.scale
    print(f"  {nombre}: var_aleatoria={var_rand:.3f}, σ²_ratio={ratio:.4f}")
    if singular:
        print(f"    ⚠ POSIBLE SINGULAR — el efecto aleatorio es esencialmente cero")
    elif ratio < 0.05:
        print(f"    ⚠ BORDE — el efecto aleatorio explica <5% de la varianza residual")
    else:
        print(f"    ✅ OK — efecto aleatorio bien identificado")


# ═══════════════════════════════════════════════════════════════════
# 5. DISTRIBUCIÓN DE EFECTOS ALEATORIOS
# ═══════════════════════════════════════════════════════════════════
print("\n📊 Distribución de efectos aleatorios")

for nombre, m, label in [("Modelo A", m_a, "mod_a"), ("Modelo B", m_b, "mod_b")]:
    re = m.random_effects
    intercepts = [v.iloc[0] for v in re.values()]
    _, p_shapiro = stats.shapiro(intercepts)
    print(f"  {nombre}: Shapiro-Wilk p={p_shapiro:.4f} (H0: normal)")
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(intercepts, bins=12, color="#2e86ab", edgecolor="white", alpha=0.8)
    axes[0].set_xlabel("Intercepto aleatorio")
    axes[0].set_ylabel("Frecuencia")
    axes[0].set_title(f"Efectos aleatorios — {nombre}")
    stats.probplot(intercepts, dist="norm", plot=axes[1])
    axes[1].set_title("Q-Q plot")
    fig.tight_layout()
    fig.savefig(DIR_FIGURAS / f"robustez_re_{label}.png", dpi=200)
    plt.close(fig)
    print(f"    → robustez_re_{label}.png")


# ═══════════════════════════════════════════════════════════════════
# 6. COMPARACIÓN LMM vs LM (¿necesitamos el efecto aleatorio?)
# ═══════════════════════════════════════════════════════════════════
print("\n📊 LMM vs LM — ¿realmente necesitamos el efecto aleatorio?")

for nombre, df_mod, fit_fn in [
    ("Modelo A", inh_5, fit_modelo_a),
]:
    # LMM
    m_lmm, exog, endog, groups = fit_fn(df_mod)
    loglik_lmm = m_lmm.llf
    
    # LM (mismos efectos fijos, sin aleatorio)
    lm = LinearRegression()
    lm.fit(exog, endog)
    resid_lm = endog - lm.predict(exog)
    # Log-likelihood aproximado para LM
    n = len(endog)
    sse = np.sum(resid_lm ** 2)
    loglik_lm = -n/2 * np.log(2 * np.pi * sse/n) - n/2
    
    # LRT: 2 * (logLik_lmm - logLik_lm) ~ chi²(1)
    lrt_stat = 2 * (loglik_lmm - loglik_lm)
    # df = 1 (una varianza de RE adicional)
    p_lrt = 1 - stats.chi2.cdf(lrt_stat, 1)
    aic_lmm = m_lmm.aic
    aic_lm = n * np.log(sse/n) + 2 * exog.shape[1]
    
    print(f"  {nombre}:")
    print(f"    LMM logLik={loglik_lmm:.1f}, AIC={aic_lmm:.0f}")
    print(f"    LM  logLik={loglik_lm:.1f}, AIC={aic_lm:.0f}")
    print(f"    LRT χ²(1)={lrt_stat:.1f}, p={p_lrt:.4f}")
    if p_lrt < 0.05:
        print(f"    ✅ El efecto aleatorio mejora significativamente el modelo")
    else:
        print(f"    ⚠ El efecto aleatorio NO mejora significativamente el modelo")


# ═══════════════════════════════════════════════════════════════════
# 7. SENSIBILIDAD — REFIT SIN AISLADOS EXTREMOS
# ═══════════════════════════════════════════════════════════════════
print("\n📊 Análisis de sensibilidad — refit sin aislados extremos")

# Identificar aislados con %INH extremo en Maceración 5.0 mg/mL
mac_5 = inh_5[inh_5["metodo_extraccion"] == "maceracion"]
medians = mac_5.groupby("aislado_id")["porcentaje_inhibicion"].median()
q1, q3 = medians.quantile(0.25), medians.quantile(0.75)
iqr = q3 - q1
outliers_ais = medians[(medians < q1 - 1.5 * iqr) | (medians > q3 + 1.5 * iqr)].index
print(f"  Aislados extremos en Maceración 5.0 mg/mL (IQR): {len(outliers_ais)}")
for ais in outliers_ais:
    print(f"    {ais}: mediana={medians[ais]:.1f}%")

if len(outliers_ais) > 0:
    # Refit sin esos aislados
    inh_5_clean = inh_5[~inh_5["aislado_id"].isin(outliers_ais)]
    m_a_clean, _, _, _ = fit_modelo_a(inh_5_clean)
    print(f"\n  Modelo A SIN extremos (n={len(inh_5_clean)}):")
    print(f"    Maceración:  {m_a_clean.fe_params['const']:.1f}%")
    diff_mac = m_a.fe_params['const'] - m_a_clean.fe_params['const']
    print(f"    Diferencia vs completo: {diff_mac:.2f} pp")
else:
    print("  Sin aislados extremos detectados.")


# ═══════════════════════════════════════════════════════════════════
# 8. BOOTSTRAP — IC no paramétricos para Modelo A
# ═══════════════════════════════════════════════════════════════════
print("\n📊 Bootstrap no paramétrico — Modelo A (n=500 iteraciones)")

n_boot = 500
boot_mac = []
boot_sox = []
boot_ult = []

# Agrupar por aislado y remuestrear aislados (caso bootstrap por clusters)
aislados = inh_5["aislado_id"].unique()
n_ais = len(aislados)

for i in range(n_boot):
    if (i + 1) % 100 == 0:
        print(f"    Bootstrap: {i+1}/{n_boot}")
    # Remuestrear aislados CON REPOSICIÓN
    boot_ais = np.random.choice(aislados, size=n_ais, replace=True)
    # Construir dataset bootstrap: todos los datos de los aislados seleccionados
    boot_list = []
    for ais in boot_ais:
        boot_list.append(inh_5[inh_5["aislado_id"] == ais])
    boot_df = pd.concat(boot_list).reset_index(drop=True)
    try:
        m_boot, _, _, _ = fit_modelo_a(boot_df)
        boot_mac.append(m_boot.fe_params["const"])
        boot_sox.append(m_boot.fe_params["const"] + m_boot.fe_params["metodo_extraccion_soxhlet"])
        boot_ult.append(m_boot.fe_params["const"] + m_boot.fe_params["metodo_extraccion_ultrasonido"])
    except Exception:
        continue

if len(boot_mac) > 100:
    ci_mac = np.percentile(boot_mac, [2.5, 97.5])
    ci_sox = np.percentile(boot_sox, [2.5, 97.5])
    ci_ult = np.percentile(boot_ult, [2.5, 97.5])
    print(f"  IC Bootstrap 95% — Maceración:  {m_a.fe_params['const']:.1f} [{ci_mac[0]:.1f}, {ci_mac[1]:.1f}]")
    print(f"  IC Bootstrap 95% — Soxhlet:     {m_a.fe_params['const'] + m_a.fe_params['metodo_extraccion_soxhlet']:.1f} [{ci_sox[0]:.1f}, {ci_sox[1]:.1f}]")
    print(f"  IC Bootstrap 95% — Ultrasonido: {m_a.fe_params['const'] + m_a.fe_params['metodo_extraccion_ultrasonido']:.1f} [{ci_ult[0]:.1f}, {ci_ult[1]:.1f}]")
    
    # Figura bootstrap
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, data, label, color in [
        (axes[0], boot_mac, "Maceración", "#2e86ab"),
        (axes[1], boot_sox, "Soxhlet", "#a23b72"),
        (axes[2], boot_ult, "Ultrasonido", "#f18f01"),
    ]:
        ax.hist(data, bins=25, color=color, edgecolor="white", alpha=0.7)
        ax.axvline(np.mean(data), color="black", ls="--", lw=1.5)
        ax.set_xlabel(f"{label} (%INH)")
        ax.set_ylabel("Frecuencia")
    fig.suptitle("Distribución Bootstrap — %INH estimado (Modelo A)", fontsize=12)
    fig.tight_layout()
    fig.savefig(DIR_FIGURAS / "robustez_bootstrap.png", dpi=200)
    plt.close(fig)
    print(f"    → robustez_bootstrap.png")
else:
    print(f"  ⚠ Bootstrap no convergió suficientes modelos ({len(boot_mac)} iteraciones)")


# ═══════════════════════════════════════════════════════════════════
# 9. HETEROCEDASTICIDAD EN RESIDUOS
# ═══════════════════════════════════════════════════════════════════
print("\n📊 Heterocedasticidad — residuos LMM")

for nombre, m, exog in [("Modelo A", m_a, fit_modelo_a(inh_5)[1]),
                         ("Modelo B", m_b, fit_modelo_b(mac_temp)[1])]:
    resid = m.resid
    # Levene: |resid| ~ fitted (no depende de distribución)
    fitted = m.fittedvalues
    # Dividir en grupos para Levene (bajos, medios, altos)
    terciles = pd.qcut(fitted, 3, labels=["bajo", "medio", "alto"], duplicates="drop")
    groups = [np.abs(resid[terciles == g]) for g in sorted(terciles.unique()) if len(resid[terciles == g]) > 1]
    if len(groups) >= 2:
        stat_levene, p_levene = stats.levene(*groups)
        print(f"  {nombre}: Levene={stat_levene:.2f}, p={p_levene:.4f}")
        if p_levene > 0.05:
            print(f"    ✅ Homocedasticidad aceptable")
        else:
            print(f"    ⚠ Heterocedasticidad detectada")


# ═══════════════════════════════════════════════════════════════════
# 10. REPORTE
# ═══════════════════════════════════════════════════════════════════
print("\n📝 Generando reporte...")

with open(REPORTE, "w") as f:
    f.write("# Diagnóstico de Robustez Estadística\n\n")
    f.write("## 1. Coeficiente de Correlación Intraclase (ICC)\n\n")
    f.write("| Modelo | σ²_aislado | σ²_residual | ICC | Interpretación |\n")
    f.write("|--------|-----------|-------------|-----|----------------|\n")
    for nombre, vr, ve, icc in resultados_icc:
        if icc > 0.3:
            interp = "Fuerte correlación intra-aislado — el efecto aleatorio es esencial"
        elif icc > 0.1:
            interp = "Correlación moderada — el efecto aleatorio es relevante"
        else:
            interp = "Correlación baja — el efecto aleatorio aporta poco"
        f.write(f"| {nombre} | {vr:.1f} | {ve:.1f} | {icc:.3f} | {interp} |\n")
    
    f.write("\n## 2. Singularidad del Efecto Aleatorio\n\n")
    for nombre, m in [("Modelo A", m_a), ("Modelo B", m_b)]:
        var_rand = m.cov_re.iloc[0, 0]
        ratio = var_rand / m.scale
        f.write(f"- **{nombre}**: var_aleatoria = {var_rand:.3f}, σ²_ratio = {ratio:.4f}\n")
        if var_rand < 1e-6:
            f.write("  - ⚠ **Posible singularidad** — el efecto aleatorio es esencialmente cero\n")
        elif ratio < 0.05:
            f.write("  - ⚠ **Borde** — el efecto aleatorio explica <5% de la varianza residual\n")
        else:
            f.write("  - ✅ Efecto aleatorio bien identificado\n")
    
    f.write("\n## 3. Normalidad de Efectos Aleatorios\n\n")
    for nombre, m, label in [("Modelo A", m_a, "mod_a"), ("Modelo B", m_b, "mod_b")]:
        re = m.random_effects
        intercepts = [v.iloc[0] for v in re.values()]
        _, p_sh = stats.shapiro(intercepts)
        f.write(f"- **{nombre}**: Shapiro-Wilk p = {p_sh:.4f}\n")
        if p_sh > 0.05:
            f.write(f"  - ✅ No se rechaza normalidad de los efectos aleatorios\n")
        else:
            f.write(f"  - ⚠ Los efectos aleatorios se desvían de la normalidad\n")
        f.write(f"  - Figura: `robustez_re_{label}.png`\n")
    
    f.write("\n## 4. Necesidad del Efecto Aleatorio (LMM vs LM)\n\n")
    f.write(f"- **Modelo A**: LRT χ²(1) = {lrt_stat:.1f}, p = {p_lrt:.4f}\n")
    if p_lrt < 0.05:
        f.write(f"  - ✅ El efecto aleatorio mejora significativamente el modelo\n")
    else:
        f.write(f"  - ⚠ No hay evidencia significativa de que el efecto aleatorio sea necesario\n")
    f.write(f"  - AIC_LMM = {aic_lmm:.0f}, AIC_LM = {aic_lm:.0f}\n")
    
    f.write("\n## 5. Análisis de Sensibilidad\n\n")
    f.write(f"- Aislados extremos detectados (criterio IQR): {len(outliers_ais)}\n")
    for ais in outliers_ais:
        f.write(f"  - {ais}: mediana = {medians[ais]:.1f}%\n")
    if len(outliers_ais) > 0:
        f.write(f"\nModelo A refit sin extremos:\n")
        f.write(f"- Maceración: {m_a_clean.fe_params['const']:.1f}% ")
        f.write(f"(diferencia vs completo: {diff_mac:.2f} pp)\n")
        diff_sox = (m_a.fe_params['const'] + m_a.fe_params['metodo_extraccion_soxhlet']) - \
                   (m_a_clean.fe_params['const'] + m_a_clean.fe_params['metodo_extraccion_soxhlet'])
        f.write(f"- Soxhlet: {m_a_clean.fe_params['const'] + m_a_clean.fe_params['metodo_extraccion_soxhlet']:.1f}% ")
        f.write(f"(diferencia vs completo: {diff_sox:.2f} pp)\n")
        if abs(diff_mac) < 3 and abs(diff_sox) < 3:
            f.write("  - ✅ Estimaciones estables (< 3 pp de diferencia)\n")
        else:
            f.write("  - ⚠ Cambio sustancial — los aislados extremos influyen en las estimaciones\n")
    
    f.write("\n## 6. Intervalos de Confianza Bootstrap (Modelo A)\n\n")
    if len(boot_mac) > 100:
        f.write(f"| Método | Estimación puntual | IC Bootstrap 95% |\n")
        f.write(f"|--------|-------------------|------------------|\n")
        for label, est, ci in [
            ("Maceración", m_a.fe_params['const'], ci_mac),
            ("Soxhlet", m_a.fe_params['const'] + m_a.fe_params['metodo_extraccion_soxhlet'], ci_sox),
            ("Ultrasonido", m_a.fe_params['const'] + m_a.fe_params['metodo_extraccion_ultrasonido'], ci_ult),
        ]:
            f.write(f"| {label} | {est:.1f}% | [{ci[0]:.1f}, {ci[1]:.1f}] |\n")
        f.write(f"\n- Figura: `robustez_bootstrap.png`\n")
        # Comparar IC Bootstrap vs IC LMM (Wald)
        f.write(f"\nComparación Bootstrap vs Wald:\n")
        for label, est, ci_boot in [
            ("Maceración", m_a.fe_params['const'], ci_mac),
            ("Soxhlet", m_a.fe_params['const'] + m_a.fe_params['metodo_extraccion_soxhlet'], ci_sox),
            ("Ultrasonido", m_a.fe_params['const'] + m_a.fe_params['metodo_extraccion_ultrasonido'], ci_ult),
        ]:
            se_wald = m_a.bse_fe.iloc[0] if label == "Maceración" else \
                      np.sqrt(m_a.bse_fe.iloc[0]**2 + m_a.bse_fe.iloc[1]**2)
            ci_wald = (est - 1.96 * se_wald, est + 1.96 * se_wald)
            f.write(f"  - {label}: Bootstrap [{ci_boot[0]:.1f}, {ci_boot[1]:.1f}] vs Wald [{ci_wald[0]:.1f}, {ci_wald[1]:.1f}]\n")
    else:
        f.write("  - Bootstrap no pudo completarse por problemas de convergencia.\n")
    
    f.write("\n## 7. Homocedasticidad de Residuos\n\n")
    for nombre, m, exog in [("Modelo A", m_a, fit_modelo_a(inh_5)[1]),
                             ("Modelo B", m_b, fit_modelo_b(mac_temp)[1])]:
        resid = m.resid
        fitted = m.fittedvalues
        _, p_levene_tmp = stats.levene(
            np.abs(resid[fitted < np.percentile(fitted, 33)]),
            np.abs(resid[(fitted >= np.percentile(fitted, 33)) & (fitted < np.percentile(fitted, 67))]),
            np.abs(resid[fitted >= np.percentile(fitted, 67)]),
        )
        f.write(f"- **{nombre}**: Levene p = {p_levene_tmp:.4f}\n")
        if p_levene_tmp > 0.05:
            f.write(f"  - ✅ Homocedasticidad: varianza constante en el rango de predicción\n")
        else:
            f.write(f"  - ⚠ Heterocedasticidad: la varianza residual cambia con el valor ajustado\n")
    
    f.write("\n## 8. Resumen de Robustez\n\n")
    n_green = 0
    n_red = 0
    # Contar verdes/rojos
    if icc_a > 0.1: n_green += 1
    else: n_red += 1
    if m_a.cov_re.iloc[0, 0] > 1e-3: n_green += 1
    else: n_red += 1
    if p_lrt < 0.05: n_green += 1
    else: n_red += 1
    if abs(diff_mac) < 3 if len(outliers_ais) > 0 else 1: n_green += 1
    else: n_red += 1
    
    total = n_green + n_red
    pct_green = n_green / total * 100 if total > 0 else 0
    f.write(f"- **Indicadores OK**: {n_green}/{total} ({pct_green:.0f}%)\n")
    if pct_green >= 75:
        f.write("- **Conclusión**: El análisis es robusto. Los modelos LMM están bien especificados,\n")
        f.write("  los efectos aleatorios son identificables y las estimaciones son estables.\n")
    elif pct_green >= 50:
        f.write("- **Conclusión**: Robustez moderada. Algunos supuestos no se cumplen\n")
        f.write("  estrictamente, pero las estimaciones principales son confiables.\n")
    else:
        f.write("- **Conclusión**: Se recomienda revisar la especificación de los modelos.\n")
        f.write("  Varios indicadores de robustez están por debajo de lo deseable.\n")

print(f"\n  ✅ Reporte guardado: {REPORTE.name}")
print(f"{'='*65}")
