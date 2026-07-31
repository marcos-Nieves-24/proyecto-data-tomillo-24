"""
09_eda.py — Análisis Exploratorio de Datos (EDA).

Genera figuras y tablas resumen para entender la estructura, distribución
y calidad de los datos antes del modelado inferencial.

Salidas:
  - dca/resultados/figuras/   → gráficos en PNG
  - dca/resultados/tablas/    → tablas resumen en CSV
  - dca/resultados/reportes/  → reporte de EDA en MD
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import numpy as np
from config import DIR_DATOS, DIR_TABLAS, DIR_FIGURAS, DIR_REPORTES, SEMILLA_ALEATORIA, COLOR_MET, LABEL_MET, COLOR_CONC, setup_figure_style

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

np.random.seed(SEMILLA_ALEATORIA)

# ─── Configuración estética ─────────────────────────────────────────
setup_figure_style(dpi=300)
COLOR_MET = COLOR_MET  # alias for backward compat in this script
DIR_FIGURAS.mkdir(parents=True, exist_ok=True)


# ═════════════════════════════════════════════════════════════════════
# 1. CARGAR DATOS
# ═════════════════════════════════════════════════════════════════════
print("=" * 65)
print("  ANÁLISIS EXPLORATORIO DE DATOS — Tomillo × Fusarium")
print("=" * 65)

crecimiento = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")
conidias = pd.read_csv(DIR_TABLAS / "conidias.csv")
rendimiento = pd.read_csv(DIR_TABLAS / "rendimiento_extraccion.csv")

CREC = crecimiento.copy()
CONI = conidias.copy()

# Normalizar método en rendimiento (no pasó por step 4)
rendimiento["metodo_extraccion"] = (
    rendimiento["metodo_extraccion"].str.strip().str.lower()
    .str.normalize("NFKD").str.encode("ascii", errors="ignore").str.decode("ascii")
)

# ─── Construir flag para análisis por método ───────────────────────
CREC["solamente_mac"] = CREC["metodo_extraccion"] == "maceracion"

print(f"\n📊 Datos cargados:")
print(f"  Crecimiento micelial: {len(CREC)} filas, {CREC['aislado_id'].nunique()} aislados")
print(f"  Conidias:            {len(CONI)} filas, {CONI['aislado_id'].nunique()} aislados")
print(f"  Rendimiento:         {len(rendimiento)} filas")


# ═════════════════════════════════════════════════════════════════════
# 2. RESUMEN NUMÉRICO
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  2. RESUMEN NUMÉRICO")
print("=" * 65)

def resumen_numerico(df, col, nombre, grupo=None):
    """Tabla resumen con n, media, mediana, desvío, min, max, %NaN."""
    if grupo:
        res = df.groupby(grupo)[col].agg(
            n="count", n_nulo=lambda x: x.isna().sum(), media="mean",
            mediana="median", std="std", min="min", max="max"
        ).reset_index()
    else:
        vals = df[col]
        res = pd.DataFrame([{
            "variable": nombre, "n": int(vals.notna().sum()),
            "n_nulo": int(vals.isna().sum()), "media": vals.mean(),
            "mediana": vals.median(), "std": vals.std(),
            "min": vals.min(), "max": vals.max()
        }])
    return res

# Crecimiento
print("\n── Crecimiento micelial (mm) ──")
print(resumen_numerico(CREC, "crecimiento_mm", "Crecimiento (mm)",
                        grupo=["metodo_extraccion", "concentracion_mg_ml"]).to_string())

# %INH
print("\n── Inhibición (%) ──")
mask_inh = ~CREC["es_control"] & CREC["porcentaje_inhibicion"].notna()
inh_vals = CREC[mask_inh]
print(resumen_numerico(inh_vals, "porcentaje_inhibicion", "Inhibición (%)",
                        grupo=["metodo_extraccion", "concentracion_mg_ml"]).to_string())

# Conidias
print("\n── Conidias (log₁₀/mL) ──")
print(resumen_numerico(CONI, "conidias_log10", "Conidias (log₁₀)",
                        grupo=["metodo_extraccion", "concentracion_mg_ml"]).to_string())

# Rendimiento
print("\n── Rendimiento (%) ──")
print(resumen_numerico(rendimiento, "rendimiento_pct", "Rendimiento (%)",
                        grupo=["metodo_extraccion"]).to_string())


# ═════════════════════════════════════════════════════════════════════
# 3. DISTRIBUCIONES
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  3. GENERANDO GRÁFICOS DE DISTRIBUCIÓN")
print("=" * 65)

# 3a. Histograma: crecimiento mm (tratamiento vs control)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, (label, mask) in zip(axes, [
    ("Control (0 mg/mL)", CREC["es_control"]),
    ("Tratamiento", ~CREC["es_control"])
]):
    datos = CREC.loc[mask & CREC["crecimiento_mm"].notna(), "crecimiento_mm"]
    ax.hist(datos, bins=25, color="#2e86ab" if label == "Control (0 mg/mL)" else "#a23b72",
            edgecolor="white", alpha=0.8)
    ax.set_xlabel("Crecimiento (mm)")
    ax.set_ylabel("Frecuencia")
    ax.set_title(f"Crecimiento micelial — {label}")
    ax.axvline(datos.median(), color="red", ls="--", label=f"Mediana={datos.median():.0f}")
    ax.legend()
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_distribucion_crecimiento.png", dpi=300)
plt.close(fig)
print("  ✅ eda_distribucion_crecimiento.png")

# 3b. Histograma: %INH
fig, ax = plt.subplots(figsize=(10, 5))
datos = inh_vals["porcentaje_inhibicion"]
ax.hist(datos, bins=30, color="#2e86ab", edgecolor="white", alpha=0.8)
ax.set_xlabel("Inhibición (%)")
ax.set_ylabel("Frecuencia")
ax.set_title("Distribución del % de inhibición — todos los métodos")
ax.axvline(datos.median(), color="red", ls="--", label=f"Mediana={datos.median():.1f}%")
ax.axvline(0, color="gray", ls=":", alpha=0.5)
ax.legend()
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_distribucion_inhibicion.png", dpi=300)
plt.close(fig)
print("  ✅ eda_distribucion_inhibicion.png")

# 3c. Histograma: conidias log10
fig, ax = plt.subplots(figsize=(10, 5))
datos = CONI.loc[CONI["conidias_log10"].notna(), "conidias_log10"]
ax.hist(datos, bins=30, color="#a23b72", edgecolor="white", alpha=0.8)
ax.set_xlabel("log₁₀(conidias/mL)")
ax.set_ylabel("Frecuencia")
ax.set_title("Distribución de conidias (log₁₀)")
ax.axvline(datos.median(), color="red", ls="--", label=f"Mediana={datos.median():.2f}")
ax.legend()
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_distribucion_conidias.png", dpi=300)
plt.close(fig)
print("  ✅ eda_distribucion_conidias.png")

# 3d. Rendimiento por método
fig, ax = plt.subplots(figsize=(10, 5.5))
colores = [COLOR_MET[m] for m in rendimiento["metodo_extraccion"]]
ax.bar(range(len(rendimiento)), rendimiento["rendimiento_pct"], color=colores, edgecolor="white")
ax.set_xticks(range(len(rendimiento)))
ax.set_xticklabels([f"{m}\n(rép {r})" for m, r in
                     zip(rendimiento["metodo_extraccion"], rendimiento["replica_biologica"])],
                    rotation=45, ha="right")
ax.set_ylabel("Rendimiento (%)")
ax.set_title("Rendimiento de extracción por método")
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_rendimiento.png", dpi=300)
plt.close(fig)
print("  ✅ eda_rendimiento.png")


# ═════════════════════════════════════════════════════════════════════
# 4. BOXPLOTS
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  4. GENERANDO BOXPLOTS")
print("=" * 65)

# 4a. %INH ~ método
fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=inh_vals, x="metodo_extraccion", y="porcentaje_inhibicion",
            palette=COLOR_MET, ax=ax)
sns.stripplot(data=inh_vals, x="metodo_extraccion", y="porcentaje_inhibicion",
              color="black", alpha=0.15, size=3, ax=ax)
ax.set_xlabel("Método de extracción")
ax.set_ylabel("Inhibición (%)")
ax.set_title("Inhibición por método de extracción")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_boxplot_inhibicion_metodo.png", dpi=300)
plt.close(fig)
print("  ✅ eda_boxplot_inhibicion_metodo.png")

# 4b. %INH ~ concentración (solo MACERACIÓN, que tiene múltiples conc)
fig, ax = plt.subplots(figsize=(8, 5))
mac = inh_vals[inh_vals["metodo_extraccion"] == "maceracion"].copy()
orden_conc = sorted(mac["concentracion_mg_ml"].unique())
sns.boxplot(data=mac, x="concentracion_mg_ml", y="porcentaje_inhibicion",
            order=orden_conc, palette=COLOR_CONC, ax=ax)
sns.stripplot(data=mac, x="concentracion_mg_ml", y="porcentaje_inhibicion",
              order=orden_conc, color="black", alpha=0.15, size=3, ax=ax)
ax.set_xlabel("Concentración (mg/mL)")
ax.set_ylabel("Inhibición (%)")
ax.set_title("Inhibición por concentración — Maceración")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_boxplot_inhibicion_concentracion_mac.png", dpi=300)
plt.close(fig)
print("  ✅ eda_boxplot_inhibicion_concentracion_mac.png")

# 4c. %INH ~ concentración para SOX y ULT (solo 2 niveles)
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
fig.subplots_adjust(wspace=0.3)
for ax, metodo in zip(axes, ["soxhlet", "ultrasonido"]):
    sub = inh_vals[inh_vals["metodo_extraccion"] == metodo]
    sns.boxplot(data=sub, x="concentracion_mg_ml", y="porcentaje_inhibicion",
                color=COLOR_MET[metodo], ax=ax)
    sns.stripplot(data=sub, x="concentracion_mg_ml", y="porcentaje_inhibicion",
                  color="black", alpha=0.15, size=3, ax=ax)
    ax.set_xlabel("Concentración (mg/mL)")
    ax.set_ylabel("Inhibición (%)")
    ax.set_title(f"Inhibición — {metodo}")
    ax.axhline(0, color="gray", ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_boxplot_inhibicion_soxhlet_ultra.png", dpi=300)
plt.close(fig)
print("  ✅ eda_boxplot_inhibicion_soxhlet_ultra.png")

# 4d. %INH ~ aislado (top 10 más susceptibles y top 10 menos)
fig, ax = plt.subplots(figsize=(14, 6))
med_por_aislado = inh_vals.groupby("aislado_id")["porcentaje_inhibicion"].median().sort_values()
top_extremos = list(med_por_aislado.head(10).index) + list(med_por_aislado.tail(10).index)
sub = inh_vals[inh_vals["aislado_id"].isin(top_extremos)]
sns.boxplot(data=sub, x="aislado_id", y="porcentaje_inhibicion",
            order=top_extremos, palette="RdYlBu_r", ax=ax)
ax.set_xlabel("Aislado")
ax.set_ylabel("Inhibición (%)")
ax.set_title("Inhibición por aislado — 10 más y 10 menos susceptibles (mediana)")
ax.tick_params(axis="x", rotation=45)
ax.axhline(0, color="gray", ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_boxplot_inhibicion_aislados_extremos.png", dpi=300)
plt.close(fig)
print("  ✅ eda_boxplot_inhibicion_aislados_extremos.png")

# 4e. Crecimiento del control por aislado
fig, ax = plt.subplots(figsize=(12, 5))
ctrl = CREC[CREC["es_control"] & CREC["crecimiento_mm"].notna()].drop_duplicates(
    subset=["aislado_id", "metodo_extraccion"]
)
sns.boxplot(data=ctrl, x="aislado_id", y="crecimiento_mm", hue="metodo_extraccion",
            palette=COLOR_MET, ax=ax)
ax.set_xlabel("Aislado")
ax.set_ylabel("Crecimiento control (mm)")
ax.set_title("Crecimiento del control por aislado y método")
ax.tick_params(axis="x", rotation=45)
ax.legend(title="Método")
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_control_por_aislado.png", dpi=300)
plt.close(fig)
print("  ✅ eda_control_por_aislado.png")


# ═════════════════════════════════════════════════════════════════════
# 5. MAPA DE DATOS FALTANTES
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  5. MAPA DE DATOS FALTANTES")
print("=" * 65)

def matriz_faltantes(df, nombre, ax):
    """Matriz de valores no nulos por (aislado, método, conc)."""
    matriz = df.pivot_table(
        index="aislado_id", columns=["metodo_extraccion", "concentracion_mg_ml"],
        values=nombre, aggfunc=lambda x: int(x.notna().any())
    )
    sns.heatmap(matriz, cmap="RdYlGn", cbar_kws={"label": "Dato presente"},
                linewidths=0.5, linecolor="white", ax=ax, vmin=0, vmax=1)
    ax.set_title(f"Datos presentes: {nombre}")
    ax.set_xlabel("")
    ax.set_ylabel("Aislado")
    return matriz

fig, axes = plt.subplots(1, 2, figsize=(20, 12))
fig.subplots_adjust(wspace=0.15)
matriz_faltantes(CREC, "crecimiento_mm", axes[0])
matriz_faltantes(CONI, "conidias_log10", axes[1])
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_mapa_faltantes.png", dpi=300)
plt.close(fig)
print("  ✅ eda_mapa_faltantes.png")

# Tabla resumen de faltantes
faltantes_crec = CREC.groupby(["metodo_extraccion", "concentracion_mg_ml"]).agg(
    total=("crecimiento_mm", "count"),
    presentes=("crecimiento_mm", lambda x: int(x.notna().sum())),
    pct=("crecimiento_mm", lambda x: f"{100 * x.notna().sum() / len(x):.0f}%")
).reset_index()
print(faltantes_crec.to_string())
faltantes_crec.to_csv(DIR_TABLAS / "resumen_faltantes_crecimiento.csv", index=False)

faltantes_coni = CONI.groupby(["metodo_extraccion", "concentracion_mg_ml"]).agg(
    total=("conidias_log10", "count"),
    presentes=("conidias_log10", lambda x: int(x.notna().sum())),
    pct=("conidias_log10", lambda x: f"{100 * x.notna().sum() / len(x):.0f}%")
).reset_index()
print(faltantes_coni.to_string())
faltantes_coni.to_csv(DIR_TABLAS / "resumen_faltantes_conidias.csv", index=False)


# ═════════════════════════════════════════════════════════════════════
# 6. INHIBICIÓN COMPLETA Y NEGATIVA
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  6. INHIBICIÓN COMPLETA Y NEGATIVA")
print("=" * 65)

# Inhibición completa (crecimiento = 0)
completa = inh_vals[inh_vals["inhibicion_completa"]]
print(f"\n  Inhibición completa (crecimiento = 0 mm): {len(completa)} casos")
print(f"  Aislados afectados: {sorted(completa['aislado_id'].unique())}")
print(f"  Por método y concentración:")
print(completa.groupby(["metodo_extraccion", "concentracion_mg_ml"]).size().to_string())

# Inhibición negativa (crecimiento > control)
negativa = inh_vals[inh_vals["inhibicion_negativa"]]
print(f"\n  Inhibición negativa (crecimiento > control): {len(negativa)} casos")
print(f"  Aislados afectados: {sorted(negativa['aislado_id'].unique())}")
print(f"  Rango de %INH negativo: {negativa['porcentaje_inhibicion'].min():.1f}% "
      f"a {negativa['porcentaje_inhibicion'].max():.1f}%")
print(f"  Por método y concentración:")
print(negativa.groupby(["metodo_extraccion", "concentracion_mg_ml"]).size().to_string())

# 6a. Barplot: conteo de inhibición completa por aislado
fig, ax = plt.subplots(figsize=(14, 5))
conteo_completa = completa.groupby("aislado_id").size().sort_values(ascending=False)
conteo_completa.plot(kind="bar", ax=ax, color="#2e86ab", edgecolor="white")
ax.set_xlabel("Aislado")
ax.set_ylabel("Casos de inhibición completa")
ax.set_title("Inhibición completa por aislado (crecimiento = 0 mm)")
ax.tick_params(axis="x", rotation=45)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_inhibicion_completa.png", dpi=300)
plt.close(fig)
print("  ✅ eda_inhibicion_completa.png")


# ═════════════════════════════════════════════════════════════════════
# 7. COMPARACIÓN %INH CALCULADO VS HOJA ORIGINAL
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  7. VALIDACIÓN: %INH CALCULADO VS HOJA ORIGINAL")
print("=" * 65)

if "diferencia_con_hoja" in CREC.columns:
    # NOTA: %INH_hoja solo es válido para C1 (5 mg/mL). Para C2 y C3,
    # la columna en el Excel contiene un valor diferente (no es %INH).
    # Ver docs/05_DATA_QUALITY_ISSUES.md → DQ09.
    diffs_validas = CREC[
        (CREC["concentracion_mg_ml"] == 5.0)
        & CREC["diferencia_con_hoja"].notna()
        & ~CREC["es_control"]
    ]["diferencia_con_hoja"]

    print(f"  Comparación solo para C1 (5 mg/mL) — {len(diffs_validas)} observaciones:")
    print(f"    |diferencia|: media={diffs_validas.abs().mean():.3f}, "
          f"max={diffs_validas.abs().max():.3f}, "
          f"mediana={diffs_validas.abs().median():.3f}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(diffs_validas, bins=30, color="#2e86ab", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Diferencia (%INH calculado − %INH hoja)")
    ax.set_ylabel("Frecuencia")
    ax.set_title("Validación: %INH calculado vs. hoja (solo C1 = 5 mg/mL)")
    ax.axvline(0, color="red", ls="--")
    fig.tight_layout()
    fig.savefig(DIR_FIGURAS / "eda_validacion_inh_hoja.png", dpi=300)
    plt.close(fig)
    print("  ✅ eda_validacion_inh_hoja.png (solo C1 = 5 mg/mL)")


# ═════════════════════════════════════════════════════════════════════
# 8. DIAGNÓSTICO DE SUPUESTOS
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  8. DIAGNÓSTICO DE SUPUESTOS")
print("=" * 65)

# 8a. Normalidad: crecimiento (tratamiento, solo MAC)
print("\n── Normalidad (Shapiro-Wilk) ──")
for var_name, var_col, df_src in [
    ("Crecimiento (mm) — tratamiento", "crecimiento_mm", CREC),
    ("Inhibición (%)", "porcentaje_inhibicion", inh_vals),
    ("Conidias (log₁₀)", "conidias_log10", CONI),
]:
    vals = df_src[var_col].dropna()
    if len(vals) > 5000:
        vals = vals.sample(5000, random_state=SEMILLA_ALEATORIA)
    if len(vals) >= 3:
        stat, p = stats.shapiro(vals)
        print(f"  {var_name}: W={stat:.4f}, p={p:.6f} "
              f"{'❌ NO normal' if p < 0.05 else '✅ Normal (α=0.05)'}")

# 8b. QQ plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.subplots_adjust(wspace=0.35)
for ax, (var_name, var_col, df_src, color) in zip(axes, [
    ("Crecimiento (mm)", "crecimiento_mm", CREC, "#2e86ab"),
    ("Inhibición (%)", "porcentaje_inhibicion", inh_vals, "#a23b72"),
    ("Conidias (log₁₀)", "conidias_log10", CONI, "#f18f01"),
]):
    vals = df_src[var_col].dropna()
    if len(vals) > 2000:
        vals = vals.sample(2000, random_state=SEMILLA_ALEATORIA)
    stats.probplot(vals, dist="norm", plot=ax)
    ax.get_lines()[0].set_markerfacecolor(color)
    ax.get_lines()[0].set_markeredgecolor(color)
    ax.get_lines()[0].set_markersize(3)
    ax.get_lines()[1].set_color("red")
    ax.set_title(f"Q-Q plot: {var_name}")
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_qq_plots.png", dpi=300)
plt.close(fig)
print("  ✅ eda_qq_plots.png")

# 8c. Homocedasticidad visual: boxplots ya generados en sección 4
#     Test de Levene para %INH ~ método
print("\n── Homocedasticidad (Levene) ──")
for var_name, var_col, df_src, grupo in [
    ("Inhibición (%) ~ método", "porcentaje_inhibicion", inh_vals, "metodo_extraccion"),
]:
    grupos = [g[var_col].dropna().values for _, g in df_src.groupby(grupo)]
    if all(len(g) >= 2 for g in grupos):
        stat, p = stats.levene(*grupos)
        print(f"  {var_name}: Levene={stat:.4f}, p={p:.6f} "
              f"{'❌ Heterocedástico' if p < 0.05 else '✅ Homocedástico'}")


# ═════════════════════════════════════════════════════════════════════
# 9. INTERACCIONES PRELIMINARES
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  9. INTERACCIONES PRELIMINARES")
print("=" * 65)

# 9a. Interacción método × concentración (solo MAC)
fig, ax = plt.subplots(figsize=(8, 5))
interaction_data = inh_vals.groupby(
    ["metodo_extraccion", "concentracion_mg_ml"]
)["porcentaje_inhibicion"].agg(["mean", "sem"]).reset_index()

for metodo in ["maceracion", "soxhlet", "ultrasonido"]:
    sub = interaction_data[interaction_data["metodo_extraccion"] == metodo]
    ax.errorbar(sub["concentracion_mg_ml"], sub["mean"], yerr=sub["sem"],
                label=metodo, color=COLOR_MET[metodo],
                marker="o", capsize=4, linewidth=2)
ax.set_xlabel("Concentración (mg/mL)")
ax.set_ylabel("Inhibición media (%)")
ax.set_title("Interacción: método × concentración")
ax.legend(title="Método")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_interaccion_metodo_conc.png", dpi=300)
plt.close(fig)
print("  ✅ eda_interaccion_metodo_conc.png")

# 9b. Perfil de inhibición por aislado (solo MAC, todas las concentraciones)
fig, ax = plt.subplots(figsize=(12, 6))
mac_perfil = inh_vals[inh_vals["metodo_extraccion"] == "maceracion"].copy()
perfil = mac_perfil.groupby(["aislado_id", "concentracion_mg_ml"])[
    "porcentaje_inhibicion"
].mean().reset_index()

# Seleccionar ~8 aislados representativos
aislados_destacados = ["HC3", "HC5", "HC17", "H11G", "FU1", "FUSARIUM JULIAN H20",
                        "H4B", "H8N"]
for aislado in aislados_destacados:
    sub = perfil[perfil["aislado_id"] == aislado]
    ax.plot(sub["concentracion_mg_ml"], sub["porcentaje_inhibicion"],
            marker="o", label=aislado, linewidth=2)
ax.set_xlabel("Concentración (mg/mL)")
ax.set_ylabel("Inhibición media (%)")
ax.set_title("Perfil de inhibición por aislado — Maceración")
ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
fig.subplots_adjust(right=0.8)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_perfil_aislados_mac.png", dpi=300)
plt.close(fig)
print("  ✅ eda_perfil_aislados_mac.png")


# ═════════════════════════════════════════════════════════════════════
# 10. CONIDIAS: ANÁLISIS ESPECÍFICO
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  10. ANÁLISIS DE CONIDIAS")
print("=" * 65)

# 10a. Boxplot: log₁₀ crudo por método y concentración
fig, ax = plt.subplots(figsize=(10, 5))
coni_trat = CONI[~CONI["es_control"] & CONI["conidias_log10"].notna()]
sns.boxplot(data=coni_trat, x="metodo_extraccion", y="conidias_log10",
            hue="concentracion_mg_ml", palette=COLOR_CONC, ax=ax)
ax.set_xlabel("Método de extracción")
ax.set_ylabel("log₁₀(conidias/mL)")
ax.set_title("Conidias — valores crudos en log₁₀")
ax.legend(title="Conc. (mg/mL)")
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_boxplot_conidias_log10.png", dpi=300)
plt.close(fig)
print("  ✅ eda_boxplot_conidias_log10.png")

# 10b. Boxplot: %INH de conidias sobre ESCALA CRUDA
fig, ax = plt.subplots(figsize=(10, 5))
coni_inh = CONI[~CONI["es_control"] & CONI["porcentaje_inhibicion"].notna()]
sns.boxplot(data=coni_inh, x="metodo_extraccion", y="porcentaje_inhibicion",
            hue="concentracion_mg_ml", palette=COLOR_CONC, ax=ax)
ax.set_xlabel("Método de extracción")
ax.set_ylabel("Reducción de conidias (%) — escala cruda")
ax.set_title("Inhibición de conidias sobre conteos crudos")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
ax.legend(title="Conc. (mg/mL)", bbox_to_anchor=(1.02, 1), loc="upper left")
fig.subplots_adjust(right=0.82)
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_boxplot_inhibicion_conidias_crudo.png", dpi=300)
plt.close(fig)
print("  ✅ eda_boxplot_inhibicion_conidias_crudo.png")

# 10c. Comparación %INH crudo vs log10
fig, ax = plt.subplots(figsize=(8, 6))
comp = CONI[CONI["porcentaje_inhibicion"].notna() & CONI["porcentaje_inhibicion_log10"].notna()]
ax.scatter(comp["porcentaje_inhibicion_log10"], comp["porcentaje_inhibicion"],
           alpha=0.3, c="#a23b72", edgecolors="none")
ax.plot([-100, 100], [-100, 100], "r--", alpha=0.5, label="Identidad")
ax.set_xlabel("%INH sobre log₁₀ (escala de la hoja)")
ax.set_ylabel("%INH sobre conteos crudos")
ax.set_title("Comparación: %INH de conidias en escala cruda vs log₁₀")
ax.legend()
fig.tight_layout()
fig.savefig(DIR_FIGURAS / "eda_comparacion_inh_conidias_crudo_vs_log10.png", dpi=300)
plt.close(fig)
print("  ✅ eda_comparacion_inh_conidias_crudo_vs_log10.png")

# 10d. Correlación crecimiento vs conidias (crudo)
fig, ax = plt.subplots(figsize=(8, 6))
merged = CREC[CREC["crecimiento_mm"].notna() & ~CREC["es_control"]].merge(
    CONI[CONI["conidias_log10"].notna() & ~CONI["es_control"]],
    on=["aislado_id", "metodo_extraccion", "concentracion_mg_ml", "replica_biologica"],
    suffixes=("_crec", "_con")
)
if len(merged) > 0:
    ax.scatter(merged["crecimiento_mm"], merged["conidias_log10"],
               alpha=0.4, c="#2e86ab", edgecolors="none")
    r, p = stats.pearsonr(merged["crecimiento_mm"], merged["conidias_log10"])
    ax.set_xlabel("Crecimiento (mm)")
    ax.set_ylabel("log₁₀(conidias/mL)")
    ax.set_title(f"Correlación crecimiento vs. conidias (r={r:.3f}, p={p:.4f})")
    fig.tight_layout()
    fig.savefig(DIR_FIGURAS / "eda_correlacion_crecimiento_conidias.png", dpi=300)
    plt.close(fig)
    print("  ✅ eda_correlacion_crecimiento_conidias.png")
else:
    plt.close(fig)
    print("  ⚠ Sin datos coincidentes para correlación crecimiento vs conidias")


# ═════════════════════════════════════════════════════════════════════
# 11. RESUMEN FINAL
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  11. RESUMEN — HALLAZGOS DEL EDA")
print("=" * 65)

hallazgos = []

# Hallazgo 1: Efecto piso
pct_completa = 100 * len(completa) / len(inh_vals)
hallazgos.append(f"• Inhibición completa en {len(completa)}/{len(inh_vals)} "
                 f"mediciones ({pct_completa:.1f}%) — efecto piso significativo")

# Hallazgo 2: Inhibición negativa
pct_negativa = 100 * len(negativa) / len(inh_vals)
hallazgos.append(f"• Inhibición negativa en {len(negativa)}/{len(inh_vals)} "
                 f"mediciones ({pct_negativa:.1f}%) — posible error de medición o variación natural")

# Hallazgo 3: Desbalance estructural
hallazgos.append("• Solo MACERACIÓN tiene diseño factorial completo (4 concentraciones)")
hallazgos.append("• SOXHLET y ULTRASONIDO tienen solo 2 concentraciones → sin dosis-respuesta")

# Hallazgo 4: Faltantes
for nombre, df_src in [("crecimiento", CREC), ("conidias", CONI)]:
    total = len(df_src)
    nulos = df_src["crecimiento_mm" if nombre == "crecimiento" else "conidias_log10"].isna().sum()
    hallazgos.append(f"• {nombre}: {nulos}/{total} filas sin dato ({100*nulos/total:.0f}%)")

# Hallazgo 5: Escala de conidias — crudo vs log10
if "porcentaje_inhibicion" in CONI.columns and "porcentaje_inhibicion_log10" in CONI.columns:
    comp_coni = CONI[CONI["porcentaje_inhibicion"].notna() & CONI["porcentaje_inhibicion_log10"].notna()]
    media_diff = (comp_coni["porcentaje_inhibicion"] - comp_coni["porcentaje_inhibicion_log10"]).mean()
    hallazgos.append(f"• Conidias: %INH en crudo difiere en {media_diff:.0f} puntos promedio del %INH en log₁₀ "
                     f"— la escala log₁₀ subestima drásticamente la reducción real")
    # Mostrar que SOX/ULT tienen poca reducción en crudo
    for metodo in ["soxhlet", "ultrasonido", "maceracion"]:
        sub = comp_coni[comp_coni["metodo_extraccion"] == metodo]
        media_crudo = sub["porcentaje_inhibicion"].mean()
        hallazgos.append(f"  • {metodo}: reducción media en crudo = {media_crudo:.1f}% "
                         f"(n={len(sub)} observaciones)")

# Hallazgo 6: Normalidad
for var_name, var_col, df_src in [
    ("Crecimiento", "crecimiento_mm", CREC),
    ("Inhibición", "porcentaje_inhibicion", inh_vals),
    ("Conidias (log₁₀)", "conidias_log10", CONI),
]:
    vals = df_src[var_col].dropna()
    if len(vals) >= 3:
        if len(vals) > 5000:
            vals = vals.sample(5000, random_state=SEMILLA_ALEATORIA)
        _, p = stats.shapiro(vals)
        hallazgos.append(f"• {var_name}: {'NO normal' if p < 0.05 else 'Normal'} "
                         f"(Shapiro-Wilk p={p:.4f})")

for h in hallazgos:
    print(f"  {h}")

# ─── Guardar reporte ────────────────────────────────────────────────
with open(DIR_REPORTES / "02_eda.md", "w", encoding="utf-8") as f:
    f.write("# Reporte de Análisis Exploratorio\n\n")
    f.write(f"**Fecha:** 2026-07-28\n\n")
    f.write(f"## Resumen de datos\n\n")
    f.write(f"- Crecimiento micelial: {len(CREC)} observaciones, {CREC['aislado_id'].nunique()} aislados\n")
    f.write(f"- Conidias: {len(CONI)} observaciones, {CONI['aislado_id'].nunique()} aislados\n")
    f.write(f"- Rendimiento: {len(rendimiento)} observaciones\n\n")
    f.write(f"## Hallazgos principales\n\n")
    for h in hallazgos:
        f.write(f"{h}\n\n")
    f.write(f"## Figuras generadas\n\n")
    import os
    for fig_file in sorted(os.listdir(DIR_FIGURAS)):
        if fig_file.endswith(".png"):
            f.write(f"- `{fig_file}`\n")

print(f"\n  ✅ Reporte guardado: {DIR_REPORTES / '02_eda.md'}")
print(f"\n{'='*65}")
print("  EDA COMPLETO — todas las figuras en dca/resultados/figuras/")
print(f"{'='*65}")
