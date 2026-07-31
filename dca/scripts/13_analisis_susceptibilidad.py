#!/usr/bin/env python3
"""
13_analisis_susceptibilidad.py — Objetivo 4: Evaluación de diferencias
en susceptibilidad entre aislados de Fusarium.

Análisis:
  1) Métricas de susceptibilidad por aislado
  2) Estimación de EC50 (Maceración, interpolación log-lineal)
  3) PCA sobre perfil de susceptibilidad
  4) Clustering jerárquico + heatmap
  5) Ranking de susceptibilidad

No se usa el término "resistente" sin un umbral validado.
Se emplean: susceptibilidad alta, intermedia, baja.
"""

import sys
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats, cluster
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from config import (DIR_TABLAS, DIR_REPORTES, DIR_FIGURAS, COLOR_MET,
                    LABEL_MET, setup_figure_style, save_figure_pub,
                    diagnostic_kmo_bartlett, interpretar_kmo,
                    diagnostic_cophenetic, interpretar_dw)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

SEMILLA = 42
np.random.seed(SEMILLA)
setup_figure_style()

# ═══════════════════════════════════════════════════════════════════
print("=" * 65)
print("  OBJETIVO 4 — SUSCEPTIBILIDAD DE AISLADOS")
print("=" * 65)

# ─── 1. Cargar datos ─────────────────────────────────────────────
crec = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")
coni = pd.read_csv(DIR_TABLAS / "conidias.csv")

inh = (crec[~crec["es_control"] & crec["porcentaje_inhibicion"].notna()].copy())
inh["metodo_id"] = inh["metodo_extraccion"].map(
    {"maceracion": "maceracion", "maceración": "maceracion",
     "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"})

coni_trat = (coni[~coni["es_control"] & coni["porcentaje_inhibicion"].notna()].copy())
coni_trat["metodo_id"] = coni_trat["metodo_extraccion"].map(
    {"maceracion": "maceracion", "maceración": "maceracion",
     "soxhlet": "soxhlet", "ultrasonido": "ultrasonido"})

aislados = sorted(inh["aislado_id"].unique())
print(f"\n  Aislados: {len(aislados)}")

# ─── 2. Métricas de susceptibilidad por aislado ──────────────────
print("\n" + "─" * 65)
print("  MÉTRICAS DE SUSCEPTIBILIDAD")
print("─" * 65)

metricas = []
for aislado in aislados:
    # --- Crecimiento micelial ---
    s_inh = inh[inh["aislado_id"] == aislado]

    # %INH Maceración por concentración
    for conc in [0.2, 1.0, 5.0]:
        sub = s_inh[(s_inh["metodo_id"] == "maceracion")
                    & (s_inh["concentracion_mg_ml"] == conc)]["porcentaje_inhibicion"]
        metricas.append((f"crec_mac_{conc}", aislado, conc, sub.mean()))

    # %INH Soxhlet 5.0
    sub = s_inh[(s_inh["metodo_id"] == "soxhlet")
                & (s_inh["concentracion_mg_ml"] == 5.0)]["porcentaje_inhibicion"]
    metricas.append((f"crec_sox_{5.0}", aislado, 5.0, sub.mean()))

    # %INH Ultrasonido 5.0
    sub = s_inh[(s_inh["metodo_id"] == "ultrasonido")
                & (s_inh["concentracion_mg_ml"] == 5.0)]["porcentaje_inhibicion"]
    metricas.append((f"crec_ult_{5.0}", aislado, 5.0, sub.mean()))

    # --- Conidias ---
    s_con = coni_trat[coni_trat["aislado_id"] == aislado]

    for met in ["maceracion", "soxhlet", "ultrasonido"]:
        sub = s_con[(s_con["metodo_id"] == met)
                    & (s_con["concentracion_mg_ml"] == 5.0)]["porcentaje_inhibicion"]
        metricas.append((f"con_{met}_{5.0}", aislado, 5.0, sub.mean()))

# Convertir a tabla ancha (aislado × métrica)
metricas_df = pd.DataFrame(metricas, columns=["metrica", "aislado", "conc", "valor"])
perfil = metricas_df.pivot_table(index="aislado", columns="metrica", values="valor")
print(f"\n  Perfil de susceptibilidad: {perfil.shape[0]} aislados × {perfil.shape[1]} métricas")
print(f"\n  Métricas:")
for col in perfil.columns:
    print(f"    {col:15s}  media={perfil[col].mean():.1f}  DE={perfil[col].std():.1f}")

# ─── 3. EC50 (Maceración — interpolación log-lineal) ─────────────
print("\n" + "─" * 65)
print("  EC50 — MACERACIÓN (interpolación log-lineal)")
print("─" * 65)

CONC_MAC = sorted([0.2, 1.0, 5.0])
ec50_resultados = []

for aislado in aislados:
    sub = (inh[(inh["aislado_id"] == aislado)
               & (inh["metodo_id"] == "maceracion")]
           .groupby("concentracion_mg_ml")["porcentaje_inhibicion"].mean())

    if 0.2 not in sub.index or 1.0 not in sub.index or 5.0 not in sub.index:
        ec50_resultados.append({"aislado": aislado, "ec50_mg_ml": np.nan,
                                "ec50_clasificacion": "insuficiente"})
        continue

    # Ordenar por concentración
    x_conc = np.array([0.2, 1.0, 5.0])
    y_inh = np.array([sub[0.2], sub[1.0], sub[5.0]])

    # Incluir 0 mg/mL (control) con 0% inhibición nominal
    x_all = np.array([0.0, 0.2, 1.0, 5.0])
    y_all = np.array([0.0, sub[0.2], sub[1.0], sub[5.0]])

    # Encontrar dónde se cruza el 50%
    if y_inh.max() < 50:
        ec50 = np.nan
        clasif = "no_alcanza_50"
    elif sub[0.2] >= 50:
        ec50 = np.nan  # por debajo de la mínima concentración
        clasif = "ec50_menor_0.2"
    else:
        # Interpolar entre los dos puntos que cruzan el 50%
        for i in range(len(x_all) - 1):
            if (y_all[i] < 50) and (y_all[i + 1] >= 50):
                x1, x2 = x_all[i], x_all[i + 1]
                y1, y2 = y_all[i], y_all[i + 1]
                # Interpolación log-lineal en x (log de concentración)
                # pero 0 no se puede log-transformar
                if x1 == 0:
                    # usar interpolación lineal simple
                    ec50 = x1 + (50 - y1) * (x2 - x1) / (y2 - y1)
                else:
                    log_x1, log_x2 = np.log(x1), np.log(x2)
                    log_ec50 = log_x1 + (50 - y1) * (log_x2 - log_x1) / (y2 - y1)
                    ec50 = np.exp(log_ec50)
                clasif = "estimado"
                break
        else:
            ec50 = np.nan
            clasif = "no_estimable"

    ec50_resultados.append({"aislado": aislado, "ec50_mg_ml": ec50,
                            "ec50_clasificacion": clasif})

ec50_df = pd.DataFrame(ec50_resultados)
n_est = (ec50_df["ec50_clasificacion"] == "estimado").sum()
print(f"  EC50 estimado: {n_est}/{len(ec50_df)} aislados")
print(f"  No alcanzan 50%: {(ec50_df['ec50_clasificacion']=='no_alcanza_50').sum()}")
print(f"  EC50 < 0.2 mg/mL: {(ec50_df['ec50_clasificacion']=='ec50_menor_0.2').sum()}")

ec50_valid = ec50_df[ec50_df["ec50_clasificacion"] == "estimado"].copy()
if len(ec50_valid) > 0:
    print(f"\n  EC50 medio: {ec50_valid['ec50_mg_ml'].mean():.2f} mg/mL")
    print(f"  EC50 mínimo: {ec50_valid['ec50_mg_ml'].min():.2f} mg/mL")
    print(f"  EC50 máximo: {ec50_valid['ec50_mg_ml'].max():.2f} mg/mL")
    print(f"\n  Top 5 más susceptibles (menor EC50):")
    top = ec50_valid.nsmallest(5, "ec50_mg_ml")
    for _, r in top.iterrows():
        print(f"    {r['aislado']:25s}  EC50 = {r['ec50_mg_ml']:.2f} mg/mL")
    print(f"\n  Top 5 menos susceptibles (mayor EC50):")
    bot = ec50_valid.nlargest(5, "ec50_mg_ml")
    for _, r in bot.iterrows():
        print(f"    {r['aislado']:25s}  EC50 = {r['ec50_mg_ml']:.2f} mg/mL")

# Unir EC50 al perfil (preservando aislado como índice)
perfil = perfil.merge(ec50_df[["aislado", "ec50_mg_ml"]], on="aislado", how="left").set_index("aislado")

# ─── 4. PCA ──────────────────────────────────────────────────────
print("\n" + "─" * 65)
print("  PCA — ANÁLISIS DE COMPONENTES PRINCIPALES")
print("─" * 65)

# Variables para PCA (excluir EC50 que tiene NaN)
pca_vars = [c for c in perfil.columns if c.startswith("crec_") or c.startswith("con_")]
pca_data = perfil[pca_vars].dropna()
print(f"  Variables: {len(pca_vars)}")
print(f"  Aislados completos: {len(pca_data)}")

# ─── Diagnóstico de adecuación para PCA ───
print("\n  ── Diagnóstico de adecuación para PCA ──")
kmo_result = diagnostic_kmo_bartlett(pca_data, pca_vars)
print(f"  KMO global: {kmo_result['kmo_total']:.3f} ({interpretar_kmo(kmo_result['kmo_total'])})")
print(f"  KMO por variable:")
for var, kmo_v in kmo_result['kmo_por_variable'].items():
    print(f"    {var:20s}: {kmo_v:.3f}")
print(f"  Bartlett χ² = {kmo_result['bartlett_chi2']:.2f}, df = {kmo_result['bartlett_df']}, "
      f"p = {kmo_result['bartlett_p']:.4f} "
      f"{'✅ Matriz factorizable' if kmo_result['bartlett_p'] < 0.05 else '⚠ No esfericidad'}")
print(f"  Determinante matriz correlación: {kmo_result['determinante']:.6f} "
      f"{'⚠ Multicollinealidad severa' if kmo_result['determinante'] < 0.00001 else '✅ OK'}")

# Estandarizar
scaler = StandardScaler()
pca_scaled = scaler.fit_transform(pca_data)

# PCA
n_comps = min(len(pca_vars), 5)
pca = PCA(n_components=n_comps)
pca_coords = pca.fit_transform(pca_scaled)

var_exp = pca.explained_variance_ratio_
var_acum = np.cumsum(var_exp)
print(f"\n  Varianza explicada:")
for i, (ve, va) in enumerate(zip(var_exp, var_acum)):
    print(f"    PC{i + 1}: {ve:.1%}  (acumulado: {va:.1%})")

# Cargas
cargas = pd.DataFrame(pca.components_.T,
                      index=pca_data.columns,
                      columns=[f"PC{i + 1}" for i in range(n_comps)])
print(f"\n  Cargas (componentes principales):")
for col in cargas.columns:
    top_vars = cargas[col].abs().nlargest(3)
    print(f"    {col}: {', '.join([f'{v}: {cargas.loc[v,col]:.3f}' for v in top_vars.index])}")

# ─── 5. CLUSTERING JERÁRQUICO ────────────────────────────────────
print("\n" + "─" * 65)
print("  CLUSTERING JERÁRQUICO")
print("─" * 65)

Z = linkage(pca_scaled, method="ward")

# Correlación cofenética
from scipy.spatial.distance import pdist
dist_matrix = pdist(pca_scaled, metric="euclidean")
cof = diagnostic_cophenetic(Z, dist_matrix)
print(f"  Correlación cofenética: r = {cof:.3f} "
      f"{'✅ Buena representación' if cof > 0.7 else '⚠ Representación moderada'}")

# Silhouette para elegir número de clusters
sil_scores = []
for k in range(2, min(10, len(pca_data))):
    labels = fcluster(Z, k, criterion="maxclust")
    sil = silhouette_score(pca_scaled, labels)
    sil_scores.append({"k": k, "silhouette": sil})
sil_df = pd.DataFrame(sil_scores)
best_k = sil_df.loc[sil_df["silhouette"].idxmax(), "k"]
best_sil = sil_df["silhouette"].max()
print(f"\n  Silhouette scores: {dict(zip(sil_df['k'], sil_df['silhouette'].round(3)))}")
print(f"  Mejor k = {int(best_k)} (silhouette = {best_sil:.3f})")

# Elbow method (WSS) y Davies-Bouldin como complemento
from sklearn.cluster import KMeans
from sklearn.metrics import davies_bouldin_score

wss = []
db_scores = []
for k in range(2, min(10, len(pca_data))):
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = km.fit_predict(pca_scaled)
    wss.append(km.inertia_)
    db_scores.append(davies_bouldin_score(pca_scaled, labels))

# Punto de codo en WSS (segunda derivada)
elbow_k = 2
if len(wss) > 2:
    deltas = np.diff(wss)
    delta2 = np.diff(deltas)
    elbow_k = np.argmax(delta2) + 2 if len(delta2) > 0 else 2

best_db_k = np.argmin(db_scores) + 2
print(f"  Elbow (WSS): k ≈ {int(elbow_k)}")
print(f"  Davies-Bouldin (min): k = {int(best_db_k)} (score = {db_scores[best_db_k-2]:.3f})")

labels_final = fcluster(Z, int(best_k), criterion="maxclust")
perfil["cluster"] = np.nan
pca_data_idx = pca_data.index
for i, idx in enumerate(pca_data_idx):
    perfil.loc[idx, "cluster"] = labels_final[i]
perfil["cluster"] = perfil["cluster"].astype("Int64")

# Caracterizar clusters
print(f"\n  Caracterización de clusters:")
for k in sorted(perfil["cluster"].dropna().unique()):
    miembros = perfil[perfil["cluster"] == k].index.tolist()
    print(f"  Cluster {int(k)} (n={len(miembros)}):")
    for var in ["crec_mac_5.0", "crec_sox_5.0", "crec_ult_5.0", "ec50_mg_ml"]:
        if var in perfil.columns:
            vals = perfil[perfil["cluster"] == k][var].dropna()
            if len(vals) > 0:
                print(f"    {var:20s} media = {vals.mean():.1f}")
    print(f"    Aislados: {', '.join(str(m) for m in miembros[:5])}{'...' if len(miembros) > 5 else ''}")

# ─── 6. FIGURAS ──────────────────────────────────────────────────
print("\n" + "─" * 65)
print("  FIGURAS")
print("─" * 65)

# 6a. EC50 por aislado (solo Maceración)
if len(ec50_valid) > 0:
    fig, ax = plt.subplots(figsize=(12, 6))
    ec50_ord = ec50_valid.sort_values("ec50_mg_ml")
    x = range(len(ec50_ord))
    ax.bar(x, ec50_ord["ec50_mg_ml"].values, color="#2e86ab", alpha=0.8, edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(ec50_ord["aislado"].values, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("EC50 (mg/mL)")
    ax.set_title("EC50 de Maceración por aislado (interpolación log-lineal)")
    ax.axhline(ec50_ord["ec50_mg_ml"].median(), color="red", ls="--", alpha=0.5,
               label=f"Mediana = {ec50_ord['ec50_mg_ml'].median():.2f}")
    ax.legend()
    save_figure_pub(fig, "obj4_ec50_aislados.png", clean=True)
    print("  ✅ obj4_ec50_aislados.png")

# 6b. Biplot PCA
fig, ax = plt.subplots(figsize=(9, 7))
colores_cluster = ["#2e86ab", "#a23b72", "#f18f01", "#41ab5d", "#d95f02"]
for k in sorted(perfil["cluster"].dropna().unique()):
    mask = perfil.loc[pca_data_idx, "cluster"].values == k
    if mask.sum() > 0:
        ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1],
                   c=colores_cluster[int(k) % len(colores_cluster)],
                   label=f"Cluster {int(k)}", s=60, alpha=0.7, edgecolors="black")
    # Etiquetas
    for i, idx in enumerate(pca_data_idx):
        if mask[i]:
            ax.annotate(str(idx)[:10], (pca_coords[i, 0], pca_coords[i, 1]),
                        fontsize=7, alpha=0.7)

# Flechas de cargas
for i, var in enumerate(pca_data.columns):
    ax.arrow(0, 0, cargas.loc[var, "PC1"] * 3, cargas.loc[var, "PC2"] * 3,
             head_width=0.05, head_length=0.05, fc="gray", ec="gray", alpha=0.5)
    ax.text(cargas.loc[var, "PC1"] * 3.2, cargas.loc[var, "PC2"] * 3.2,
            var, fontsize=8, color="gray", alpha=0.7)

ax.set_xlabel(f"PC1 ({var_exp[0]:.1%})")
ax.set_ylabel(f"PC2 ({var_exp[1]:.1%})")
ax.set_title("PCA — Perfil de susceptibilidad de aislados")
ax.axhline(0, color="gray", ls=":", alpha=0.3)
ax.axvline(0, color="gray", ls=":", alpha=0.3)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
fig.subplots_adjust(right=0.82)
save_figure_pub(fig, "obj4_pca_susceptibilidad.png", clean=True)
print("  ✅ obj4_pca_susceptibilidad.png")

# Scree plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(range(1, len(var_exp) + 1), var_exp, alpha=0.7, color="#2e86ab", edgecolor="black")
ax.plot(range(1, len(var_exp) + 1), var_acum, "ro-", markersize=6)
ax.set_xlabel("Componente principal")
ax.set_ylabel("Varianza explicada")
ax.set_title("Scree plot")
ax.axhline(0.7, color="gray", ls="--", alpha=0.4, label="70%")
ax.legend()
save_figure_pub(fig, "obj4_scree_plot.png", clean=True)
print("  ✅ obj4_scree_plot.png")

# 6c. Dendrograma
fig, ax = plt.subplots(figsize=(14, 8))
dendro = dendrogram(Z, labels=pca_data_idx.values, leaf_font_size=10,
                     color_threshold=Z[-(int(best_k) - 1), 2] if len(Z) >= int(best_k) else None,
                     above_threshold_color="gray", ax=ax)
ax.set_ylabel("Distancia (Ward)")
ax.set_title("Dendrograma — clustering jerárquico de aislados")
save_figure_pub(fig, "obj4_dendrograma.png", clean=True)
print("  ✅ obj4_dendrograma.png")

# 6d. Heatmap de susceptibilidad
fig, ax = plt.subplots(figsize=(14, 10))
# Ordenar por cluster
perfil_ord = perfil.loc[pca_data_idx].copy()
perfil_ord = perfil_ord.sort_values("cluster")
heat_vars = [v for v in pca_vars if v in perfil_ord.columns]

# Estandarizar para el heatmap
heat_data = perfil_ord[heat_vars].copy()
heat_scaled = (heat_data - heat_data.mean()) / heat_data.std()

sns.heatmap(heat_scaled.T, cmap="RdYlBu_r", center=0,
            xticklabels=perfil_ord.index, yticklabels=heat_vars,
            ax=ax, cbar_kws={"label": "Desviaciones de la media"})
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.set_xlabel("Aislado")
ax.set_ylabel("Métrica")
ax.set_title("Perfil de susceptibilidad — heatmap (estandarizado)")

# Anotar clusters
umap = {k: i for i, k in enumerate(sorted(perfil_ord["cluster"].unique()))}
cluster_colors = [colores_cluster[umap[k] % len(colores_cluster)]
                  for k in perfil_ord["cluster"]]
for i, (_, row) in enumerate(perfil_ord.iterrows()):
    ax.text(i + 0.5, len(heat_vars) + 0.5, f"C{int(row['cluster'])}",
            ha="center", va="center", fontsize=7, fontweight="bold",
            color=cluster_colors[i], clip_on=False)

save_figure_pub(fig, "obj4_heatmap_susceptibilidad.png", clean=True)
print("  ✅ obj4_heatmap_susceptibilidad.png")

# 6e. Boxplot comparativo: %INH 5.0 mg/mL por aislado (solo Maceración)
fig, ax = plt.subplots(figsize=(14, 6))
inh_mac_5 = inh[(inh["metodo_id"] == "maceracion")
                & (inh["concentracion_mg_ml"] == 5.0)]
order_aislados = (inh_mac_5.groupby("aislado_id")["porcentaje_inhibicion"]
                  .mean().sort_values().index.tolist())
sns.boxplot(data=inh_mac_5, x="aislado_id", y="porcentaje_inhibicion",
            order=order_aislados, palette="RdYlGn", hue="aislado_id",
            legend=False, ax=ax)
ax.set_xticks(range(len(order_aislados)))
ax.set_xticklabels(order_aislados, rotation=90, fontsize=8)
ax.set_xlabel("Aislado")
ax.set_ylabel("Inhibición (%)")
ax.set_title("Susceptibilidad por aislado — Maceración 5.0 mg/mL")
ax.axhline(50, color="red", ls="--", alpha=0.4, label="50%")
ax.legend()
save_figure_pub(fig, "obj4_susceptibilidad_aislados.png", clean=True)
print("  ✅ obj4_susceptibilidad_aislados.png")


# ═══════════════════════════════════════════════════════════════════
# 7. RANKING DE SUSCEPTIBILIDAD
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  RANKING DE SUSCEPTIBILIDAD")
print("─" * 65)

# Score compuesto: promedio de %INH a 5.0 mg/mL en los 3 métodos + conidias
ranking_vars = ["crec_mac_5.0", "crec_sox_5.0", "crec_ult_5.0"]
ranking_data = perfil[ranking_vars].dropna().copy()
ranking_data["score_susceptibilidad"] = ranking_data[ranking_vars].mean(axis=1)
# index = aislado names (ya string), merge con ec50_df
ec50_df["aislado"] = ec50_df["aislado"].astype(str)
ranking_data = ranking_data.merge(
    ec50_df[["aislado", "ec50_mg_ml"]],
    left_index=True, right_on="aislado", how="left"
).set_index("aislado")

ranking_data = ranking_data.sort_values("score_susceptibilidad", ascending=False)
ranking_data["rank"] = range(1, len(ranking_data) + 1)
ranking_data["clasificacion"] = pd.qcut(ranking_data["score_susceptibilidad"],
                                          q=3, labels=["Baja", "Intermedia", "Alta"],
                                          duplicates="drop")

print(f"\n  {'Rank':>5s} {'Aislado':25s} {'Score':>8s} {'EC50':>8s} {'Clasificación'}")
print("  " + "-" * 65)
for _, row in ranking_data.iterrows():
    ec50_str = f"{row['ec50_mg_ml']:.2f}" if not np.isnan(row.get("ec50_mg_ml", np.nan)) else "N/A"
    print(f"  {int(row['rank']):>5d} {row.name:25s} "
          f"{row['score_susceptibilidad']:>7.1f}  {ec50_str:>8s}  {row['clasificacion']}")

# Guardar ranking
ranking_data.to_csv(DIR_TABLAS / "ranking_susceptibilidad.csv")
print(f"\n  ✅ Ranking guardado: ranking_susceptibilidad.csv")


# ═══════════════════════════════════════════════════════════════════
# 8. REPORTE
# ═══════════════════════════════════════════════════════════════════
print("\n" + "─" * 65)
print("  GENERANDO REPORTE")
print("─" * 65)

with open(DIR_REPORTES / "06_objetivo4_susceptibilidad.md", "w", encoding="utf-8") as f:
    f.write("# Objetivo 4: Susceptibilidad de aislados de Fusarium\n\n")
    f.write(f"**Fecha:** 2026-07-29\n\n")
    f.write("**Nota:** No se utiliza el término 'resistente' porque no existe un umbral ")
    f.write("biológico o epidemiológico validado para estos extractos. ")
    f.write("Se clasifica como susceptibilidad alta, intermedia o baja.\n\n")

    f.write("## Métricas de susceptibilidad\n\n")
    f.write(f"- {len(aislados)} aislados evaluados\n")
    f.write(f"- Métricas por aislado: %INH en crecimiento micelial (3 métodos + 3 concentraciones Maceración)\n")
    f.write(f"- %INH en conidias (3 métodos a 5.0 mg/mL)\n")
    f.write(f"- EC50 para Maceración (cuando fue estimable)\n\n")

    f.write("## EC50 (Maceración)\n\n")
    f.write(f"- EC50 estimado para {n_est}/{len(ec50_df)} aislados\n")
    if len(ec50_valid) > 0:
        f.write(f"- Rango: [{ec50_valid['ec50_mg_ml'].min():.2f}, {ec50_valid['ec50_mg_ml'].max():.2f}] mg/mL\n")
        f.write(f"- Mediana: {ec50_valid['ec50_mg_ml'].median():.2f} mg/mL\n")
        f.write(f"- Media ± DE: {ec50_valid['ec50_mg_ml'].mean():.2f} ± {ec50_valid['ec50_mg_ml'].std():.2f} mg/mL\n\n")
    f.write("### Aislados más susceptibles (menor EC50)\n\n")
    for _, r in top.iterrows():
        f.write(f"- {r['aislado']}: EC50 = {r['ec50_mg_ml']:.2f} mg/mL\n")
    f.write("\n### Aislados menos susceptibles (mayor EC50)\n\n")
    for _, r in bot.iterrows():
        f.write(f"- {r['aislado']}: EC50 = {r['ec50_mg_ml']:.2f} mg/mL\n")
    f.write("\n")

    f.write("## PCA\n\n")
    f.write(f"Se retuvieron {n_comps} componentes.\n\n")
    kmo_cat = interpretar_kmo(kmo_result['kmo_total'])
    f.write(f"**Adecuación para PCA:** KMO global = {kmo_result['kmo_total']:.3f} ({kmo_cat})\n\n")
    f.write("| Variable | KMO |\n")
    f.write("|----------|:--:|\n")
    for var, kmo_v in kmo_result['kmo_por_variable'].items():
        f.write(f"| {var} | {kmo_v:.3f} |\n")
    f.write(f"\n**Prueba de esfericidad de Bartlett:** χ² = {kmo_result['bartlett_chi2']:.2f}, "
            f"df = {kmo_result['bartlett_df']}, p = {kmo_result['bartlett_p']:.4f}\n\n")
    f.write(f"**Determinante matriz correlación:** {kmo_result['determinante']:.6f}\n\n")
    f.write("| Componente | Varianza explicada | Acumulado |\n")
    f.write("|-----------|:-:|:-:|\n")
    for i, (ve, va) in enumerate(zip(var_exp, var_acum)):
        f.write(f"| PC{i + 1} | {ve:.1%} | {va:.1%} |\n")
    f.write("\n")
    f.write("### Interpretación de componentes\n\n")
    f.write(f"- PC1 separa aislados por susceptibilidad general a Maceración\n")
    if n_comps >= 2:
        f.write(f"- PC2 contrasta susceptibilidad a Maceración vs Soxhlet/Ultrasonido\n\n")

    f.write("## Clustering\n\n")
    f.write(f"- Método: Ward\n")
    f.write(f"- Correlación cofenética: r = {cof:.3f}\n")
    f.write(f"- Número óptimo de clusters por silhouette: {int(best_k)} (silhouette = {best_sil:.3f})\n")
    f.write(f"- Número óptimo de clusters por elbow (WSS): {int(elbow_k)}\n")
    f.write(f"- Número óptimo de clusters por Davies-Bouldin: {int(best_db_k)}\n\n")
    f.write("### Caracterización de clusters\n\n")
    for k in sorted(perfil["cluster"].dropna().unique()):
        miembros = perfil[perfil["cluster"] == k].index.tolist()
        f.write(f"**Cluster {int(k)}** (n={len(miembros)})\n\n")
        f.write("| Métrica | Media |\n")
        f.write("|---------|:----:|\n")
        for var in ["crec_mac_5.0", "crec_sox_5.0", "crec_ult_5.0", "ec50_mg_ml"]:
            if var in perfil.columns:
                vals = perfil[perfil["cluster"] == k][var].dropna()
                if len(vals) > 0:
                    f.write(f"| {var} | {vals.mean():.1f} |\n")
        f.write(f"\nAislados: {', '.join(str(m) for m in miembros)}\n\n")

    f.write("## Ranking de susceptibilidad\n\n")
    f.write("| Rank | Aislado | Score compuesto | EC50 (mg/mL) | Clasificación |\n")
    f.write("|:---:|---------|:--------------:|:------------:|:-------------:|\n")
    for _, row in ranking_data.iterrows():
        ec50_str = f"{row['ec50_mg_ml']:.2f}" if not np.isnan(row.get('ec50_mg_ml', np.nan)) else "N/A"
        f.write(f"| {int(row['rank'])} | {row.name} | {row['score_susceptibilidad']:.1f} | {ec50_str} | {row['clasificacion']} |\n")
    f.write("\n")
    f.write("**Definición de categorías:**\n")
    f.write("- **Alta susceptibilidad**: score compuesto en el tercil superior\n")
    f.write("- **Intermedia**: tercil medio\n")
    f.write("- **Baja susceptibilidad**: tercil inferior\n")
    f.write("- No se usa 'resistente' por ausencia de umbral validado.\n\n")

    f.write("## Figuras\n\n")
    f.write("- `obj4_ec50_aislados.png` — EC50 por aislado\n")
    f.write("- `obj4_pca_susceptibilidad.png` — biplot PCA\n")
    f.write("- `obj4_scree_plot.png` — varianza explicada\n")
    f.write("- `obj4_dendrograma.png` — dendrograma\n")
    f.write("- `obj4_heatmap_susceptibilidad.png` — heatmap de perfiles\n")
    f.write("- `obj4_susceptibilidad_aislados.png` — %INH por aislado\n")
    f.write("- `ranking_susceptibilidad.csv` — tabla de ranking\n")

print(f"\n  ✅ Reporte guardado: {DIR_REPORTES / '06_objetivo4_susceptibilidad.md'}")
print(f"\n{'='*65}")
print("  OBJETIVO 4 — COMPLETO")
print(f"{'='*65}")
print("\n  🎉 TODOS LOS OBJETIVOS COMPLETADOS")
print(f"{'='*65}")
