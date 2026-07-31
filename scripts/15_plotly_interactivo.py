"""
15_plotly_interactivo.py — Figuras interactivas con Plotly.

Genera HTML autónomos con gráficos interactivos (hover, zoom, toggle)
para exploración visual de los resultados.

Requiere: plotly >= 5.0
Salida:  resultados/figuras_interactivas/*.html
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from config import (
    DIR_RESULTADOS, DIR_FIGURAS, DIR_TABLAS, COLOR_MET, LABEL_MET
)
import warnings
warnings.filterwarnings("ignore")

# ─── Plotly ────────────────────────────────────────────────────────
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

DIR_INTERACTIVO = DIR_RESULTADOS / "figuras_interactivas"
DIR_INTERACTIVO.mkdir(parents=True, exist_ok=True)

# ─── Paleta consistente ────────────────────────────────────────────
COLOR_MET_LIST = [COLOR_MET["maceracion"], COLOR_MET["soxhlet"], COLOR_MET["ultrasonido"]]
COLOR_CONC_PLOTLY = ["#b3b3b3", "#7ba0b4", "#4a7c9b", "#1a4d6b"]
ORDER_MET = ["maceracion", "soxhlet", "ultrasonido"]

# Template base
TEMPLATE = "plotly_white"


# ═══════════════════════════════════════════════════════════════════
# 1. CARGA DE DATOS
# ═══════════════════════════════════════════════════════════════════
print("Cargando datos...")
crec = pd.read_csv(DIR_TABLAS / "crecimiento_micelial.csv")
coni = pd.read_csv(DIR_TABLAS / "conidias.csv")
rend = pd.read_csv(DIR_TABLAS / "rendimiento_extraccion.csv")
ranking = pd.read_csv(DIR_TABLAS / "ranking_susceptibilidad.csv")

# Datos de tratamiento (sin controles)
inh = crec[~crec["es_control"] & crec["porcentaje_inhibicion"].notna()].copy()
inh["metodo_label"] = inh["metodo_extraccion"].map(LABEL_MET)
inh["conc_label"] = "Control"
inh.loc[inh["concentracion_mg_ml"] > 0, "conc_label"] = (
    inh.loc[inh["concentracion_mg_ml"] > 0, "concentracion_mg_ml"].astype(str) + " mg/mL"
)
inh_5 = inh[inh["concentracion_mg_ml"] == 5.0].copy()

# Conidias tratamiento
coni_trat = coni[~coni["es_control"] & coni["conidias_log10"].notna()].copy()
coni_trat["metodo_label"] = coni_trat["metodo_extraccion"].map(LABEL_MET)


# ═══════════════════════════════════════════════════════════════════
# 2. FUNCIONES AUXILIARES
# ═══════════════════════════════════════════════════════════════════
def save(fig, nombre, width=900, height=550):
    """Guardar figura como HTML autónomo."""
    path = DIR_INTERACTIVO / nombre
    fig.write_html(str(path), include_plotlyjs="cdn", 
                   full_html=True, auto_open=False)
    print(f"  ✅ {path.name} ({width}×{height})")
    return path


def hover_text_aislado(row, extra_cols=None):
    """Texto informativo para hover en figuras de aislados."""
    base = (
        f"<b>{row.get('aislado_id', '')}</b><br>"
        f"{row.get('metodo_label', row.get('metodo_extraccion', ''))}"
    )
    if "replica_biologica" in row.index and pd.notna(row["replica_biologica"]):
        base += f" | R{int(row['replica_biologica'])}"
    if "concentracion_mg_ml" in row.index:
        base += f"<br>{row['concentracion_mg_ml']} mg/mL"
    if "porcentaje_inhibicion" in row.index and pd.notna(row["porcentaje_inhibicion"]):
        base += f"<br>%INH: {row['porcentaje_inhibicion']:.1f}%"
    if extra_cols:
        for col, fmt in extra_cols:
            if col in row.index and pd.notna(row[col]):
                val = row[col]
                if fmt == ".2f":
                    base += f"<br>{col}: {val:.2f}"
                elif fmt == ".1f":
                    base += f"<br>{col}: {val:.1f}"
                elif fmt == "int":
                    base += f"<br>{col}: {int(val)}"
    return base


# ═══════════════════════════════════════════════════════════════════
# 3. FIGURAS
# ═══════════════════════════════════════════════════════════════════

# ─── 3a. %INH por método a 5.0 mg/mL (boxplot + strip) ──────────
print("\nGenerando figuras interactivas...")

fig = go.Figure()
for i, met in enumerate(ORDER_MET):
    sub = inh_5[inh_5["metodo_extraccion"] == met]
    hover = sub.apply(lambda r: hover_text_aislado(r, [("crecimiento_mm", ".1f")]), axis=1)
    
    # Box
    fig.add_trace(go.Box(
        y=sub["porcentaje_inhibicion"], name=LABEL_MET[met],
        marker_color=COLOR_MET[met], boxmean="sd",
        legendgroup=met, text=hover, hovertemplate="%{text}<br>%INH: %{y:.1f}%<extra></extra>",
        boxpoints=False, width=0.5,
    ))
    # Strip
    fig.add_trace(go.Scatter(
        y=sub["porcentaje_inhibicion"],
        x=[LABEL_MET[met]] * len(sub),
        mode="markers", marker=dict(color=COLOR_MET[met], size=6, opacity=0.4),
        showlegend=False, legendgroup=met,
        text=hover, hovertemplate="%{text}<extra></extra>",
    ))

fig.update_layout(
    title="Inhibición de crecimiento micelial a 5.0 mg/mL por método de extracción",
    yaxis_title="Inhibición (%)",
    xaxis_title="Método de extracción",
    template=TEMPLATE,
    width=700, height=500,
)
save(fig, "inh_metodo_5mg_ml.html")


# ─── 3b. Dosis-respuesta Maceración ───────────────────────────────
mac = inh[inh["metodo_extraccion"] == "maceracion"].copy()
# Boxplot por concentración
fig = go.Figure()
hover_mac = mac.apply(lambda r: hover_text_aislado(r), axis=1)

for i, conc in enumerate(sorted(mac["concentracion_mg_ml"].unique())):
    sub = mac[mac["concentracion_mg_ml"] == conc]
    label = f"{conc} mg/mL" if conc > 0 else "Control"
    fig.add_trace(go.Box(
        y=sub["porcentaje_inhibicion"], name=label,
        marker_color=COLOR_CONC_PLOTLY[i],
        boxmean="sd", boxpoints=False, width=0.5,
        text=sub.apply(lambda r: hover_text_aislado(r), axis=1),
        hovertemplate="%{text}<br>%INH: %{y:.1f}%<extra></extra>",
    ))

fig.update_layout(
    title="Maceración — efecto dosis-respuesta",
    yaxis_title="Inhibición (%)",
    xaxis_title="Concentración",
    template=TEMPLATE,
    width=700, height=500,
)
save(fig, "dosis_respuesta_maceracion.html")


# ─── 3c. Rendimiento de extracción ────────────────────────────────
fig = go.Figure()
for i, met in enumerate(ORDER_MET):
    sub = rend[rend["metodo_extraccion"] == met]
    label = LABEL_MET.get(met, met)
    fig.add_trace(go.Bar(
        name=label, x=[label],
        y=[sub["rendimiento_pct"].mean()],
        error_y=dict(type="data", array=[sub["rendimiento_pct"].std()], visible=True),
        marker_color=COLOR_MET[met],
        text=f"{sub['rendimiento_pct'].mean():.1f}%",
        textposition="outside",
        hovertemplate=f"<b>{label}</b><br>Rendimiento: %{{y:.1f}}%<br>"
                      f"n={len(sub)}<extra></extra>",
    ))

fig.update_layout(
    title="Rendimiento de extracción por método",
    yaxis_title="Rendimiento (%)",
    xaxis_title="",
    template=TEMPLATE,
    width=600, height=450,
    showlegend=False,
)
save(fig, "rendimiento_extraccion.html")


# ─── 3d. PCA biplot de susceptibilidad ────────────────────────────
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Reconstruir PCA desde el ranking
pca_vars = ["crec_mac_5.0", "crec_sox_5.0", "crec_ult_5.0"]
pca_data = ranking[pca_vars].dropna()
scaler = StandardScaler()
X_scaled = scaler.fit_transform(pca_data)
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)
expl_var = pca.explained_variance_ratio_

# Colorear por clasificación
clases = ranking.loc[pca_data.index, "clasificacion"]
paleta_cluster = {"Alta": "#2e86ab", "Intermedia": "#f18f01", "Baja": "#a23b72"}

fig = go.Figure()
for cls in ["Alta", "Intermedia", "Baja"]:
    mask = clases == cls
    if mask.sum() == 0:
        continue
    fig.add_trace(go.Scatter(
        x=X_pca[mask, 0], y=X_pca[mask, 1],
        mode="markers+text", name=cls,
        marker=dict(size=10, color=paleta_cluster[cls], line=dict(width=1, color="black")),
        text=ranking.loc[pca_data.index[mask], "aislado"],
        textposition="top center", textfont=dict(size=9),
        hovertemplate="<b>%{text}</b><br>"
                     f"PC1: %{{x:.2f}}<br>PC2: %{{y:.2f}}<br>"
                     f"Clasificación: {cls}<extra></extra>",
    ))

# Vectores de carga
for i, var in enumerate(pca_vars):
    fig.add_annotation(
        x=pca.components_[0, i] * 3,
        y=pca.components_[1, i] * 3,
        text=var, showarrow=True,
        arrowhead=2, ax=0, ay=0,
        font=dict(size=11),
    )

fig.update_layout(
    title=f"PCA — Susceptibilidad de aislados (PC1={expl_var[0]*100:.1f}%, PC2={expl_var[1]*100:.1f}%)",
    xaxis_title=f"PC1 ({expl_var[0]*100:.1f}%)",
    yaxis_title=f"PC2 ({expl_var[1]*100:.1f}%)",
    template=TEMPLATE,
    width=800, height=600,
)
save(fig, "pca_susceptibilidad.html")


# ─── 3e. EC50 por aislado ─────────────────────────────────────────
ec50_ord = ranking.dropna(subset=["ec50_mg_ml"]).sort_values("ec50_mg_ml")

fig = go.Figure()
fig.add_trace(go.Bar(
    y=ec50_ord["aislado"], x=ec50_ord["ec50_mg_ml"],
    orientation="h",
    marker=dict(
        color=ec50_ord["ec50_mg_ml"],
        colorscale="Blues", reversescale=True,
        line=dict(width=0.5, color="gray"),
    ),
    text=ec50_ord["ec50_mg_ml"].round(2),
    textposition="outside",
    hovertemplate="<b>%{y}</b><br>EC50: %{x:.2f} mg/mL<extra></extra>",
))

fig.update_layout(
    title="EC50 por aislado — Maceración",
    xaxis_title="EC50 (mg/mL)",
    yaxis_title="",
    template=TEMPLATE,
    width=700, height=max(400, 15 * len(ec50_ord)),
    margin=dict(l=120, r=40),
)
save(fig, "ec50_aislados.html")


# ─── 3f. Correlación crecimiento vs conidias ──────────────────────
# Merge a nivel de tratamiento 5.0 mg/mL
crec_5 = crec[~crec["es_control"] & (crec["concentracion_mg_ml"] == 5.0)][
    ["aislado_id", "metodo_extraccion", "replica_biologica", "porcentaje_inhibicion"]
].rename(columns={"porcentaje_inhibicion": "inh_crecimiento"})
coni_5 = coni[~coni["es_control"] & (coni["concentracion_mg_ml"] == 5.0)][
    ["aislado_id", "metodo_extraccion", "replica_biologica", "porcentaje_inhibicion_log10"]
].rename(columns={"porcentaje_inhibicion_log10": "inh_conidias"})
merged = pd.merge(crec_5, coni_5, on=["aislado_id", "metodo_extraccion", "replica_biologica"]).dropna()

fig = go.Figure()
for met in ORDER_MET:
    sub = merged[merged["metodo_extraccion"] == met]
    fig.add_trace(go.Scatter(
        x=sub["inh_crecimiento"], y=sub["inh_conidias"],
        mode="markers", name=LABEL_MET[met],
        marker=dict(color=COLOR_MET[met], size=7, opacity=0.6),
        text=sub.apply(lambda r: f"<b>{r['aislado_id']}</b><br>"
                                 f"Crec: {r['inh_crecimiento']:.1f}%<br>"
                                 f"Con: {r['inh_conidias']:.1f}%", axis=1),
        hovertemplate="%{text}<extra></extra>",
    ))

# Línea de tendencia global
from numpy.polynomial import polynomial as P
x_all = merged["inh_crecimiento"]
y_all = merged["inh_conidias"]
mask_t = ~np.isnan(x_all) & ~np.isnan(y_all)
if mask_t.sum() > 5:
    coefs = np.polyfit(x_all[mask_t], y_all[mask_t], 1)
    r_val = np.corrcoef(x_all[mask_t], y_all[mask_t])[0, 1]
    x_line = np.linspace(x_all[mask_t].min(), x_all[mask_t].max(), 100)
    fig.add_trace(go.Scatter(
        x=x_line, y=np.polyval(coefs, x_line),
        mode="lines", name=f"Tendencia (r={r_val:.3f})",
        line=dict(color="black", dash="dash", width=1.5),
    ))

fig.update_layout(
    title="Correlación: inhibición crecimiento vs inhibición conidias (5.0 mg/mL)",
    xaxis_title="Inhibición crecimiento (%)",
    yaxis_title="Inhibición conidias — log₁₀ (%)",
    template=TEMPLATE,
    width=700, height=550,
)
save(fig, "correlacion_crecimiento_conidias.html")


# ─── 3g. Heatmap de susceptibilidad (interactivo) ─────────────────
perfil = ranking.set_index("aislado")[pca_vars].dropna()

fig = go.Figure(data=go.Heatmap(
    z=perfil.values,
    x=perfil.columns,
    y=perfil.index,
    colorscale="YlOrRd",
    text=np.round(perfil.values, 1),
    texttemplate="%{text}",
    textfont=dict(size=9),
    hovertemplate="<b>%{y}</b><br>%{x}: %{z:.1f}%<extra></extra>",
))

fig.update_layout(
    title="Perfil de susceptibilidad por aislado (%INH a 5.0 mg/mL)",
    xaxis_title="Variable",
    yaxis_title="Aislado",
    template=TEMPLATE,
    width=650, height=max(400, 14 * len(perfil)),
    yaxis=dict(tickfont=dict(size=9)),
)
save(fig, "heatmap_susceptibilidad.html")


# ─── 3h. Distribución de %INH ─────────────────────────────────────
fig = make_subplots(rows=1, cols=3, subplot_titles=[
    "Maceración", "Soxhlet", "Ultrasonido"
])

for i, met in enumerate(ORDER_MET):
    sub = inh_5[inh_5["metodo_extraccion"] == met]["porcentaje_inhibicion"]
    fig.add_trace(go.Histogram(
        x=sub, nbinsx=15, name=LABEL_MET[met],
        marker_color=COLOR_MET[met], opacity=0.75,
        showlegend=False,
        hovertemplate="%INH: %{x:.1f}%<br>Frecuencia: %{y}<extra></extra>",
    ), row=1, col=i + 1)

fig.update_layout(
    title="Distribución de %INH por método (5.0 mg/mL)",
    template=TEMPLATE,
    width=900, height=400,
)
fig.update_xaxes(title_text="Inhibición (%)")
fig.update_yaxes(title_text="Frecuencia")
save(fig, "distribucion_inh_metodo.html", width=900, height=400)


# ─── 3i. Perfil de aislados en Maceración ──────────────────────────
mac_perfil = inh[inh["metodo_extraccion"] == "maceracion"].pivot_table(
    index="aislado_id", columns="concentracion_mg_ml",
    values="porcentaje_inhibicion", aggfunc="mean"
)
conc_order = sorted(mac_perfil.columns)

fig = go.Figure()
for ais in mac_perfil.index:
    fig.add_trace(go.Scatter(
        x=[str(c) for c in conc_order],
        y=mac_perfil.loc[ais],
        mode="lines+markers", name=str(ais),
        line=dict(width=1), marker=dict(size=4),
        opacity=0.5,
        hovertemplate=f"<b>{ais}</b><br>%{{x}} mg/mL<br>%INH: %{{y:.1f}}%<extra></extra>",
    ))

# Línea de la media
mean_line = mac_perfil.mean()
fig.add_trace(go.Scatter(
    x=[str(c) for c in conc_order],
    y=mean_line,
    mode="lines+markers", name="Media global",
    line=dict(color="black", width=3), marker=dict(size=8, color="black"),
))

fig.update_layout(
    title="Perfil individual de aislados — Maceración (media por concentración)",
    xaxis_title="Concentración (mg/mL)",
    yaxis_title="Inhibición media (%)",
    template=TEMPLATE,
    width=800, height=500,
    legend=dict(font=dict(size=7), y=1, itemsizing="constant"),
)
save(fig, "perfil_aislados_maceracion.html")


# ─── 3j. Ranking de susceptibilidad ────────────────────────────────
rank_plot = ranking.sort_values("rank").head(20)

fig = go.Figure()
fig.add_trace(go.Bar(
    y=rank_plot["aislado"], x=rank_plot["score_susceptibilidad"],
    orientation="h",
    marker=dict(
        color=rank_plot["score_susceptibilidad"],
        colorscale="RdYlGn", reversescale=True,
    ),
    text=rank_plot["score_susceptibilidad"].round(2),
    textposition="outside",
    hovertemplate="<b>%{y}</b><br>Score: %{x:.2f}<br>"
                  f"Clasificación: %{{customdata}}<extra></extra>",
    customdata=rank_plot["clasificacion"],
))

fig.update_layout(
    title="Ranking de susceptibilidad (top 20 aislados)",
    xaxis_title="Score de susceptibilidad",
    yaxis_title="",
    template=TEMPLATE,
    width=700, height=500,
    margin=dict(l=120, r=40),
)
save(fig, "ranking_susceptibilidad.html")


# ═══════════════════════════════════════════════════════════════════
print(f"\n🎉 {DIR_INTERACTIVO} — {len(list(DIR_INTERACTIVO.glob('*.html')))} figuras interactivas generadas.")
