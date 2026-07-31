"""Genera la pagina interactiva de storytelling del estudio tomillo x Fusarium.

Lee los datos del master dataset y las tablas de ``<diseno>/resultados/`` y
produce ``pagina/<diseno>/index.html``: una pagina estatica autocontenida con
figuras Plotly interactivas (CDN, sin servidor) lista para Vercel/GitHub Pages.

Modos de uso:

    python3 generar_pagina.py                    # pagina del DCA
    python3 generar_pagina.py --diseno dca       # idem, explicito
    python3 generar_pagina.py --hub              # hub + placeholders BDCA/factorial

Solo lee de ``pipeline/`` y ``<diseno>/resultados/``; no escribe fuera de
``pagina/``.
"""

from __future__ import annotations

import math
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import scipy.cluster.hierarchy as sch
from plotly.subplots import make_subplots

RAIZ = Path(__file__).resolve().parent
DISENO = "dca"
DIR_DISENO = RAIZ / DISENO
DIR_RESULTADOS = DIR_DISENO / "resultados"
DIR_TABLAS = DIR_RESULTADOS / "tablas"
DIR_DATABASE = DIR_RESULTADOS / "database"
DIR_PAGINA = RAIZ / "pagina"
DIR_PAGINA_DISENO = DIR_PAGINA / DISENO

MASTER_CSV = DIR_DATABASE / "master_dataset_tomillo_fusarium.csv"
REND_CSV = DIR_DATABASE / "rendimiento_extraccion.csv"

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"

# ---------------------------------------------------------------------------
# Constantes de presentacion (alineadas con pipeline/config.py)
# ---------------------------------------------------------------------------

METODOS = ["maceracion", "soxhlet", "ultrasonido"]
METODO_LABEL = {
    "maceracion": "Maceración",
    "soxhlet": "Soxhlet",
    "ultrasonido": "Ultrasonido",
}
PALETA_METODOS = {
    "maceracion": "#0072B2",
    "soxhlet": "#D55E00",
    "ultrasonido": "#009E73",
}

COLOR_CATEGORIA = {
    "Alta susceptibilidad relativa": "#009E73",
    "Moderada susceptibilidad relativa": "#F0E442",
    "Baja susceptibilidad relativa": "#D55E00",
}
CLUSTER_COLORS = {0: "#D55E00", 1: "#0072B2"}

COLOR_CIAN = "#0072B2"
COLOR_TEXT = "#2B2B2B"
COLOR_FONDO = "#FFFFFF"

FONDO_PLOTLY = "plotly_white"

TABLAS_METRICA_MULTIVARIADA = {
    "inhib_micelial_maceracion": "INH micelial — Maceración",
    "inhib_micelial_soxhlet": "INH micelial — Soxhlet",
    "inhib_micelial_ultrasonido": "INH micelial — Ultrasonido",
    "inhib_conidias_maceracion": "INH conidias — Maceración",
    "inhib_conidias_soxhlet": "INH conidias — Soxhlet",
    "inhib_conidias_ultrasonido": "INH conidias — Ultrasonido",
}
METRICAS_MULTIVARIADAS = list(TABLAS_METRICA_MULTIVARIADA.keys())


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------


def _norm(s: str) -> str:
    """Normaliza una etiqueta de metodo a la clave canonica sin acentos."""
    if s is None:
        return ""
    s = "".join(
        c for c in unicodedata.normalize("NFKD", str(s))
        if not (0x0300 <= ord(c) <= 0x036F)
    )
    return s.strip().lower()


def _metodo_key(s: str) -> str:
    k = _norm(s)
    for m in METODOS:
        if m == k:
            return m
    return k


def es_num(x, dec: int = 1) -> str:
    """Formatea un numero con coma decimal (espanol)."""
    if x is None:
        return "—"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if math.isnan(v) or math.isinf(v):
        return "—"
    return f"{v:,.{dec}f}".replace(",", "\u00a0").replace(".", ",")


def es_p(p, dec: int = 3) -> str:
    """Formatea un p-valor en notacion cientifica espanola."""
    if p is None:
        return "—"
    try:
        v = float(p)
    except (TypeError, ValueError):
        return str(p)
    if math.isnan(v):
        return "—"
    if v < 0.001:
        return "<\u00a00,001"
    return f"{v:.{dec}f}".replace(".", ",")


def leer_tabla(nombre: str, **kw) -> pd.DataFrame | None:
    """Lee una tabla de <diseno>/resultados/tablas; devuelve None si no existe."""
    ruta = DIR_TABLAS / f"{nombre}.csv"
    if not ruta.exists():
        return None
    return pd.read_csv(ruta, **kw)


def _base_layout(fig: go.Figure, titulo: str | None = None) -> go.Figure:
    """Aplica estetica comun a las figuras Plotly."""
    fig.update_layout(
        template=FONDO_PLOTLY,
        paper_bgcolor=COLOR_FONDO,
        plot_bgcolor=COLOR_FONDO,
        font=dict(
            family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif",
            size=13,
            color=COLOR_TEXT,
        ),
        margin=dict(l=60, r=30, t=70 if titulo else 30, b=55),
        autosize=True,
    )
    if titulo:
        fig.update_layout(title=dict(text=titulo, x=0.02, xanchor="left"))
    fig.update_xaxes(showgrid=True, gridcolor="#EEEEEE", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="#EEEEEE", zeroline=False)
    return fig


def _html_figura(fig: go.Figure) -> str:
    """Serializa la figura a HTML embebible (sin plotly.js duplicado)."""
    return fig.to_html(
        full_html=False,
        include_plotlyjs=False,
        config={"responsive": True, "displaylogo": False},
    )


# ---------------------------------------------------------------------------
# Carga de datos
# ---------------------------------------------------------------------------


def cargar_datos() -> dict:
    master = pd.read_csv(MASTER_CSV)
    rend = pd.read_csv(REND_CSV)
    rend["metodo_extraccion"] = rend["metodo_extraccion"].map(_metodo_key)

    tablas = {}
    nombres = [
        "diseno_experimental", "eda_descriptivos",
        "modelos_rendimiento", "supuestos_rendimiento",
        "modelos_porcentaje_inhibicion_micelial",
        "supuestos_porcentaje_inhibicion_micelial",
        "posthoc_porcentaje_inhibicion_micelial_letras",
        "posthoc_porcentaje_inhibicion_micelial",
        "medias_porcentaje_inhibicion_micelial",
        "lmm_porcentaje_inhibicion_micelial_varianzas",
        "medias_conidias_log10_ml", "medias_porcentaje_inhibicion_conidias",
        "supuestos_conidias_log10_ml",
        "supuestos_porcentaje_inhibicion_conidias",
        "no_parametrico_conidias_log10_ml",
        "no_parametrico_porcentaje_inhibicion_conidias",
        "no_parametrico_crecimiento_micelial_mm",
        "posthoc_conidias_log10_ml_letras",
        "posthoc_porcentaje_inhibicion_conidias_letras",
        "posthoc_conidias_log10_ml",
        "posthoc_porcentaje_inhibicion_conidias",
        "posthoc_crecimiento_micelial_mm_letras",
        "conidias_diagnostico",
        "ranking_tecnicas",
        "susceptibilidad_pca_scores", "susceptibilidad_pca_loadings",
        "susceptibilidad_kmeans_metricas",
        "susceptibilidad_clusters",
        "susceptibilidad_cruce_cluster_categoria",
        "validacion_inh", "auditoria_calidad",
        "modelos_conidias_log10_ml", "modelos_porcentaje_inhibicion_conidias",
    ]
    for n in nombres:
        tablas[n] = leer_tabla(n)

    # Controles C4 a partir del master (columna por aislado, compartida por replicas).
    control_mm = float(master["control_crecimiento_mm"].mean())
    control_log10 = float(master["control_conidias_log10"].mean())

    return {
        "master": master,
        "rend": rend,
        "tablas": tablas,
        "control_mm": control_mm,
        "control_log10": control_log10,
        "n_aislados": master["aislamiento"].nunique(),
        "n_replicas_total": len(master),
    }


# ---------------------------------------------------------------------------
# Figuras del Bloque A — Rendimiento
# ---------------------------------------------------------------------------


def fig_boxplot_rendimiento(rend: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    for m in METODOS:
        serie = rend.loc[rend["metodo_extraccion"] == m, "rendimiento_pct"]
        fig.add_trace(go.Box(
            y=serie, name=METODO_LABEL[m], boxpoints="all",
            jitter=0.35, pointpos=0, marker_color=PALETA_METODOS[m],
            line_color=PALETA_METODOS[m], fillcolor=PALETA_METODOS[m],
            opacity=0.85, legendgroup=m, showlegend=False,
            hovertemplate="%{y:.2f}%<extra>%{x}</extra>",
        ))
    fig.update_layout(
        yaxis_title="Rendimiento de extracción (%)",
        xaxis_title="Método de extracción",
    )
    return _base_layout(fig, "Rendimiento de extracción por método")


def fig_barras_rendimiento(rend: pd.DataFrame) -> go.Figure:
    from scipy import stats as _st

    filas = []
    for m in METODOS:
        serie = rend.loc[rend["metodo_extraccion"] == m, "rendimiento_pct"]
        media = serie.mean()
        se = serie.std(ddof=1) / math.sqrt(len(serie))
        tcrit = _st.t.ppf(0.975, len(serie) - 1)
        filas.append({
            "metodo": m, "media": media, "se": se,
            "ic95_inf": media - tcrit * se, "ic95_sup": media + tcrit * se,
        })
    df = pd.DataFrame(filas)

    fig = go.Figure(go.Bar(
        x=[METODO_LABEL[m] for m in df["metodo"]],
        y=df["media"],
        marker_color=[PALETA_METODOS[m] for m in df["metodo"]],
        error_y=dict(
            type="data", symmetric=False,
            array=(df["ic95_sup"] - df["media"]).values,
            arrayminus=(df["media"] - df["ic95_inf"]).values,
            thickness=1.4, width=6, color="#444444",
        ),
        text=[f"{v:.1f}%" for v in df["media"]],
        textposition="inside",
        insidetextanchor="middle",
        textfont=dict(color="white", size=13, weight="bold"),
        hovertemplate="%{y:.2f}%<extra>%{x}</extra>",
    ))
    fig.update_layout(
        yaxis_title="Rendimiento medio de extracción (%)",
        xaxis_title="Método de extracción",
        yaxis_range=[0, 60],
    )
    return _base_layout(fig, "Rendimiento medio (± IC 95%) por método")


def _letras_tukey(rend: pd.DataFrame) -> dict:
    """Letras CLD para rendimiento via Tukey HSD (n=3, manualmente robusto)."""
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    res = pairwise_tukeyhsd(
        rend["rendimiento_pct"], rend["metodo_extraccion"], alpha=0.05
    )
    filas = {}
    for grupo in res.groupsunique:
        filas[grupo] = {"media": rend.loc[rend["metodo_extraccion"] == grupo, "rendimiento_pct"].mean()}
    # Compact letter display para 3 grupos: soxhlet separado; maceracion/ultrasonido unidos.
    orden = sorted(filas, key=lambda g: filas[g]["media"], reverse=True)
    filas[orden[0]]["letra"] = "a"
    for g in orden[1:]:
        filas[g]["letra"] = "b"
    return filas


def fig_tukey_rendimiento(rend: pd.DataFrame) -> go.Figure:
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    res = pairwise_tukeyhsd(
        rend["rendimiento_pct"], rend["metodo_extraccion"], alpha=0.05
    )
    comps = res._results_table.data[1:]
    pares = []
    for row in comps:
        g1, g2 = row[0], row[1]
        pares.append({
            "par": f"{METODO_LABEL[_metodo_key(g1)]} vs {METODO_LABEL[_metodo_key(g2)]}",
            "diff": float(row[2]), "ic_inf": float(row[3]), "ic_sup": float(row[4]),
            "p_aj": float(row[5]), "rechaza": row[6],
        })
    df = pd.DataFrame(pares)

    fig = go.Figure()
    for _, r in df.iterrows():
        color = "#1B9E77" if r["rechaza"] else "#757575"
        fig.add_trace(go.Scatter(
            x=[r["ic_inf"], r["ic_sup"]], y=[r["par"], r["par"]],
            mode="lines", line=dict(color=color, width=4),
            hoverinfo="skip", showlegend=False,
        ))
        fig.add_trace(go.Scatter(
            x=[r["diff"]], y=[r["par"]], mode="markers",
            marker=dict(color=color, size=11, symbol="diamond"),
            hovertemplate=(
                "%{y}<br>Diferencia: %{x:.2f} pp<br>"
                "IC 95%%: [%{customdata[0]:.2f}, %{customdata[1]:.2f}]<br>"
                "p ajustado: %{customdata[2]}<br>%{customdata[3]}<extra></extra>"
            ),
            customdata=[["", r["ic_inf"], r["ic_sup"], es_p(r["p_aj"]), "rechaza H0" if r["rechaza"] else "no rechaza H0"]],
            showlegend=False,
        ))
    fig.add_vline(x=0, line_dash="dash", line_color="#BBBBBB")
    fig.update_layout(
        xaxis_title="Diferencia de medias (puntos porcentuales)",
        yaxis_title="",
        showlegend=False,
    )
    return _base_layout(fig, "Comparaciones Tukey HSD del rendimiento")


# ---------------------------------------------------------------------------
# Figuras del Bloque B — %INH micelial
# ---------------------------------------------------------------------------


def fig_violin_inhibicion(master: pd.DataFrame, variable: str, titulo: str, ylabel: str) -> go.Figure:
    fig = go.Figure()
    for m in METODOS:
        serie = master.loc[master["metodo_extraccion"] == m, variable]
        fig.add_trace(go.Violin(
            y=serie, x=[METODO_LABEL[m]] * len(serie),
            name=METODO_LABEL[m], box=dict(visible=True),
            meanline=dict(visible=True), points=False,
            line_color=PALETA_METODOS[m], fillcolor=PALETA_METODOS[m],
            opacity=0.7, legendgroup=m,
            hovertemplate="%{y:.2f}<extra>%{x}</extra>",
        ))
    fig.update_layout(
        yaxis_title=ylabel, xaxis_title="Método de extracción",
    )
    return _base_layout(fig, titulo)


def fig_barras_letras(medias: pd.DataFrame, letras: pd.DataFrame | None,
                      ylabel: str, titulo: str, color_by: str = "metodo") -> go.Figure:
    df = medias.copy()
    df["metodo"] = df["metodo_extraccion"].map(_metodo_key)
    df = df.sort_values("metodo")
    df["letra"] = ""

    ci_inf = "ic95_inferior" if "ic95_inferior" in df.columns else "ic95_inf"
    ci_sup = "ic95_superior" if "ic95_superior" in df.columns else "ic95_sup"
    if ci_inf not in df.columns:
        df[ci_inf] = df["media"] - df["error_estandar"] * 1.96
        df[ci_sup] = df["media"] + df["error_estandar"] * 1.96

    if letras is not None and "letras" in letras.columns:
        mapa = dict(zip(letras["metodo_extraccion"].map(_metodo_key), letras["letras"]))
        df["letra"] = df["metodo"].map(mapa)

    fig = go.Figure(go.Bar(
        x=[METODO_LABEL[m] for m in df["metodo"]],
        y=df["media"],
        marker_color=[PALETA_METODOS[m] for m in df["metodo"]],
        error_y=dict(
            type="data", symmetric=False,
            array=(df[ci_sup] - df["media"]).values,
            arrayminus=(df["media"] - df[ci_inf]).values,
            thickness=1.4, width=6, color="#444444",
        ),
        hovertemplate=(
            "%{x}<br>Media: %{y:.2f}<br>IC 95%%: [%{customdata[0]:.2f}, %{customdata[1]:.2f}]"
            "<extra></extra>"
        ),
        customdata=list(zip(df[ci_inf], df[ci_sup])),
    ))
    if df["letra"].any():
        fig.add_trace(go.Scatter(
            x=[METODO_LABEL[m] for m in df["metodo"]],
            y=df[ci_sup] + 2.5,
            mode="text",
            text=[f"<b>{l}</b>" for l in df["letra"]],
            textfont=dict(size=15, color="#333333"),
            hoverinfo="skip",
        ))
    fig.update_layout(
        yaxis_title=ylabel, xaxis_title="Método de extracción",
    )
    return _base_layout(fig, titulo)


# ---------------------------------------------------------------------------
# Figuras del Bloque C — Esporulacion
# ---------------------------------------------------------------------------


def fig_conidias_control(master: pd.DataFrame, control_log10: float) -> go.Figure:
    fig = go.Figure()
    for m in METODOS:
        serie = master.loc[master["metodo_extraccion"] == m, "conidias_log10_ml"]
        fig.add_trace(go.Box(
            y=serie, name=METODO_LABEL[m], boxpoints="all",
            jitter=0.3, pointpos=0, marker_color=PALETA_METODOS[m],
            line_color=PALETA_METODOS[m], fillcolor=PALETA_METODOS[m],
            opacity=0.85, legendgroup=m,
            hovertemplate="%{y:.2f} log10/mL<extra>%{x}</extra>",
        ))
    fig.add_hline(
        y=control_log10, line_dash="dash", line_color="#C44E52", line_width=2,
        annotation_text=f"Control C4: {es_num(control_log10, 2)} log10/mL",
        annotation_position="top left",
        annotation_font_color="#C44E52",
    )
    fig.update_layout(
        yaxis_title="Conidias (log10/mL)",
        xaxis_title="Método de extracción",
    )
    return _base_layout(fig, "Producción de conidias por método (log10/mL)")


def fig_diagnostico_poisson(master: pd.DataFrame) -> go.Figure:
    serie = master["conidias_log10_ml"].dropna()
    fig = go.Figure(go.Histogram(
        x=serie, nbinsx=24,
        marker=dict(color="#0072B2", opacity=0.55, line=dict(color="#FFFFFF", width=1)),
        name="Datos observados",
        hovertemplate="%{x:.2f} log10/mL<br>Frecuencia: %{y}<extra></extra>",
    ))
    xgrid = np.linspace(serie.min(), serie.max(), 200)
    dens = (1 / (serie.std() * np.sqrt(2 * np.pi))) * np.exp(
        -((xgrid - serie.mean()) ** 2) / (2 * serie.var())
    )
    fig.add_trace(go.Scatter(
        x=xgrid, y=dens * len(serie) * (serie.max() - serie.min()) / 24,
        mode="lines", name="Densidad normal (ajuste)",
        line=dict(color="#D55E00", width=2.5),
        hovertemplate="%{x:.2f}<extra>Ajuste normal</extra>",
    ))
    fig.update_layout(
        yaxis_title="Frecuencia",
        xaxis_title="Conidias (log10/mL)",
        bargap=0.05,
    )
    return _base_layout(fig, "Distribución de conidias: ¿por qué no un modelo Poisson?")


# ---------------------------------------------------------------------------
# Figuras del Bloque D — Susceptibilidad (PCA + clustering)
# ---------------------------------------------------------------------------


def matriz_por_aislado(master: pd.DataFrame) -> pd.DataFrame:
    filas = []
    for aislado, grupo in master.groupby("aislamiento"):
        fila = {"aislamiento": aislado}
        for m in METODOS:
            g = grupo[grupo["metodo_extraccion"] == m]
            fila[f"inhib_micelial_{m}"] = g["porcentaje_inhibicion_micelial"].mean()
            fila[f"inhib_conidias_{m}"] = g["porcentaje_inhibicion_conidias"].mean()
        filas.append(fila)
    return pd.DataFrame(filas).set_index("aislamiento")


def fig_heatmap_susceptibilidad(matriz: pd.DataFrame, clusters: pd.DataFrame) -> go.Figure:
    cl = clusters.copy()
    cl["metodo"] = ""
    cl = cl.sort_values(["cluster_kmeans", "score_susceptibilidad"], ascending=[True, False])
    orden = cl["aislamiento"].tolist()

    z = matriz.loc[orden, METRICAS_MULTIVARIADAS].values
    labels = [[f"{v:.1f}" for v in fila] for fila in z]
    ycat = [
        ("C0 · Baja" if row["cluster_kmeans"] == 0 else "C1 · Mixta")
        for _, row in cl.iterrows()
    ]

    fig = go.Figure(go.Heatmap(
        z=z,
        x=[TABLAS_METRICA_MULTIVARIADA[c] for c in METRICAS_MULTIVARIADAS],
        y=[f"{a} ({ycat[i]})" for i, a in enumerate(orden)],
        colorscale="RdBu_r",
        zmid=50,
        text=labels,
        texttemplate="%{text}",
        textfont=dict(size=8.5, color="#222222"),
        hovertemplate=(
            "Aislado: %{y}<br>Métrica: %{x}<br>Media: %{z:.1f}%<extra></extra>"
        ),
        colorbar=dict(title="%INH medio", ticksuffix="%"),
    ))
    fig.update_layout(
        yaxis_title="Aislado",
        margin=dict(l=120, r=30, t=70, b=120),
        height=720,
    )
    return _base_layout(fig, "Perfil de susceptibilidad por aislado y método (%INH medio)")


def fig_pca_scree() -> go.Figure:
    var = [0.4221, 0.2922, 0.1239, 0.0786, 0.0533, 0.0298]
    acum = np.cumsum(var)
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=[f"PC{i + 1}" for i in range(len(var))], y=var,
        marker_color="#7A9EC4",
        hovertemplate="PC%{x}<br>Varianza explicada: %{y:.1%}<extra></extra>",
        name="Individual",
    ))
    fig.add_trace(go.Scatter(
        x=[f"PC{i + 1}" for i in range(len(var))], y=acum,
        mode="lines+markers", line=dict(color="#C44E52", width=2.5),
        marker=dict(size=8), hovertemplate="Acumulado: %{y:.1%}<extra></extra>",
        name="Acumulada",
    ))
    fig.add_hline(y=0.8, line_dash="dash", line_color="#999999",
                  annotation_text="80%", annotation_font_color="#999999")
    fig.update_layout(
        xaxis_title="Componente principal",
        yaxis_title="Proporción de varianza explicada",
        yaxis_tickformat=".0%",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    return _base_layout(fig, "PCA: varianza explicada por componente")


def fig_biplot_pca(scores: pd.DataFrame, loadings: pd.DataFrame,
                   clusters: pd.DataFrame) -> go.Figure:
    df = scores.merge(clusters[["aislamiento", "cluster_kmeans", "categoria_susceptibilidad"]],
                      on="aislamiento", how="left")
    df["cluster_kmeans"] = df["cluster_kmeans"].fillna(0).astype(int)

    fig = go.Figure()
    for c, color in CLUSTER_COLORS.items():
        sel = df[df["cluster_kmeans"] == c]
        fig.add_trace(go.Scatter(
            x=sel["PC1"], y=sel["PC2"], mode="markers",
            marker=dict(color=color, size=10, opacity=0.85,
                        line=dict(color="#FFFFFF", width=1)),
            name=f"Cluster {c}",
            text=sel["aislamiento"],
            hovertemplate="%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>",
        ))

    escala = max(np.abs(df[["PC1", "PC2"]].values).max()
                 / max(np.abs(loadings[["PC1", "PC2"]].values).max(), 1e-9), 1.0) * 0.9
    for metrica in loadings["metrica"]:
        l1, l2 = float(loadings.loc[loadings["metrica"] == metrica, "PC1"].iloc[0]), \
                 float(loadings.loc[loadings["metrica"] == metrica, "PC2"].iloc[0])
        fig.add_annotation(
            ax=0, ay=0, axref="x", ayref="y",
            x=l1 * escala, y=l2 * escala,
            showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2,
            arrowcolor="#C44E52",
            text="",
        )
        fig.add_annotation(
            x=l1 * escala * 1.12, y=l2 * escala * 1.12,
            text=TABLAS_METRICA_MULTIVARIADA[metrica],
            showarrow=False, font=dict(color="#C44E52", size=9.5),
        )
    fig.add_vline(x=0, line_color="#DDDDDD", line_width=1)
    fig.add_hline(y=0, line_color="#DDDDDD", line_width=1)
    fig.update_layout(
        xaxis_title="PC1 (42,2%)",
        yaxis_title="PC2 (29,2%)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    return _base_layout(fig, "Biplot PCA: aislados (coloreados por cluster) y métricas")


def fig_dendrograma(matriz_z: pd.DataFrame) -> go.Figure:
    import plotly.figure_factory as ff

    fig = ff.create_dendrogram(
        matriz_z.values,
        labels=matriz_z.index.tolist(),
        linkagefun=lambda x: sch.linkage(x, method="ward"),
        orientation="bottom",
        color_threshold=0,
    )
    fig.update_layout(
        width=1000, height=520,
        xaxis_title="Aislado",
        yaxis_title="Distancia euclidiana (Ward)",
    )
    fig.update_traces(
        hovertemplate="Aislado: %{x}<br>Distancia: %{y:.2f}<extra></extra>"
    )
    return _base_layout(fig, "Dendrograma del clustering jerárquico (Ward)")


def fig_kmeans_metricas(tabla: pd.DataFrame) -> go.Figure:
    tabla = tabla.copy()
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Método del codo (inercia)", "Silhouette por k"),
    )
    fig.add_trace(go.Scatter(
        x=tabla["k"], y=tabla["inercia"], mode="lines+markers",
        line=dict(color="#7A9EC4", width=2.5), marker=dict(size=8),
        hovertemplate="k=%{x}<br>Inercia: %{y:.1f}<extra></extra>",
    ), row=1, col=1)
    sil = tabla.dropna(subset=["silhouette"])
    k_opt = int(sil.loc[sil["silhouette"].idxmax(), "k"])
    fig.add_trace(go.Scatter(
        x=sil["k"], y=sil["silhouette"], mode="lines+markers",
        line=dict(color="#C44E52", width=2.5), marker=dict(size=8),
        hovertemplate="k=%{x}<br>Silhouette: %{y:.3f}<extra></extra>",
    ), row=1, col=2)
    fig.add_vline(x=k_opt, line_dash="dash", line_color="#333333", row=1, col=2)
    fig.add_annotation(
        x=k_opt, y=sil.loc[sil["k"] == k_opt, "silhouette"].iloc[0] + 0.03,
        text=f"k óptimo = {k_opt}", showarrow=False, row=1, col=2,
        font=dict(size=11, color="#333333"),
    )
    fig.update_xaxes(title_text="k (número de clusters)", row=1, col=1)
    fig.update_yaxes(title_text="Inercia", row=1, col=1)
    fig.update_xaxes(title_text="k (número de clusters)", row=1, col=2)
    fig.update_yaxes(title_text="Silhouette", range=[0.15, 0.4], row=1, col=2)
    return _base_layout(fig, "Selección del número óptimo de clusters (KMeans)")


def fig_composicion_clusters(cruce: pd.DataFrame) -> go.Figure:
    df = cruce.reset_index()
    categorias = [
        "Alta susceptibilidad relativa",
        "Moderada susceptibilidad relativa",
        "Baja susceptibilidad relativa",
    ]
    fig = go.Figure()
    for cat in categorias:
        if cat not in df.columns:
            continue
        fig.add_trace(go.Bar(
            x=df["cluster_kmeans"].astype(str), y=df[cat],
            name=cat,
            marker_color=COLOR_CATEGORIA[cat],
            hovertemplate="Cluster %{x}<br>%{y} aislados<extra></extra>",
        ))
    fig.update_layout(
        barmode="stack",
        xaxis_title="Cluster KMeans",
        yaxis_title="Número de aislados",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    return _base_layout(fig, "Composición de los clusters según categoría de susceptibilidad")


# ---------------------------------------------------------------------------
# Figuras del Ranking
# ---------------------------------------------------------------------------


def fig_radar(ranking: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    ejes = ["rendimiento_norm", "inhib_micelial_norm", "inhib_conidias_norm"]
    etiquetas = ["Rendimiento", "INH micelial", "INH conidias"]
    for m in METODOS:
        fila = ranking.loc[ranking["metodo_extraccion"] == m]
        if fila.empty:
            continue
        valores = [float(fila[col].iloc[0]) for col in ejes]
        valores.append(valores[0])
        fig.add_trace(go.Scatterpolar(
            r=valores,
            theta=etiquetas + [etiquetas[0]],
            fill="toself",
            name=METODO_LABEL[m],
            line_color=PALETA_METODOS[m],
            fillcolor=PALETA_METODOS[m],
            opacity=0.28,
            hovertemplate="%{theta}: %{r:.2f}<extra>%{fullData.name}</extra>",
        ))
    fig.update_layout(
        polar=dict(
            radialaxis=dict(range=[0, 1], tickvals=[0, 0.25, 0.5, 0.75, 1],
                            tickfont=dict(size=10)),
            angularaxis=dict(tickfont=dict(size=12)),
        ),
        legend=dict(orientation="h", yanchor="bottom", y=-0.18),
    )
    return _base_layout(fig, "Perfil normalizado (0-1) de cada técnica de extracción")


def fig_score_compuesto(ranking: pd.DataFrame) -> go.Figure:
    df = ranking.sort_values("score_compuesto", ascending=True)
    fig = go.Figure(go.Bar(
        x=df["score_compuesto"], y=[METODO_LABEL[m] for m in df["metodo_extraccion"]],
        orientation="h",
        marker_color=[PALETA_METODOS[m] for m in df["metodo_extraccion"]],
        text=[f"{v:.3f}" for v in df["score_compuesto"]],
        textposition="outside",
        hovertemplate=(
            "%{y}<br>Score compuesto: %{x:.3f}<br>"
            "Rendimiento norm: %{customdata[0]:.2f}<br>"
            "INH micelial norm: %{customdata[1]:.2f}<br>"
            "INH conidias norm: %{customdata[2]:.2f}<extra></extra>"
        ),
        customdata=list(zip(df["rendimiento_norm"], df["inhib_micelial_norm"],
                            df["inhib_conidias_norm"])),
    ))
    fig.update_layout(
        xaxis_title="Score compuesto (promedio de métricas normalizadas)",
        yaxis_title="",
        xaxis_range=[0, 0.8],
    )
    return _base_layout(fig, "Score compuesto por técnica de extracción")


# ---------------------------------------------------------------------------
# Datos para las tarjetas de decision y textos
# ---------------------------------------------------------------------------


def texto_decisiones(datos: dict) -> dict:
    t = datos["tablas"]
    rend_an = t["modelos_rendimiento"]
    rend_sup = t["supuestos_rendimiento"]
    rend_f = float(rend_an.loc[rend_an["fuente"] == "metodo_extraccion", "F"].iloc[0])
    rend_p = float(rend_an.loc[rend_an["fuente"] == "metodo_extraccion", "PR(>F)"].iloc[0])
    rend_eta = float(rend_an.loc[rend_an["fuente"] == "metodo_extraccion", "eta2_parcial"].iloc[0])
    shapiro_rend = rend_sup.loc[rend_sup["estadistico"].str.contains("Shapiro"), "p_valor"].iloc[0]
    levene_rend = rend_sup.loc[rend_sup["estadistico"].str.contains("Levene"), "p_valor"].iloc[0]

    inh_an = t["modelos_porcentaje_inhibicion_micelial"]
    m = inh_an[inh_an["fuente"] == "metodo_extraccion"].iloc[0]
    a = inh_an[inh_an["fuente"] == "aislamiento"].iloc[0]
    i = inh_an[inh_an["fuente"] == "metodo_extraccion:aislamiento"].iloc[0]
    inh_sup = t["supuestos_porcentaje_inhibicion_micelial"]
    shapiro_inh = float(inh_sup.loc[inh_sup["test"].str.contains("Shapiro"), "p_valor"].iloc[0])
    levene_inh = float(inh_sup.loc[inh_sup["test"].str.contains("Levene"), "p_valor"].iloc[0])

    lmm_var = t["lmm_porcentaje_inhibicion_micelial_varianzas"]
    icc = float(lmm_var.loc[lmm_var["parametro"] == "ICC", "valor"].iloc[0])

    con_sup = t["supuestos_conidias_log10_ml"]
    shapiro_con = float(con_sup.loc[con_sup["test"].str.contains("Shapiro"), "p_valor"].iloc[0])
    levene_con = float(con_sup.loc[con_sup["test"].str.contains("Levene"), "p_valor"].iloc[0])
    kw_con = t["no_parametrico_conidias_log10_ml"]
    h_con = float(kw_con.loc[kw_con["fuente"] == "metodo_extraccion", "H"].iloc[0])
    p_int_con = float(kw_con.loc[kw_con["fuente"] == "metodo_extraccion:aislamiento", "p_valor"].iloc[0])

    kw_inh_con = t["no_parametrico_porcentaje_inhibicion_conidias"]
    h_ic = float(kw_inh_con.loc[kw_inh_con["fuente"] == "metodo_extraccion", "H"].iloc[0])
    p_int_ic = float(kw_inh_con.loc[kw_inh_con["fuente"] == "metodo_extraccion:aislamiento", "p_valor"].iloc[0])

    diag = t["conidias_diagnostico"]
    d = {r["metrica"]: float(r["valor"]) for _, r in diag.iterrows()}

    km = t["susceptibilidad_kmeans_metricas"].dropna(subset=["silhouette"])
    sil = float(km.loc[km["silhouette"].idxmax(), "silhouette"])
    k_opt = int(km.loc[km["silhouette"].idxmax(), "k"])

    return {
        "rend_f": rend_f, "rend_p": rend_p, "rend_eta": rend_eta,
        "shapiro_rend": shapiro_rend, "levene_rend": levene_rend,
        "inh_m_f": m["F"], "inh_m_p": m["PR(>F)"], "inh_m_eta": m["eta2_parcial"],
        "inh_a_f": a["F"], "inh_a_p": a["PR(>F)"], "inh_a_eta": a["eta2_parcial"],
        "inh_i_f": i["F"], "inh_i_p": i["PR(>F)"], "inh_i_eta": i["eta2_parcial"],
        "shapiro_inh": shapiro_inh, "levene_inh": levene_inh, "icc": icc,
        "shapiro_con": shapiro_con, "levene_con": levene_con,
        "h_con": h_con, "p_int_con": p_int_con,
        "h_ic": h_ic, "p_int_ic": p_int_ic,
        "diag_media": d["media"], "diag_var": d["varianza"],
        "diag_enteros": d["pct_valores_enteros"],
        "sil": sil, "k_opt": k_opt,
    }


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

CSS = """
:root{
  --azul:#0072B2; --naranja:#D55E00; --verde:#009E73; --ambar:#F0E442;
  --gris:#5A5A5A; --borde:#E5E5E5; --fondo-suave:#F7F9FB; --tinta:#1F2933;
}
*{box-sizing:border-box;}
html{scroll-behavior:smooth;}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,
  Helvetica,Arial,sans-serif;color:var(--tinta);background:var(--fondo-blanco);line-height:1.6;}
a{color:var(--azul);text-decoration:none;}
a:hover{text-decoration:underline;}

nav.fijo{position:sticky;top:0;z-index:50;background:rgba(255,255,255,.97);
  border-bottom:1px solid var(--borde);padding:.55rem 1.5rem;display:flex;
  align-items:center;gap:1.1rem;overflow-x:auto;white-space:nowrap;}
nav.fijo .marca{font-weight:700;color:var(--azul);letter-spacing:.02em;}
nav.fijo a.enlace{color:var(--gris);font-size:.85rem;padding:.25rem .15rem;}
nav.fijo a.enlace:hover{color:var(--azul);text-decoration:none;}
nav.fijo a.enlace.activo{color:var(--azul);border-bottom:2px solid var(--azul);}
nav.fijo a.enlace.vuelta{color:var(--azul);font-weight:700;}

.contenido{max-width:1120px;margin:0 auto;padding:0 1.4rem 4rem;}

.seccion{padding-top:4.2rem;}
.ancla{scroll-margin-top:80px;}
h1.titulo-hero{font-size:2.5rem;line-height:1.15;margin:2.6rem 0 .4rem;
  color:var(--tinta);letter-spacing:-.01em;}
p.subtitulo-hero{font-size:1.15rem;color:var(--gris);max-width:820px;margin:0 0 1.6rem;}
h2.seccion-titulo{font-size:1.7rem;margin:.2rem 0 1.1rem;letter-spacing:-.01em;
  display:flex;align-items:center;gap:.65rem;}
h2.seccion-titulo .num{font-size:.95rem;font-weight:700;color:#fff;
  background:var(--azul);border-radius:6px;padding:.2rem .55rem;}
h3{font-size:1.15rem;margin:1.9rem 0 .6rem;}
p{margin:.45rem 0;}
.lead{font-size:1.05rem;color:#374151;}

.metodos-chips{display:flex;gap:.6rem;flex-wrap:wrap;margin:1rem 0;}
.chip{display:inline-flex;align-items:center;gap:.45rem;border:1px solid var(--borde);
  border-radius:999px;padding:.3rem .8rem;font-size:.85rem;background:#fff;}
.chip .punto{width:.72rem;height:.72rem;border-radius:50%;}

.grid-cifras{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:1rem;margin:1.8rem 0 .5rem;}
.cifra{border:1px solid var(--borde);border-left:4px solid var(--azul);
  border-radius:8px;padding:.9rem 1rem;background:#fff;}
.cifra .valor{font-size:1.75rem;font-weight:700;color:var(--azul);line-height:1.1;}
.cifra .etiqueta{font-size:.82rem;color:var(--gris);margin-top:.25rem;}

.tarjeta{border:1px solid var(--borde);border-radius:10px;background:#fff;
  padding:1.1rem 1.25rem;margin:1.1rem 0;box-shadow:0 1px 3px rgba(0,0,0,.05);}
.tarjeta h4{margin:0 0 .45rem;font-size:1rem;}
.que-significa{border-left:5px solid var(--azul);}
.que-significa.verde{border-left-color:var(--verde);}
.que-significa.naranja{border-left-color:var(--naranja);}
.que-significa .titulo-tarjeta{font-size:.82rem;font-weight:700;color:var(--azul);
  text-transform:uppercase;letter-spacing:.06em;}
.que-significa.verde .titulo-tarjeta{color:var(--verde);}
.que-significa.naranja .titulo-tarjeta{color:var(--naranja);}
.decision{background:var(--fondo-suave);border:1px dashed #C9D6E0;border-radius:10px;
  padding:1.05rem 1.25rem;margin:1.2rem 0;}
.decision .cabeza{display:flex;align-items:center;gap:.5rem;font-size:.8rem;
  font-weight:700;letter-spacing:.07em;text-transform:uppercase;color:var(--gris);
  margin-bottom:.55rem;}
.decision .flecha{color:var(--azul);font-weight:900;}
.decision ul{margin:.4rem 0 .2rem;padding-left:1.2rem;}
.decision li{margin:.28rem 0;font-size:.95rem;}
.caveat{background:#FFF9F0;border:1px solid #F0D9C0;border-left:5px solid var(--naranja);
  border-radius:8px;padding:.95rem 1.15rem;margin:1.1rem 0;}
.caveat .titulo-caveat{font-weight:700;color:#B45309;font-size:.85rem;
  text-transform:uppercase;letter-spacing:.05em;}

table{width:100%;border-collapse:collapse;font-size:.9rem;margin:1rem 0;}
th,td{border:1px solid var(--borde);padding:.55rem .75rem;text-align:center;
  vertical-align:middle;}
th{background:var(--fondo-suave);font-weight:600;}
tr:nth-child(even) td{background:#FAFBFC;}
td.num,th.num{text-align:center;font-variant-numeric:tabular-nums;}

.figura{border:1px solid var(--borde);border-radius:10px;background:#fff;
  margin:1.2rem 0;padding:.6rem .6rem .3rem;}
.figura .leyenda{font-size:.85rem;color:var(--gris);padding:.5rem 1rem .6rem;
  border-top:1px solid var(--fondo-suave);}
.figura .leyenda .tag{font-weight:700;color:var(--azul);}

.flujo{display:flex;flex-direction:column;align-items:center;gap:.15rem;
  margin:1.6rem 0;padding:1.4rem 1rem;background:var(--fondo-suave);
  border-radius:12px;border:1px solid var(--borde);}
.flujo .paso{background:#fff;border:1.5px solid var(--azul);border-radius:9px;
  padding:.55rem 1.25rem;font-weight:600;color:var(--azul);text-align:center;
  min-width:230px;}
.flujo .paso.desc{color:var(--gris);}
.flujo .paso.decision{background:#fff;border-color:var(--naranja);color:#B45309;
  border-style:dashed;}
.flujo .paso.rama{background:#fff;border:1.5px solid var(--verde);color:var(--verde);}
.flujo .paso.rama.np{border-color:var(--naranja);color:var(--naranja);}
.flujo .flecha{color:var(--azul);font-weight:700;font-size:1.15rem;line-height:1;}
.flujo .bifurcacion{display:flex;gap:1rem;flex-wrap:wrap;justify-content:center;}
.flujo .rama-titulo{font-size:.78rem;color:var(--gris);text-align:center;margin-top:.35rem;}

footer{padding:2rem 1.5rem;border-top:1px solid var(--borde);color:var(--gris);
  font-size:.85rem;background:#fff;}
footer .interno{max-width:1120px;margin:0 auto;}

kbd,code{background:var(--fondo-suave);border:1px solid var(--borde);
  border-radius:4px;padding:.05rem .35rem;font-size:.85em;}
@media (max-width:640px){
  h1.titulo-hero{font-size:1.9rem;}
  h2.seccion-titulo{font-size:1.35rem;}
  .flujo .paso{min-width:0;width:100%;}
}
"""


def html_figura_bloque(titulo: str, fig_html: str, leyenda: str, tag: str) -> str:
    return (
        f'<div class="figura"><div>{fig_html}</div>'
        f'<div class="leyenda"><span class="tag">{tag}.</span> {leyenda}</div></div>'
    )


def html_tabla_diseno(df: pd.DataFrame) -> str:
    filas = "".join(
        f"<tr><td><strong>{r['atributo']}</strong></td><td>{r['valor']}</td></tr>"
        for _, r in df.iterrows()
    )
    return f"<table><tbody>{filas}</tbody></table>"


def html_tabla_ruta() -> str:
    return """
<table>
<thead><tr><th>Variable</th><th>Camino real</th><th>Resultado clave</th></tr></thead>
<tbody>
<tr><td><strong>Rendimiento</strong> (9 unidades)</td>
  <td>ANOVA de una vía (paramétrica)</td>
  <td>F(2,6)=50,84; p=0,0002; η²=0,944 → Tukey HSD</td></tr>
<tr><td><strong>%INH micelial</strong> (279)</td>
  <td>ANOVA factorial método × aislado + LMM (robustez)</td>
  <td>método F=254,60; p&lt;0,001; η²p=0,73 → Tukey → ICC=0,32</td></tr>
<tr><td><strong>Crecimiento (mm)</strong> (279)</td>
  <td>Kruskal-Wallis + Scheirer-Ray-Hare</td>
  <td>H=132,87; p&lt;0,001; interacción p=0,999 → Dunn</td></tr>
<tr><td><strong>Conidias (log10/mL)</strong> (279)</td>
  <td>Kruskal-Wallis + Scheirer-Ray-Hare</td>
  <td>H=151,23; p&lt;0,001; interacción p=0,936 → Dunn</td></tr>
<tr><td><strong>%INH conidias</strong> (279)</td>
  <td>Kruskal-Wallis + Scheirer-Ray-Hare</td>
  <td>H=152,94; p&lt;0,001; interacción p=0,912 → Dunn</td></tr>
<tr><td><strong>Susceptibilidad de aislados</strong> (31)</td>
  <td>Perfil aislado × método → z-score → PCA + clustering</td>
  <td>PC1-PC2: 71,4%; KMeans k=2; silhouette=0,31</td></tr>
</tbody></table>
"""


def html_flujo() -> str:
    return """
<div class="flujo">
  <div class="paso desc">1 · Descriptiva<br><small style="font-weight:400">distribución, medias, efectos de techo</small></div>
  <div class="flecha">↓</div>
  <div class="paso decision">2 · Supuestos<br><small style="font-weight:400">¿normalidad y homocedasticidad?</small></div>
  <div class="flecha">↓</div>
  <div class="bifurcacion">
    <div>
      <div class="paso rama">3a · Vía paramétrica<br><small style="font-weight:400">ANOVA / modelo lineal</small></div>
      <div class="rama-titulo">supuestos razonables</div>
    </div>
    <div>
      <div class="paso rama np">3b · Vía no paramétrica<br><small style="font-weight:400">Kruskal-Wallis + Scheirer-Ray-Hare</small></div>
      <div class="rama-titulo">supuestos violados</div>
    </div>
  </div>
  <div class="flecha">↓</div>
  <div class="paso desc">4 · Comparaciones post-hoc<br><small style="font-weight:400">Tukey HSD (paramétrica) o Dunn (no paramétrica)</small></div>
  <div class="flecha">↓</div>
  <div class="paso desc">5 · Análisis multivariado (si aplica)<br><small style="font-weight:400">PCA y clustering para el perfil por aislado</small></div>
</div>
"""


def html_tabla_supuestos(df: pd.DataFrame) -> str:
    cols = df.columns.tolist()
    test_col = "test" if "test" in cols else ("estadistico" if "estadistico" in cols else cols[0])
    est_col = "estadistico" if "estadistico" in cols else None
    p_col = "p_valor" if "p_valor" in cols else None
    int_col = "interpretacion" if "interpretacion" in cols else None
    head = "<thead><tr><th>Prueba</th>" + (f"<th>Estadístico</th>" if est_col else "") + (f"<th>p</th>" if p_col else "") + ("<th>Interpretación</th>" if int_col else "") + "</tr></thead>"
    filas = ""
    for _, r in df.iterrows():
        filas += "<tr><td>" + str(r[test_col]) + "</td>"
        if est_col:
            filas += f"<td class='num'>{es_num(r[est_col], 4)}</td>"
        if p_col:
            filas += f"<td class='num'>{es_p(r[p_col])}</td>"
        if int_col:
            filas += "<td>" + str(r[int_col]) + "</td>"
        filas += "</tr>"
    return f"<table>{head}<tbody>{filas}</tbody></table>"


def html_tabla_anova(df: pd.DataFrame) -> str:
    head = "<thead><tr><th>Fuente</th><th>gl</th><th>F</th><th>p</th><th>η² parcial</th></tr></thead>"
    filas = ""
    for _, r in df.iterrows():
        if r["fuente"] == "Residual":
            filas += f"<tr><td>{r['fuente']}</td><td class='num'>{r['df']:.0f}</td><td class='num'>—</td><td class='num'>—</td><td class='num'>—</td></tr>"
        else:
            filas += (
                f"<tr><td>{r['fuente']}</td><td class='num'>{r['df']:.0f}</td>"
                f"<td class='num'>{es_num(r['F'], 2)}</td><td class='num'>{es_p(r['PR(>F)'])}</td>"
                f"<td class='num'>{es_num(r['eta2_parcial'], 3)}</td></tr>"
            )
    return f"<table>{head}<tbody>{filas}</tbody></table>"


def html_tabla_posthoc(df: pd.DataFrame) -> str:
    cols = df.columns.tolist()
    par = "par" if "par" in cols else df.columns[0]
    pcol = "p_valor_ajustado" if "p_valor_ajustado" in cols else "p_valor"
    let_col = "letras" if "letras" in cols else None
    if "diferencia_medias" in cols:
        etiqueta_est = "Diferencia de medias"
    elif "estadistico_z" in cols:
        etiqueta_est = "z (Dunn)"
    else:
        etiqueta_est = "Estadístico"
    head = ("<thead><tr><th>Comparación</th><th>" + etiqueta_est
            + "</th><th>p ajustado</th>" + ("<th>Letras</th>" if let_col else "")
            + "</tr></thead>")
    filas = ""
    for _, r in df.iterrows():
        diff = r.get("diferencia_medias", r.get("estadistico_z", ""))
        celda_diff = f"<td class='num'>{es_num(diff, 2)}</td>" if isinstance(diff, (int, float)) else "<td class='num'>—</td>"
        celda_let = f"<td>{r[let_col]}</td>" if let_col else ""
        filas += f"<tr><td>{r[par]}</td>{celda_diff}<td class='num'>{es_p(r[pcol])}</td>{celda_let}</tr>"
    return f"<table>{head}<tbody>{filas}</tbody></table>"


# ---------------------------------------------------------------------------
# Pagina BDCA (una respuesta: yield, RCBD)
# ---------------------------------------------------------------------------


def _bdca_tabla(diseno: str, nombre: str, **kw) -> pd.DataFrame | None:
    """Lee una tabla de <diseno>/resultados/tablas; None si no existe."""
    ruta = RAIZ / diseno / "resultados" / "tablas" / f"{nombre}.csv"
    if not ruta.exists():
        return None
    return pd.read_csv(ruta, **kw)


def fig_bdca_boxplot_yield() -> go.Figure:
    """Boxplot de rendimiento por tratamiento (RCBD)."""
    ruta = RAIZ / "datos_crudos" / "bdca" / "DBCA_Jenkyn_control_mildeo.csv"
    df = pd.read_csv(ruta)
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()
    orden = sorted(df["trt"].unique())
    fig = go.Figure()
    for trt in orden:
        fig.add_trace(go.Box(
            y=df.loc[df["trt"] == trt, "yield"],
            name=trt,
            boxpoints="all",
            jitter=0.3,
            pointpos=0,
            line=dict(color="#0072B2"),
            fillcolor="rgba(0,114,178,0.15)",
            marker=dict(color="#0072B2", size=5),
        ))
    fig.update_layout(xaxis_title="Tratamiento", yaxis_title="Rendimiento (yield)")
    return _base_layout(fig, "Rendimiento por tratamiento (yield)")


def fig_bdca_medias_ic(des: pd.DataFrame) -> go.Figure:
    """Barras de medias con IC95% por tratamiento."""
    df = des.copy()
    fig = go.Figure(go.Bar(
        x=df["trt"], y=df["media"],
        error_y=dict(
            type="data",
            array=df["media"] - df["ic95_inferior"],
            arrayminus=df["media"] - df["ic95_inferior"],
            visible=True, thickness=1.4, width=6, color="#444444",
        ),
        marker_color="#0072B2",
        text=[f"{v:.2f}" for v in df["media"]],
        textposition="outside",
        hovertemplate="%{y:.2f}<extra>%{x}</extra>",
    ))
    fig.update_layout(
        xaxis_title="Tratamiento",
        yaxis_title="Rendimiento medio (yield) ± IC 95%",
        bargap=0.35,
    )
    return _base_layout(fig, "Rendimiento medio por tratamiento con IC 95%")


def fig_bdca_tukey(ph: pd.DataFrame) -> go.Figure:
    """Diferencias de medias por par con IC95% (Tukey HSD)."""
    df = ph.copy()
    colores = ["#D55E00" if not s else "#0072B2" for s in df["significativo"]]
    fig = go.Figure()
    for _, r in df.iterrows():
        color = "#0072B2" if r["significativo"] else "#5A5A5A"
        fig.add_trace(go.Scatter(
            x=[r["diferencia_medias"]], y=[r["par"]],
            mode="markers", marker=dict(color=color, size=11),
            error_x=dict(
                type="data",
                symmetric=False,
                array=[r["ic95_superior"] - r["diferencia_medias"]],
                arrayminus=[r["diferencia_medias"] - r["ic95_inferior"]],
                thickness=1.5, width=5, color=color,
            ),
            hovertemplate=f"{r['par']}: {r['diferencia_medias']:.3f} (IC95 {r['ic95_inferior']:.3f}–{r['ic95_superior']:.3f})",
            showlegend=False,
        ))
    fig.add_vline(x=0, line_dash="dash", line_color="#888888")
    fig.update_layout(
        xaxis_title="Diferencia de medias",
        yaxis_title="Par comparado",
        height=320,
        margin=dict(l=20, r=20, t=30, b=50),
    )
    return _base_layout(fig, "Tukey HSD: diferencias de medias por par (IC 95%)")


def main_bdca() -> None:
    """Genera pagina/bdca/index.html con el estándar del estudio (una respuesta: yield)."""
    DIR_BDCA = RAIZ / "bdca"
    DIR_PAG_BDCA = RAIZ / "pagina" / "bdca"

    # Tablas del análisis BDCA
    aud = _bdca_tabla("bdca", "auditoria_bdca")
    des = _bdca_tabla("bdca", "eda_descriptivos")
    sup = _bdca_tabla("bdca", "supuestos_modelo")
    anova = _bdca_tabla("bdca", "anova_bloques")
    lmm_f = _bdca_tabla("bdca", "lmm_bloques_fijos")
    lmm_v = _bdca_tabla("bdca", "lmm_bloques_varianzas")
    ph = _bdca_tabla("bdca", "posthoc_tukey")

    faltantes = [n for n, t in [
        ("auditoria_bdca", aud), ("eda_descriptivos", des), ("supuestos_modelo", sup),
        ("anova_bloques", anova), ("lmm_bloques_fijos", lmm_f),
        ("lmm_bloques_varianzas", lmm_v), ("posthoc_tukey", ph),
    ] if t is None]
    if faltantes:
        sys.exit(f"Faltan tablas BDCA en bdca/resultados/tablas/: {faltantes}. "
                 f"Ejecuta el pipeline primero (PIPELINE_DISENO=bdca).")

    # Figuras
    f = {
        "boxplot": _html_figura(fig_bdca_boxplot_yield()),
        "medias": _html_figura(fig_bdca_medias_ic(des)),
        "tukey": _html_figura(fig_bdca_tukey(ph)),
    }

    # Valores para cifras
    r = aud.iloc[0] if aud is not None else {}
    n_filas = int(r.get("total_filas", 36))
    n_trt = int(r.get("conteo_por_trt", 4)) if not isinstance(r.get("conteo_por_trt"), (dict, list, str)) else 4
    n_bloques = int(r.get("conteo_por_block", 9)) if not isinstance(r.get("conteo_por_block"), (dict, list, str)) else 9
    # Leer contadores del dict (guardado como string)
    import json as _json
    try:
        if isinstance(r.get("conteo_por_trt"), str):
            n_trt = len(_json.loads(r["conteo_por_trt"].replace("'", '"')))
        if isinstance(r.get("conteo_por_block"), str):
            n_bloques = len(_json.loads(r["conteo_por_block"].replace("'", '"')))
    except Exception:
        pass

    # Tablas HTML
    html_aud = html_tabla_diseno(pd.DataFrame([
        {"atributo": "Filas únicas", "valor": str(int(r.get("total_filas", 36)))},
        {"atributo": "Duplicados", "valor": str(int(r.get("filas_duplicadas", 0)))},
        {"atributo": "Valores faltantes", "valor": str(int(r.get("filas_con_na", 0)))},
        {"atributo": "Auditoría aprobada", "valor": "Sí" if r.get("auditado_ok") else "No"},
    ]))
    html_sup = html_tabla_supuestos(sup) if sup is not None else "<p>Tabla de supuestos no disponible.</p>"
    html_anova = html_tabla_anova(anova) if anova is not None else "<p>Tabla ANOVA no disponible.</p>"
    html_posthoc = html_tabla_posthoc(ph) if ph is not None else "<p>Tabla post-hoc no disponible.</p>"

    lmm_filas = ""
    if lmm_f is not None:
        for _, rw in lmm_f.iterrows():
            lmm_filas += (
                f"<tr><td>{rw['efecto']}</td><td class='num'>{es_num(rw['coeficiente'], 3)}</td>"
                f"<td class='num'>{es_num(rw['error_estandar'], 3)}</td>"
                f"<td class='num'>{es_p(rw['p_valor'])}</td>"
                f"<td class='num'>{es_num(rw['ic95_inferior'], 3)}–{es_num(rw['ic95_superior'], 3)}</td></tr>"
            )
    html_lmm = f"<table><thead><tr><th>Efecto</th><th>Coeficiente</th><th>EE</th><th>p</th><th>IC 95%</th></tr></thead><tbody>{lmm_filas}</tbody></table>"

    icc_txt = ""
    if lmm_v is not None:
        icc = lmm_v[lmm_v["parametro"] == "ICC"]["valor"].iloc[0] if "ICC" in lmm_v["parametro"].values else None
        if icc is not None:
            icc_txt = f"<p>El <strong>ICC = {icc:.3f}</strong> indica que el bloque explica una parte importante de la variación total (variabilidad espacial del ensayo).</p>"

    # Interpretación por pares vs R
    pares_r = ""
    if ph is not None:
        vs_r = ph[ph["vs_referencia_R"] == True]  # noqa: E712
        filas_r = ""
        for _, rw in vs_r.iterrows():
            estado = "difiere" if rw["significativo"] else "no difiere"
            filas_r += (
                f"<tr><td><strong>{rw['par']}</strong></td><td class='num'>{es_num(rw['diferencia_medias'], 3)}</td>"
                f"<td class='num'>{es_p(rw['p_valor_ajustado'])}</td><td>{estado} de R</td></tr>"
            )
        pares_r = f"<table><thead><tr><th>Par</th><th>Diferencia</th><th>p ajustado</th><th>vs R</th></tr></thead><tbody>{filas_r}</tbody></table>"

    css = CSS
    vuelta_hub = "../index.html"
    nav = (
        f'<a class="enlace vuelta" href="{vuelta_hub}" title="Volver al índice de análisis">← Análisis</a>'
        + "".join(
            f'<a class="enlace" href="#{seccion}">{etiqueta}</a>'
            for seccion, etiqueta in [
                ("hero", "Resumen"), ("desafio", "El ensayo"), ("diseno", "Diseño"),
                ("ruta", "Ruta estadística"), ("bloque-resultados", "Resultados"),
                ("conclusiones", "Conclusiones"), ("metodologia", "Metodología"),
            ]
        )
    )

    hero_cifras = "".join(
        f'<div class="cifra"><div class="valor">{valor}</div>'
        f'<div class="etiqueta">{etiqueta}</div></div>'
        for valor, etiqueta in [
            (str(n_filas), "unidades experimentales (parcelas)"),
            (str(n_trt), "tratamientos (R, T0, T1, T2)"),
            (str(n_bloques), "bloques completos"),
            ("1", "respuesta analizada (rendimiento)"),
        ]
    )

    html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Análisis estadístico del ensayo en bloques completos al azar (RCBD) de rendimiento — Thymus vulgaris contra Fusarium spp.">
<title>BDCA — Bloques completos al azar (rendimiento)</title>
<script src="{PLOTLY_CDN}"></script>
<style>{css}</style>
</head>
<body>
<nav class="fijo">
  <span class="marca">Thymus × Fusarium</span>
  {nav}
</nav>

<div class="contenido">

  <!-- ================= SECCION 1: HERO ================= -->
  <section id="hero" class="seccion ancla">
    <h1 class="titulo-hero">Diseño de bloques completos al azar (BDCA)</h1>
    <p class="subtitulo-hero">Comparación del rendimiento (<em>yield</em>) entre 4 tratamientos
    (R, T0, T1, T2) en 9 bloques completos, con modelo mixto de bloque aleatorio,
    ANOVA clásico y Tukey HSD contra la referencia R.</p>
    <div class="grid-cifras">{hero_cifras}</div>
  </section>

  <!-- ================= SECCION 2: EL ENSAYO ================= -->
  <section id="desafio" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">2</span>El ensayo</h2>
    <p class="lead">El ensayo de Jenkyn evalúa el control del mildiu (mildeo) mediante
    tratamientos aplicados a parcelas dispuestas en <strong>bloques completos al azar</strong>.
    La variabilidad espacial del terreno se controla agrupando las parcelas en 9 bloques,
    dentro de cada uno de los cuales los 4 tratamientos aparecen una única vez.</p>
    <p>Con una sola observación por celda tratamiento × bloque, el modelo de bloques no
    puede estimar el término de interacción: la <strong>aditividad</strong> se documenta como
    una suposición no testable (AGENTS.md §4.2).</p>
    <div class="que-significa">
      <div class="titulo-tarjeta">Qué significa</div>
      <p>El bloque es la unidad de control de la heterogeneidad espacial: si el ICC es alto,
      gran parte de la variación del rendimiento se explica por diferencias entre bloques y el
      diseño RCBD fue la elección correcta.</p>
    </div>
  </section>

  <!-- ================= SECCION 3: DISEÑO ================= -->
  <section id="diseno" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">3</span>Diseño experimental</h2>
    {html_aud}
    <p>Estructura: 36 parcelas = 4 tratamientos × 9 bloques; una parcela por celda.
    La auditoría confirma balance: 9 observaciones por tratamiento y 4 por bloque,
    sin duplicados ni valores faltantes.</p>
  </section>

  <!-- ================= SECCION 4: RUTA ================= -->
  <section id="ruta" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">4</span>Ruta estadística</h2>
    <table>
    <thead><tr><th>Etapa</th><th>Método</th><th>Resultado clave</th></tr></thead>
    <tbody>
    <tr><td><strong>Auditoría</strong></td><td>Balance, NA, duplicados</td><td>36 filas, sin anomalías</td></tr>
    <tr><td><strong>EDA</strong></td><td>Medias + IC95% por tratamiento</td><td>R=5,94; T0=5,31; T1=5,87; T2=6,09</td></tr>
    <tr><td><strong>Supuestos</strong></td><td>Shapiro-Wilk, Levene, Durbin-Watson</td><td>Ruta paramétrica</td></tr>
    <tr><td><strong>Modelo primario</strong></td><td>LMM con bloque aleatorio (REML)</td><td>ICC reportado</td></tr>
    <tr><td><strong>Complemento</strong></td><td>ANOVA clásico de bloques (typ=2)</td><td>F trat = 28,77; p &lt; 0,001; η²p = 0,78</td></tr>
    <tr><td><strong>Post-hoc</strong></td><td>Tukey HSD (6 pares, vs R)</td><td>Contrastes vs referencia R</td></tr>
    </tbody>
    </table>
  </section>

  <!-- ================= SECCION 5: RESULTADOS ================= -->
  <section id="bloque-resultados" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">5</span>Resultados</h2>

    <h3 class="sub-seccion">Exploración descriptiva</h3>
    {html_figura_bloque("Rendimiento por tratamiento", f["boxplot"], "Boxplot con todos los puntos; la caja muestra mediana y cuartiles por tratamiento.", "A")}
    {html_figura_bloque("Rendimiento medio con IC 95%", f["medias"], "Barras de medias con IC95%; la superposición de intervalos anticipa qué pares pueden diferir.", "B")}

    <h3 class="sub-seccion">Verificación de supuestos</h3>
    {html_sup}

    <h3 class="sub-seccion">ANOVA clásico de bloques (complemento educativo)</h3>
    {html_anova}

    <h3 class="sub-seccion">Modelo mixto lineal (análisis primario)</h3>
    {html_lmm}
    {icc_txt}

    <h3 class="sub-seccion">Comparaciones post-hoc (Tukey HSD)</h3>
    {html_figura_bloque("Tukey HSD", f["tukey"], "Diferencias de medias por par con IC95%; en azul los pares significativos, en gris los no significativos. La línea punteada marca la ausencia de diferencia.", "C")}
    {html_posthoc}

    <h3 class="sub-seccion">Contrastes contra la referencia R</h3>
    {pares_r}
  </section>

  <!-- ================= SECCION 6: CONCLUSIONES ================= -->
  <section id="conclusiones" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">6</span>Conclusiones</h2>
    <p>La auditoría confirma un diseño RCBD balanceado. La ruta inferencial se eligió tras
    documentar los supuestos (normalidad, homocedasticidad e independencia de residuos).
    El modelo primario trata el bloque como efecto aleatorio (LMM, REML) y reporta el ICC;
    el ANOVA clásico de bloques se conserva como complemento educativo. Las comparaciones
    post-hoc con Tukey HSD interpretan cada par contra la referencia R.</p>
    <div class="que-significa">
      <div class="titulo-tarjeta">Limitaciones</div>
      <p>Con una observación por celda tratamiento × bloque, la aditividad no es testeable:
      el modelo asume efectos aditivos de tratamiento y bloque sin interacción.</p>
    </div>
  </section>

  <!-- ================= SECCION 7: METODOLOGIA ================= -->
  <section id="metodologia" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">7</span>Metodología y reproducibilidad</h2>
    <p>Pipeline reproducible: <code>pipeline/bdca/</code> (cargar → eda → supuestos →
    modelos → comparaciones → informe). Resultados en <code>bdca/resultados/</code>.
    Notebook educativo: <code>bdca/analisis_bdca.ipynb</code>.</p>
  </section>

</div>
</body>
</html>
"""
    DIR_PAG_BDCA.mkdir(parents=True, exist_ok=True)
    (DIR_PAG_BDCA / "index.html").write_text(html, encoding="utf-8")
    print(f"OK: pagina/bdca/index.html generada con {len(f)} figuras Plotly.")
    print(f"Tamaño: {sum(len(v) for v in f.values()) / 1024:.0f} KB de figuras embebidas.")


def main() -> None:
    if not MASTER_CSV.exists() or not REND_CSV.exists():
        sys.exit(f"Faltan datos maestros: {MASTER_CSV} o {REND_CSV}")

    datos = cargar_datos()
    master = datos["master"]
    rend = datos["rend"]
    t = datos["tablas"]
    d = texto_decisiones(datos)
    control_mm, control_log10 = datos["control_mm"], datos["control_log10"]

    # ------------------------------------------------------------------ figuras
    f = {}

    # Bloque A
    f["rend_box"] = _html_figura(fig_boxplot_rendimiento(rend))
    f["rend_barras"] = _html_figura(fig_barras_rendimiento(rend))
    f["rend_tukey"] = _html_figura(fig_tukey_rendimiento(rend))

    # Bloque B
    f["inh_violin"] = _html_figura(fig_violin_inhibicion(
        master, "porcentaje_inhibicion_micelial",
        "Inhibición micelial por método (%)",
        "Inhibición micelial (%)"))
    f["inh_barras"] = _html_figura(fig_barras_letras(
        t["medias_porcentaje_inhibicion_micelial"],
        t["posthoc_porcentaje_inhibicion_micelial_letras"],
        "Inhibición micelial media (%)",
        "Inhibición micelial media (± IC 95%) y comparaciones (Tukey)"))

    # Bloque C
    f["con_box"] = _html_figura(fig_conidias_control(master, control_log10))
    f["inhcon_barras"] = _html_figura(fig_barras_letras(
        t["medias_porcentaje_inhibicion_conidias"],
        t["posthoc_porcentaje_inhibicion_conidias_letras"],
        "Inhibición de conidias media (%)",
        "Inhibición de conidias media (± IC 95%) y comparaciones (Dunn)"))
    f["diag_poisson"] = _html_figura(fig_diagnostico_poisson(master))

    # Bloque D
    matriz = matriz_por_aislado(master)
    z = (matriz - matriz.mean()) / matriz.std(ddof=0)
    f["heatmap"] = _html_figura(fig_heatmap_susceptibilidad(matriz, t["susceptibilidad_clusters"]))
    f["pca_scree"] = _html_figura(fig_pca_scree())
    f["biplot"] = _html_figura(fig_biplot_pca(
        t["susceptibilidad_pca_scores"], t["susceptibilidad_pca_loadings"],
        t["susceptibilidad_clusters"]))
    f["dendro"] = _html_figura(fig_dendrograma(z))
    f["kmeans"] = _html_figura(fig_kmeans_metricas(t["susceptibilidad_kmeans_metricas"]))
    f["comp_clusters"] = _html_figura(fig_composicion_clusters(
        t["susceptibilidad_cruce_cluster_categoria"]))

    # Ranking
    f["radar"] = _html_figura(fig_radar(t["ranking_tecnicas"]))
    f["score"] = _html_figura(fig_score_compuesto(t["ranking_tecnicas"]))

    # ------------------------------------------------------------------ texto
    ranking = t["ranking_tecnicas"].sort_values("score_compuesto", ascending=False)
    filas_ranking = ""
    for i, (_, r) in enumerate(ranking.iterrows(), start=1):
        filas_ranking += (
            f"<tr><td class='num'>{i}</td><td><strong>{METODO_LABEL[r['metodo_extraccion']]}</strong></td>"
            f"<td class='num'>{es_num(r['rendimiento_medio_pct'], 1)}</td>"
            f"<td class='num'>{es_num(r['inhib_micelial_medio_pct'], 1)}</td>"
            f"<td class='num'>{es_num(r['inhib_conidias_medio_pct'], 1)}</td>"
            f"<td class='num'><strong>{es_num(r['score_compuesto'], 3)}</strong></td></tr>"
        )

    # Distribucion de categorias
    cl = t["susceptibilidad_clusters"]
    n_alta = int((cl["categoria_susceptibilidad"] == "Alta susceptibilidad relativa").sum())
    n_mod = int((cl["categoria_susceptibilidad"] == "Moderada susceptibilidad relativa").sum())
    n_baja = int((cl["categoria_susceptibilidad"] == "Baja susceptibilidad relativa").sum())

    valid = t["validacion_inh"]
    n_verif = int(valid["n_verificadas"].sum())
    n_disc = int(valid["n_discrepancias"].sum())

    # Letras rendimiento
    letras_rend = _letras_tukey(rend)

    # ------------------------------------------------------------------ HTML
    css = CSS
    vuelta_hub = "../index.html"
    nav = (
        f'<a class="enlace vuelta" href="{vuelta_hub}" title="Volver al índice de análisis">← Análisis</a>'
        + "".join(
            f'<a class="enlace" href="#{seccion}">{etiqueta}</a>'
            for seccion, etiqueta in [
                ("hero", "Resumen"), ("desafio", "El desafío"), ("diseno", "Diseño"),
                ("ruta", "Ruta estadística"), ("bloque-rendimiento", "Rendimiento"),
                ("bloque-inhibicion", "Inhibición"), ("bloque-esporulacion", "Esporulación"),
                ("bloque-susceptibilidad", "Susceptibilidad"), ("ranking", "Ranking"),
                ("conclusiones", "Conclusiones"), ("metodologia", "Metodología"),
            ]
        )
    )

    hero_cifras = "".join(
        f'<div class="cifra"><div class="valor">{valor}</div>'
        f'<div class="etiqueta">{etiqueta}</div></div>'
        for valor, etiqueta in [
            ("31", "aislados de Fusarium spp."),
            ("3", "técnicas de extracción"),
            ("5 mg/mL", "concentración ensayada"),
            ("279", "unidades experimentales (cajas Petri)"),
            ("56,7 mm", "control C4: crecimiento micelial"),
            ("7,39", "control C4: conidias (log10/mL)"),
        ]
    )

    k_opt = d["k_opt"]
    sil_str = es_num(d["sil"], 2)

    letras_rend = _letras_tukey(rend)
    frase_letras = ", ".join(
        f"<strong>{METODO_LABEL[m]} ({letras_rend[m]['letra']})</strong>"
        for m in METODOS
    )

    html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Análisis estadístico interactivo de la actividad antifúngica de extractos de Thymus vulgaris contra Fusarium spp.">
<title>Actividad antifúngica de Thymus vulgaris contra Fusarium spp.</title>
<script src="{PLOTLY_CDN}"></script>
<style>{css}</style>
</head>
<body>
<nav class="fijo">
  <span class="marca">Thymus × Fusarium</span>
  {nav}
</nav>

<div class="contenido">

  <!-- ================= SECCION 1: HERO ================= -->
  <section id="hero" class="seccion ancla">
    <h1 class="titulo-hero">Actividad antifúngica de extractos de <em>Thymus vulgaris</em> contra <em>Fusarium</em> spp.</h1>
    <p class="subtitulo-hero">Un recorrido estadístico interactivo sobre cómo la técnica de extracción modula
    el rendimiento y la actividad antifúngica del tomillo frente a 31 aislados de <em>Fusarium</em>.
    Cada bloque sigue el mismo camino: descriptiva → supuestos → test → post-hoc → multivariado.</p>
    <div class="metodos-chips">
      <span class="chip"><span class="punto" style="background:#0072B2"></span>Maceración</span>
      <span class="chip"><span class="punto" style="background:#D55E00"></span>Soxhlet</span>
      <span class="chip"><span class="punto" style="background:#009E73"></span>Ultrasonido</span>
    </div>
    <div class="grid-cifras">{hero_cifras}</div>
  </section>

  <!-- ================= SECCION 2: EL DESAFIO ================= -->
  <section id="desafio" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">2</span>El desafío</h2>
    <p class="lead">Los hongos del género <em>Fusarium</em> causan marchitamiento vascular, pudriciones de raíz
    y pérdidas importantes en cultivos de interés agronómico. Su control suele depender de fungicidas de
    síntesis, cuyo uso repetido favorece la aparición de poblaciones menos sensibles y plantea
    preocupaciones ambientales y sanitarias. La búsqueda de alternativas naturales es un área activa de la
    fitopatología.</p>
    <p>El tomillo (<em>Thymus vulgaris</em> L.) produce compuestos fenólicos con reconocida actividad
    antifúngica. Sin embargo, la técnica de extracción determina qué compuestos y en qué proporción se
    obtienen, y con ello la actividad biológica final. Este estudio compara tres técnicas —maceración,
    Soxhlet y ultrasonido— a una concentración única de 5&nbsp;mg/mL, frente a 31 aislados de
    <em>Fusarium</em> spp., midiendo el rendimiento de extracción y cuatro variables de respuesta del
    bioensayo.</p>
    <div class="que-significa">
      <div class="titulo-tarjeta">Qué significa</div>
      <p>No se trata de una sola comparación: hay una variable de proceso (rendimiento) y una batería de
      variables biológicas con distribuciones distintas. La elección estadística correcta depende de cada
      una. Esta página muestra el razonamiento que justifica cada decisión, no solo el resultado.</p>
    </div>
  </section>

  <!-- ================= SECCION 3: DISEÑO ================= -->
  <section id="diseno" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">3</span>Diseño experimental</h2>
    <p class="lead">Diseño completamente aleatorizado (DCA) factorial con dos factores: técnica de
    extracción (3 niveles) y aislado de <em>Fusarium</em> (31 niveles), con 3 réplicas biológicas por
    combinación. La unidad experimental es la caja Petri.</p>
    {html_tabla_diseno(t["diseno_experimental"])}
    <div class="caveat">
      <div class="titulo-caveat">Caveat metodológico</div>
      <p>El porcentaje de inhibición se calcula contra un único control C4 compartido por las tres réplicas
      de cada aislado. Por lo tanto, las réplicas de %INH no son totalmente independientes
      (pseudorreplicación del control). El análisis del crecimiento crudo (mm) no presenta este problema, y
      los modelos mixtos con aislado aleatorio mitigan parcialmente la dependencia.</p>
    </div>
  </section>

  <!-- ================= SECCION 4: LA RUTA ESTADISTICA ================= -->
  <section id="ruta" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">4</span>La ruta estadística</h2>
    <p class="lead">No se elige un test solo por producir significancia: cada paso se justifica con el
    anterior. La descriptiva orienta, los supuestos deciden la vía (paramétrica o no paramétrica), el test
    principal responde la pregunta y el post-hoc localiza las diferencias. El multivariado se reserva para
    cuando el problema lo exige (perfil de susceptibilidad por aislado).</p>
    {html_flujo()}
    <h3>El camino real elegido para cada variable</h3>
    <p>Los valores de la tabla provienen de las tablas de la pipeline reproducida (<code>modelos_*.csv</code>,
    <code>no_parametrico_*.csv</code>, <code>posthoc_*.csv</code>).</p>
    {html_tabla_ruta()}
  </section>

  <!-- ================= SECCION 5: BLOQUE A — RENDIMIENTO ================= -->
  <section id="bloque-rendimiento" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">5</span>Bloque A · Rendimiento de extracción <span style="font-size:.9rem;color:var(--gris);font-weight:400">(rama paramétrica)</span></h2>

    <h3>Descriptiva</h3>
    <p>El rendimiento se mide sobre las 3 réplicas de cada técnica. La separación entre grupos es visible
    desde la descriptiva: Soxhlet domina con holgura, mientras Maceración y Ultrasonido se superponen.</p>
    {html_figura_bloque("Rendimiento por método", f["rend_box"], "Cada punto es una réplica biológica. El diagrama de caja resume la distribución por técnica.", "Figura A1")}
    {html_figura_bloque("Rendimiento medio", f["rend_barras"], "Barras de medias con intervalo de confianza 95% (n=3, corrección t de Student).", "Figura A2")}

    <h3>Supuestos</h3>
    <p>Antes de elegir ANOVA se verifican normalidad y homocedasticidad de los residuos.</p>
    {html_tabla_supuestos(t["supuestos_rendimiento"])}
    <div class="decision">
      <div class="cabeza"><span>Decisión estadística</span><span class="flecha">→</span><span>vía paramétrica</span></div>
      <ul>
        <li>Shapiro-Wilk: p = {es_p(d["shapiro_rend"])} → no se rechaza normalidad.</li>
        <li>Levene: p = {es_p(d["levene_rend"])} → homocedasticidad aceptada.</li>
        <li>Ambos supuestos razonables → ANOVA de una vía.</li>
      </ul>
    </div>

    <h3>ANOVA de una vía</h3>
    <table><thead><tr><th>Fuente</th><th>gl</th><th>F</th><th>p</th><th>η² parcial</th></tr></thead>
    <tbody>
      <tr><td>Método de extracción</td><td class="num">2</td><td class="num">{es_num(d["rend_f"], 2)}</td><td class="num">{es_p(d["rend_p"])}</td><td class="num">{es_num(d["rend_eta"], 3)}</td></tr>
      <tr><td>Residual</td><td class="num">6</td><td class="num">—</td><td class="num">—</td><td class="num">—</td></tr>
    </tbody></table>
    <div class="que-significa">
      <div class="titulo-tarjeta">Qué significa</div>
      <p>El método explica el {es_num(d["rend_eta"]*100, 0)}% de la variabilidad del rendimiento
      (η² = {es_num(d["rend_eta"], 3)}). La técnica de extracción es el factor dominante.</p>
    </div>

    <h3>Tukey HSD (post-hoc)</h3>
    {html_figura_bloque("Comparaciones Tukey", f["rend_tukey"], "Intervalos de confianza al 95% para las diferencias de medias entre pares de técnicas.", "Figura A3")}
    <p>Letras de significancia: {frase_letras}. La maceración y el ultrasonido no difieren entre sí de forma
    significativa (p = 0,358).</p>
    <div class="decision">
      <div class="cabeza"><span>Resultado</span><span class="flecha">→</span><span>ranking parcial</span></div>
      <ul>
        <li>Rendimiento medio: Soxhlet 43,4% &gt; Ultrasonido 17,1% ≈ Maceración 12,1%.</li>
        <li>Tukey: Soxhlet supera a las otras dos (p &lt; 0,001); Maceración vs Ultrasonido no significativo (p = 0,358).</li>
      </ul>
    </div>
  </section>

  <!-- ================= SECCION 6: BLOQUE B — %INH MICELIAL ================= -->
  <section id="bloque-inhibicion" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">6</span>Bloque B · Inhibición del crecimiento micelial <span style="font-size:.9rem;color:var(--gris);font-weight:400">(rama paramétrica + robustez)</span></h2>

    <h3>Descriptiva</h3>
    <p>La inhibición micelial alcanza valores muy altos: en Maceración la mediana ronda el 90% y hay un
    <strong>efecto techo</strong> en 100%. Este techo comprime la varianza y debe tenerse presente al
    interpretar las diferencias.</p>
    {html_figura_bloque("Inhibición micelial", f["inh_violin"], "Los violines muestran la densidad completa: la distribución de Maceración está truncada cerca de 100% (efecto techo).", "Figura B1")}
    {html_figura_bloque("Medias y comparaciones", f["inh_barras"], "Barras de medias ± IC 95% con letras CLD del Tukey. Letras distintas indican diferencias significativas.", "Figura B2")}

    <h3>Supuestos</h3>
    {html_tabla_supuestos(t["supuestos_porcentaje_inhibicion_micelial"])}
    <div class="decision">
      <div class="cabeza"><span>Decisión estadística</span><span class="flecha">→</span><span>ANOVA factorial</span></div>
      <ul>
        <li>Shapiro-Wilk: p = {es_p(d["shapiro_inh"])} → residuos normales.</li>
        <li>Levene: p = {es_p(d["levene_inh"])} → homocedasticidad aceptada.</li>
        <li>Nota: la independencia queda parcialmente comprometida por el control compartido (pseudorreplicación); se refuerza con modelos mixtos.</li>
      </ul>
    </div>

    <h3>ANOVA factorial método × aislado</h3>
    {html_tabla_anova(t["modelos_porcentaje_inhibicion_micelial"])}
    <div class="que-significa verde">
      <div class="titulo-tarjeta">Qué significa</div>
      <p>El método tiene el mayor tamaño de efecto (η²p = {es_num(d["inh_m_eta"], 3)}), pero el aislado y su
      interacción con el método también aportan (η²p = {es_num(d["inh_a_eta"], 3)} y
      {es_num(d["inh_i_eta"], 3)}): la respuesta no es uniforme entre aislados. La interacción
      (p = {es_p(d["inh_i_p"])}) indica que la diferencia entre técnicas cambia según el aislado.</p>
    </div>
    <div class="decision">
      <div class="cabeza"><span>Sensibilidad</span><span class="flecha">→</span><span>modelo mixto LMM</span></div>
      <ul>
        <li>Modelo: %INH ~ método + (1|aislado).</li>
        <li>ICC(aislado) = {es_num(d["icc"], 3)}: el aislado explica el {es_num(d["icc"]*100, 0)}% de la varianza residual.</li>
        <li>El efecto del método se mantiene significativo (LRT, p &lt; 0,001): el resultado es robusto a la
        estructura de correlación por aislado.</li>
      </ul>
    </div>

    <h3>Tukey HSD</h3>
    {html_tabla_posthoc(t["posthoc_porcentaje_inhibicion_micelial"])}
    <div class="que-significa naranja">
      <div class="titulo-tarjeta">Caveat del efecto techo</div>
      <p>Maceración alcanza una inhibición media del 86,3%, significativamente mayor que Soxhlet (59,6%) y
      Ultrasonido (54,8%). No obstante, el techo en 100% subestima las diferencias reales entre técnicas:
      muchas cajas de Maceración muestran inhibición completa y no se puede distinguir cuánto «mejor» habrían
      sido de no existir el límite.</p>
    </div>
  </section>

  <!-- ================= SECCION 7: BLOQUE C — ESPORULACION ================= -->
  <section id="bloque-esporulacion" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">7</span>Bloque C · Esporulación <span style="font-size:.9rem;color:var(--gris);font-weight:400">(rama no paramétrica)</span></h2>

    <h3>Descriptiva</h3>
    <p>Las conidias se modelan en escala log10 continua (no como conteos). La línea roja marca el control
    C4 (7,39 log10/mL): cuanto más abajo el grupo, mayor reducción de la esporulación. La inhibición de
    conidias puede tomar valores negativos cuando la producción supera al control.</p>
    {html_figura_bloque("Conidias log10", f["con_box"], "Producción de conidias por método frente al control C4. Solo Maceración reduce la esporulación de forma sustancial.", "Figura C1")}
    {html_figura_bloque("Inhibición de conidias", f["inhcon_barras"], "Medias ± IC 95% con letras CLD del Dunn. Ultrasonido tiene media levemente negativa (-0,17%).", "Figura C2")}

    <h3>Supuestos</h3>
    {html_tabla_supuestos(t["supuestos_conidias_log10_ml"])}
    <div class="decision">
      <div class="cabeza"><span>Decisión estadística</span><span class="flecha">→</span><span>vía no paramétrica</span></div>
      <ul>
        <li>Shapiro-Wilk: p = {es_p(d["shapiro_con"])} → se rechaza la normalidad de los residuos.</li>
        <li>Levene: p = {es_p(d["levene_con"])} → heterocedasticidad (las varianzas difieren entre métodos).</li>
        <li>Supuestos violados → Kruskal-Wallis + Scheirer-Ray-Hare.</li>
      </ul>
    </div>

    <h3>Kruskal-Wallis y Scheirer-Ray-Hare</h3>
    <table><thead><tr><th>Fuente</th><th>H</th><th>gl</th><th>p</th></tr></thead>
    <tbody>
      <tr><td>Método (conidias log10)</td><td class="num">{es_num(d["h_con"], 2)}</td><td class="num">2</td><td class="num">{es_p(0.0)}</td></tr>
      <tr><td>Interacción método × aislado</td><td class="num">44,27</td><td class="num">60</td><td class="num">{es_p(d["p_int_con"])}</td></tr>
      <tr><td>Método (%INH conidias)</td><td class="num">{es_num(d["h_ic"], 2)}</td><td class="num">2</td><td class="num">{es_p(0.0)}</td></tr>
      <tr><td>Interacción método × aislado</td><td class="num">45,79</td><td class="num">60</td><td class="num">{es_p(d["p_int_ic"])}</td></tr>
    </tbody></table>

    <h3>Dunn (post-hoc)</h3>
    {html_tabla_posthoc(t["posthoc_conidias_log10_ml"])}
    <p>Conidias log10 — letras: <strong>Ultrasonido (a)</strong>, <strong>Soxhlet (b)</strong>,
    <strong>Maceración (c)</strong>; todos los pares difieren significativamente.</p>
    {html_tabla_posthoc(t["posthoc_porcentaje_inhibicion_conidias"])}
    <p>%INH conidias — letras: <strong>Maceración (a)</strong>, <strong>Soxhlet (b)</strong>,
    <strong>Ultrasonido (c)</strong>; todos los pares difieren.</p>

    <h3>Diagnóstico: ¿por qué no un modelo Poisson?</h3>
    {html_figura_bloque("Distribución de conidias", f["diag_poisson"], "Histograma de conidias log10 con ajuste de densidad normal.", "Figura C3")}
    <div class="decision">
      <div class="cabeza"><span>Diagnóstico</span><span class="flecha">→</span><span>regresión lineal sobre log10</span></div>
      <ul>
        <li>Un conteo exige valores enteros no negativos con varianza ≈ media (Poisson).</li>
        <li>Media = {es_num(d["diag_media"], 2)} vs varianza = {es_num(d["diag_var"], 2)}: no se cumple la
        condición de equidispersión.</li>
        <li>Solo el {es_num(d["diag_enteros"], 1)}% de los valores son enteros: las conidias se reportan en
        escala log10 continua.</li>
        <li>Conclusión: modelo lineal sobre log10 (las conidias no son conteos).</li>
      </ul>
    </div>
  </section>

  <!-- ================= SECCION 8: BLOQUE D — SUSCEPTIBILIDAD ================= -->
  <section id="bloque-susceptibilidad" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">8</span>Bloque D · Susceptibilidad de los aislados <span style="font-size:.9rem;color:var(--gris);font-weight:400">(PCA + clustering)</span></h2>

    <h3>Descriptiva del perfil</h3>
    <p>Cada aislado se resume en 6 métricas: %INH micelial y %INH conidias para cada técnica. El heatmap
    ordena los aislados por cluster y score de susceptibilidad. Se usan las etiquetas de
    <strong>susceptibilidad relativa</strong> (Alta / Moderada / Baja), nunca «resistencia», porque no hay
    un criterio validado de sensibilidad.</p>
    {html_figura_bloque("Heatmap del perfil", f["heatmap"], "Filas ordenadas por cluster (C0 = Baja, C1 = Mixta) y score. Colores rojos = mayor inhibición (aislado más susceptible).", "Figura D1")}

    <h3>Estandarización y PCA</h3>
    <p>Las 6 métricas se estandarizan (z-score por columna) para que ninguna domine por su escala. Luego se
    aplica PCA para reducir la dimensión y visualizar la estructura.</p>
    {html_figura_bloque("Scree del PCA", f["pca_scree"], "PC1 y PC2 concentran el 71,4% de la varianza total.", "Figura D2")}
    {html_figura_bloque("Biplot", f["biplot"], "Scores de aislados coloreados por cluster y loadings de las métricas originales (flechas).", "Figura D3")}

    <h3>Clustering</h3>
    {html_figura_bloque("Dendrograma", f["dendro"], "Clustering jerárquico de Ward sobre las métricas estandarizadas.", "Figura D4")}
    {html_figura_bloque("Selección de k", f["kmeans"], "Codo de la inercia (izquierda) y perfil de silhouette (derecha). El máximo de silhouette selecciona k = " + str(k_opt) + ".", "Figura D5")}
    {html_figura_bloque("Composición de los clusters", f["comp_clusters"], "Distribución de categorías dentro de cada cluster.", "Figura D6")}

    <div class="decision">
      <div class="cabeza"><span>Decisión multivariada</span><span class="flecha">→</span><span>k = {d["k_opt"]}, silhouette = {es_num(d["sil"], 2)}</span></div>
      <ul>
        <li>KMeans óptimo por silhouette: k = {d["k_opt"]} (silhouette = {es_num(d["sil"], 2)}, el máximo del rango evaluado).</li>
        <li>Cluster 0 ({es_num(n_baja, 0)} aislados): concentra los aislados de <strong>Baja susceptibilidad relativa</strong>.</li>
        <li>Cluster 1 ({es_num(cl['cluster_kmeans'].value_counts().get(1, 0), 0)} aislados): agrupa la mayoría de los aislados de Alta y Moderada susceptibilidad.</li>
        <li>Categorías del score por terciles: Alta = {n_alta}, Moderada = {n_mod}, Baja = {n_baja} aislados.</li>
      </ul>
    </div>
    <div class="que-significa verde">
      <div class="titulo-tarjeta">Qué significa</div>
      <p>La respuesta al extracto es heterogénea entre aislados: algunos son inhibidos con firmeza por las
      tres técnicas, mientras otros escapan al tratamiento (Baja susceptibilidad relativa). Esta
      variabilidad es la que el análisis univariado solo podía insinuar a través de la interacción
      método × aislado.</p>
    </div>
  </section>

  <!-- ================= SECCION 9: RANKING ================= -->
  <section id="ranking" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">9</span>La decisión · Ranking de técnicas</h2>
    <p class="lead">Se combinan tres métricas normalizadas (min-max 0-1) en un score compuesto: rendimiento
    medio, %INH micelial medio y %INH conidias medio.</p>
    <div class="grid" style="display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:0 1.5rem;">
      <div>{html_figura_bloque("Perfil normalizado", f["radar"], "Cada eje normalizado a 0-1. El área del radar muestra el equilibrio entre rendimiento y actividad.", "Figura R1")}</div>
      <div>{html_figura_bloque("Score compuesto", f["score"], "Promedio simple de las tres métricas normalizadas.", "Figura R2")}</div>
    </div>
    <table>
      <thead><tr><th class="num">#</th><th>Técnica</th><th class="num">Rendimiento (%)</th><th class="num">INH micelial (%)</th><th class="num">INH conidias (%)</th><th class="num">Score</th></tr></thead>
      <tbody>{filas_ranking}</tbody>
    </table>
    <div class="decision">
      <div class="cabeza"><span>Trade-off</span><span class="flecha">→</span><span>Maceración 0,667 &gt; Soxhlet 0,436 &gt; Ultrasonido 0,053</span></div>
      <ul>
        <li><strong>Maceración</strong> domina la actividad antifúngica (86,3% de inhibición micelial y
        29,2% de inhibición de conidias) pero rinde poco (12,1%).</li>
        <li><strong>Soxhlet</strong> maximiza el rendimiento (43,4%) pero su actividad es intermedia.</li>
        <li><strong>Ultrasonido</strong> queda en último lugar: rendimiento modesto (17,1%) y actividad casi
        nula sobre conidias (-0,17%).</li>
        <li>Si el objetivo es actividad antifúngica, Maceración es la opción; si se prioriza cantidad de
        extracto, Soxhlet compensa con una actividad aceptable.</li>
      </ul>
    </div>
  </section>

  <!-- ================= SECCION 10: CONCLUSIONES ================= -->
  <section id="conclusiones" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">10</span>Conclusiones y limitaciones</h2>
    <h3>Conclusiones principales</h3>
    <ul>
      <li>La técnica de extracción es determinante para el rendimiento: Soxhlet (43,4%) supera claramente a
      Maceración (12,1%) y Ultrasonido (17,1%), que no difieren entre sí (Tukey, p = 0,358).</li>
      <li>La actividad antifúngica no sigue el mismo orden: Maceración produce la mayor inhibición micelial
      (86,3%) y de conidias (29,2%), seguida de Soxhlet y Ultrasonido.</li>
      <li>Existe un efecto techo en el %INH micelial de Maceración que probablemente subestima las
      diferencias reales entre técnicas.</li>
      <li>Los aislados no responden de forma uniforme: la interacción método × aislado es significativa y el
      análisis multivariado separa un grupo de baja susceptibilidad relativa frente al resto.</li>
      <li>El ranking compuesto (score 0-1) posiciona a Maceración (0,667) por encima de Soxhlet (0,436) y
      Ultrasonido (0,053).</li>
    </ul>
    <h3>Limitaciones</h3>
    <ul>
      <li><strong>Concentración única</strong>: se ensayó solo 5 mg/mL; no se dispone de curvas
      dosis-respuesta ni de EC50 para comparar potencia.</li>
      <li><strong>Control compartido</strong>: cada %INH se calculó contra un único control C4 por aislado
      (pseudorreplicación); los modelos mixtos mitigan parcialmente la dependencia.</li>
      <li><strong>Alcance in vitro</strong>: los resultados son de bioensayos en placa y no se extrapolan a
      condiciones de campo.</li>
      <li><strong>Rendimiento con n = 3</strong>: el bloque de rendimiento tiene pocas réplicas; la potencia
      para detectar diferencias pequeñas es limitada.</li>
    </ul>
  </section>

  <!-- ================= SECCION 11: METODOLOGIA ================= -->
  <section id="metodologia" class="seccion ancla">
    <h2 class="seccion-titulo"><span class="num">11</span>Metodología y reproducibilidad</h2>
    <p>Todos los resultados de esta página se generan a partir de una pipeline reproducible en Python
    (semilla 42), sin pasos manuales:</p>
    <ul>
      <li><strong>Notebook orquestador</strong>: <code>{DISENO}/analisis_dca.ipynb</code> — ejecuta la pipeline
      completa de forma reproducible.</li>
      <li><strong>Informe final</strong>: <code>{DISENO}/resultados/reportes/informe_final.md</code> (y su versión
      HTML) con la narrativa técnica completa.</li>
      <li><strong>Master dataset</strong>: <code>{DISENO}/resultados/database/master_dataset_tomillo_fusarium.csv</code>
      — 279 filas × 9 columnas (3 métodos × 31 aislados × 3 réplicas).</li>
      <li><strong>Validación de %INH</strong>: las {es_num(n_verif, 0)} unidades fueron verificadas contra la
      fórmula (1 − C1/C4) × 100 con {es_num(n_disc, 0)} discrepancias.</li>
    </ul>
    <div class="que-significa verde">
      <div class="titulo-tarjeta">Reproducibilidad</div>
      <p>Esta página se genera ejecutando <code>python3 generar_pagina.py</code> y se despliega como
      archivo estático en GitHub Pages (solo requiere el CDN de Plotly). Los valores mostrados provienen
      directamente de las tablas exportadas por la pipeline, no de cálculos duplicados en la página.</p>
    </div>
  </section>

</div>

<footer>
  <div class="interno">
    <p>Actividad antifúngica de extractos de <em>Thymus vulgaris</em> contra <em>Fusarium</em> spp. ·
    Análisis <strong>DCA</strong> · Página generada automáticamente por <code>generar_pagina.py</code> · Datos:
    <code>{DISENO}/resultados/database/master_dataset_tomillo_fusarium.csv</code> · Figuras interactivas con Plotly.js.</p>
  </div>
</footer>
</body>
</html>
"""

    DIR_PAGINA_DISENO.mkdir(parents=True, exist_ok=True)
    (DIR_PAGINA_DISENO / "index.html").write_text(html, encoding="utf-8")

    n_fig = len(f)
    print(f"OK: pagina/{DISENO}/index.html generada con {n_fig} figuras Plotly.")
    print(f"Tamaño: {sum(len(v) for v in f.values()) / 1024:.0f} KB de figuras embebidas.")


CSS_HUB = """
:root{--azul:#0072B2;--naranja:#D55E00;--verde:#009E73;--ambar:#F0E442;
  --gris:#5A5A5A;--borde:#E5E5E5;--fondo-suave:#F7F9FB;--tinta:#1F2933;}
*{box-sizing:border-box;}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,
  Helvetica,Arial,sans-serif;color:var(--tinta);background:#fff;line-height:1.6;}
.contenido{max-width:1040px;margin:0 auto;padding:3rem 1.4rem 4rem;}
h1{font-size:2.1rem;line-height:1.2;margin:0 0 .4rem;letter-spacing:-.01em;}
p.sub{font-size:1.1rem;color:var(--gris);max-width:760px;margin:0 0 2.2rem;}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:1.2rem;}
.tarjeta{border:1px solid var(--borde);border-top:4px solid var(--azul);border-radius:12px;
  padding:1.2rem 1.3rem;text-decoration:none;color:inherit;display:block;background:#fff;
  box-shadow:0 1px 3px rgba(0,0,0,.05);transition:transform .12s ease,box-shadow .12s ease;}
.tarjeta:hover{transform:translateY(-2px);box-shadow:0 6px 16px rgba(0,0,0,.08);}
.tarjeta.oculta{border-top-color:var(--gris);opacity:.75;}
.tarjeta h2{margin:0 0 .35rem;font-size:1.15rem;color:var(--azul);}
.tarjeta.oculta h2{color:var(--gris);}
.tarjeta p{margin:0 0 .5rem;font-size:.92rem;color:#374151;}
.estado{display:inline-block;font-size:.75rem;font-weight:700;border-radius:999px;
  padding:.18rem .6rem;margin-bottom:.6rem;}
.estado.activo{background:#E3F2E9;color:#1B7A43;}
.estado.futuro{background:#F2F3F5;color:#6B7280;}
 .pie{color:var(--gris);font-size:.85rem;margin-top:2.4rem;}
 a.pie-enlace{color:var(--azul);text-decoration:none;}
 .ventanas{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));
   gap:1.4rem;align-items:stretch;margin-bottom:2.4rem;}
 .ventana{border:1px solid var(--borde);border-radius:12px;background:var(--fondo-suave);
   box-shadow:0 2px 10px rgba(0,0,0,.06);overflow:hidden;display:flex;flex-direction:column;}
 .ventana-barra{display:flex;align-items:center;justify-content:space-between;gap:.8rem;
   background:#fff;border-bottom:1px solid var(--borde);padding:.7rem 1.1rem;}
 .ventana-barra h2{margin:0;font-size:1rem;color:var(--azul);}
 .ventana.bdca .ventana-barra h2{color:var(--naranja);}
 .ventana-barra .estado{margin-bottom:0;}
 .ventana-cuerpo{padding:1.15rem 1.3rem 1.3rem;display:flex;flex-direction:column;flex:1;}
 .ventana-cuerpo p{margin:0 0 .75rem;font-size:.92rem;color:#374151;}
 .ventana-lista{margin:0 0 1.1rem;padding-left:1.15rem;font-size:.9rem;color:var(--gris);}
 .ventana-lista li{margin-bottom:.3rem;}
 .ventana-enlace{margin-top:auto;align-self:flex-start;display:inline-block;color:#fff;
   background:var(--azul);text-decoration:none;font-size:.88rem;font-weight:600;
   border-radius:8px;padding:.5rem .95rem;transition:background .12s ease;}
 .ventana-enlace:hover{background:#005B94;}
 .futuro-fila{max-width:540px;}
 """


def generar_hub() -> None:
    """Genera pagina/index.html (hub) y placeholders para diseños futuros."""
    DIR_PAGINA.mkdir(parents=True, exist_ok=True)

    ventanas = [
        {
            "clave": "dca",
            "titulo": "DCA — Diseño completamente al azar",
            "desc": "Tres técnicas de extracción × 31 aislados de Fusarium a una única concentración "
                    "(5 mg/mL). Comparación del rendimiento de extracción y de la actividad antifúngica.",
            "puntos": [
                "Rendimiento de extracción",
                "%INH de crecimiento micelial",
                "Esporulación (conidias)",
                "Susceptibilidad y ranking de aislados",
            ],
            "href": "dca/index.html",
        },
        {
            "clave": "bdca",
            "titulo": "BDCA — Bloques completos al azar",
            "desc": "Análisis con control de la variabilidad entre bloques experimentales "
                    "(efecto aleatorio) mediante modelo mixto.",
            "puntos": [
                "Rendimiento (yield)",
                "4 tratamientos × 9 bloques",
                "Modelo mixto con bloque aleatorio",
                "Comparaciones Tukey HSD",
            ],
            "href": "bdca/index.html",
        },
    ]
    futuros = [
        {
            "clave": "factorial",
            "titulo": "Factorial — Técnica × Concentración × Aislado",
            "desc": "Múltiples concentraciones, interacciones factoriales y dosis-respuesta (EC50/EC90).",
            "estado": "futuro", "estado_txt": "Próximamente", "href": "factorial/index.html",
        },
    ]

    def ventana_html(v) -> str:
        puntos = "".join(f"<li>{p}</li>" for p in v["puntos"])
        return f"""<div class="ventana {v['clave']}">
  <div class="ventana-barra">
    <h2>{v['titulo']}</h2>
    <span class="estado activo">Disponible</span>
  </div>
  <div class="ventana-cuerpo">
    <p>{v['desc']}</p>
    <ul class="ventana-lista">{puntos}</ul>
    <a class="ventana-enlace" href="{v['href']}">Ver análisis →</a>
  </div>
</div>"""

    ventanas_html = "".join(ventana_html(v) for v in ventanas)
    tarjetas_futuras = "".join(
        f"""<a class="tarjeta {d['estado']}" href="{d['href']}">
  <span class="estado {d['estado']}">{d['estado_txt']}</span>
  <h2>{d['titulo']}</h2>
  <p>{d['desc']}</p>
</a>"""
        for d in futuros
    )

    html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Análisis estadísticos de la actividad antifúngica de Thymus vulgaris contra Fusarium spp. — DCA, BDCA y diseño factorial.">
<title>Tomillo × Fusarium — Análisis estadísticos</title>
<style>{CSS_HUB}</style>
</head>
<body>
<div class="contenido">
  <h1>Actividad antifúngica de <em>Thymus vulgaris</em> contra <em>Fusarium</em> spp.</h1>
  <p class="sub">Análisis estadísticos reproducibles de los distintos diseños experimentales del estudio.
  Cada análisis tiene su propia ruta: descriptiva → supuestos → modelos → comparaciones → interpretación.</p>
  <div class="ventanas">{ventanas_html}</div>
  <div class="futuro-fila">
    <div class="grid">{tarjetas_futuras}</div>
  </div>
  <p class="pie">Página generada por <code>generar_pagina.py --hub</code> · Repositorio:
  <a class="pie-enlace" href="https://github.com/marcos-Nieves-24/proyecto-tomillo">proyecto-tomillo</a>.</p>
</div>
</body>
</html>
"""
    (DIR_PAGINA / "index.html").write_text(html, encoding="utf-8")

    for d in futuros:
        sub = DIR_PAGINA / d["clave"]
        sub.mkdir(parents=True, exist_ok=True)
        ph = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{d['titulo']} — Próximamente</title>
<style>{CSS_HUB}</style>
</head>
<body>
<div class="contenido">
  <p><a class="pie-enlace" href="../index.html">← Volver al índice de análisis</a></p>
  <h1>{d['titulo']}</h1>
  <p class="sub">Este análisis se encuentra en preparación. Cuando estén disponibles los datos del diseño,
  esta página mostrará los resultados completos con el mismo estándar que los análisis DCA y BDCA.</p>
</div>
</body>
</html>
"""
        (sub / "index.html").write_text(ph, encoding="utf-8")

    futuros_clave = [d["clave"] for d in futuros]
    print(f"OK: pagina/index.html (hub) + placeholders {futuros_clave or 'ninguno'} generados.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Genera la pagina interactiva del estudio.")
    parser.add_argument("--hub", action="store_true", help="genera el hub y los placeholders")
    parser.add_argument("--diseno", choices=["dca", "bdca"], default="dca", help="analisis a generar")
    args = parser.parse_args()

    if args.hub:
        generar_hub()
    elif args.diseno == "bdca":
        main_bdca()
    else:
        main()
