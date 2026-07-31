"""Fase 10: analisis multivariado de susceptibilidad.

Construye una matriz por aislado (31 x 6) con las medias de %INH micelial y
%INH conidias para cada metodo, estandariza (z-score), y aplica: PCA, clus-
tering jerarquico de Ward (con coeficiente cofenetico), KMeans con k optimo
por silhouette, categorias biologicas de susceptibilidad relativa (Alta /
Moderada / Baja) y visualizaciones (scree, biplot, dendrograma, heatmap con
dendrograma y scatter de clusters).

Regla cientifica del proyecto: NUNCA se usa la palabra "resistente" sin un
criterio validado; se emplean las etiquetas de susceptibilidad relativa.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from pipeline.config import (
    METODOS,
    METODO_LABEL,
    PALETA_METODOS,
    SEMILLA_ALEATORIA,
    guardar_tabla,
    save_figure_pub,
)

VARIABLES_MULTIVARIADAS = [
    f"inhib_micelial_{m}" for m in METODOS
] + [f"inhib_conidias_{m}" for m in METODOS]

CATEGORIAS = [
    "Alta susceptibilidad relativa",
    "Moderada susceptibilidad relativa",
    "Baja susceptibilidad relativa",
]
COLOR_CATEGORIAS = {
    "Alta susceptibilidad relativa": "#009E73",
    "Moderada susceptibilidad relativa": "#F0E442",
    "Baja susceptibilidad relativa": "#D55E00",
}


def _matriz_por_aislado(df_bio: pd.DataFrame) -> pd.DataFrame:
    """Matriz aislado x 6 metricas (medias por metodo y por tipo de %INH)."""
    filas = []
    for aislado, grupo in df_bio.groupby("aislamiento"):
        fila = {"aislamiento": aislado}
        for m in METODOS:
            g = grupo[grupo["metodo_extraccion"] == m]
            fila[f"inhib_micelial_{m}"] = g["porcentaje_inhibicion_micelial"].mean()
            fila[f"inhib_conidias_{m}"] = g["porcentaje_inhibicion_conidias"].mean()
        filas.append(fila)
    matriz = pd.DataFrame(filas).set_index("aislamiento")
    matriz.columns.name = "metrica"
    return matriz


def _estandarizar(matriz: pd.DataFrame) -> pd.DataFrame:
    """Estandarizacion z-score por columna (media 0, desv. 1)."""
    return (matriz - matriz.mean()) / matriz.std(ddof=0)


def _score_susceptibilidad(matriz_z: pd.DataFrame) -> pd.Series:
    """Score compuesto: promedio de los z de las 6 metricas por aislado."""
    return matriz_z.mean(axis=1)


def _categorias_terciles(score: pd.Series) -> pd.Series:
    """Terciles del score -> etiquetas de susceptibilidad relativa."""
    ranks = score.rank(method="first")
    n = len(ranks)
    tercil = pd.qcut(ranks, 3, labels=False) if n >= 6 else pd.Series(0, index=ranks.index)
    etiquetas = tercil.map({0: CATEGORIAS[2], 1: CATEGORIAS[1], 2: CATEGORIAS[0]})
    return etiquetas.astype("str")


def _kmeans_optimo(matriz_z: pd.DataFrame, kmax=8) -> tuple[int, pd.DataFrame]:
    """Codo (inercia) + silhouette para k=2..kmax; k optimo = max silhouette."""
    inercias = []
    for k in range(1, kmax + 1):
        km = KMeans(n_clusters=k, random_state=SEMILLA_ALEATORIA, n_init=10)
        km.fit(matriz_z)
        inercias.append(km.inertia_)

    siluetas = []
    for k in range(2, kmax + 1):
        km = KMeans(n_clusters=k, random_state=SEMILLA_ALEATORIA, n_init=10)
        etiquetas = km.fit_predict(matriz_z)
        siluetas.append(silhouette_score(matriz_z, etiquetas))
    k_opt = int(range(2, kmax + 1)[np.argmax(siluetas)])

    tabla = pd.DataFrame({
        "k": list(range(1, kmax + 1)),
        "inercia": [round(v, 3) for v in inercias],
    })
    tabla.loc[tabla["k"].isin(range(2, kmax + 1)), "silhouette"] = [round(v, 4) for v in siluetas]

    return k_opt, tabla


def analisis_multivariado(df_bio: pd.DataFrame) -> dict:
    """Ejecuta el analisis multivariado completo y guarda tablas y figuras."""
    matriz = _matriz_por_aislado(df_bio)
    matriz_z = _estandarizar(matriz)
    score = _score_susceptibilidad(matriz_z)
    categoria = _categorias_terciles(score)

    # 1. PCA ----------------------------------------------------------------
    pca = PCA(n_components=len(VARIABLES_MULTIVARIADAS), random_state=SEMILLA_ALEATORIA)
    scores = pca.fit_transform(matriz_z)
    varianza = pca.explained_variance_ratio_
    loadings = pd.DataFrame(
        pca.components_.T, index=matriz_z.columns,
        columns=[f"PC{i + 1}" for i in range(pca.n_components_)],
    )
    scores_df = pd.DataFrame(
        scores, index=matriz_z.index, columns=[f"PC{i + 1}" for i in range(scores.shape[1])]
    ).reset_index()
    scores_df["aislamiento"] = scores_df["aislamiento"]
    guardar_tabla(scores_df, "susceptibilidad_pca_scores", index=False)
    guardar_tabla(loadings.reset_index().rename(columns={"index": "metrica"}),
                  "susceptibilidad_pca_loadings", index=False)

    # Scree plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(range(1, len(varianza) + 1), varianza, color="#4C72B0", alpha=0.85,
           label="Varianza explicada")
    ax.plot(range(1, len(varianza) + 1), np.cumsum(varianza), marker="o",
            color="#C44E52", label="Acumulada")
    ax.axhline(0.8, color="grey", ls="--", lw=0.9)
    ax.set_xlabel("Componente principal")
    ax.set_ylabel("Proporcion de varianza")
    ax.set_xticks(range(1, len(varianza) + 1))
    ax.legend()
    save_figure_pub(fig, "multivariado_pca_scree", titulo="PCA - Varianza explicada por componente")

    # Biplot PC1-PC2
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(scores[:, 0], scores[:, 1], s=45, alpha=0.8, color="#4C72B0", zorder=3)
    for i, aislado in enumerate(matriz_z.index):
        ax.annotate(aislado, (scores[i, 0], scores[i, 1]),
                    fontsize=7, alpha=0.85, xytext=(4, 4), textcoords="offset points")
    escala = max(np.abs(scores[:, :2]).max() / np.abs(loadings.iloc[:, :2]).max().max(), 1.0) * 0.9
    for metrica in loadings.index:
        ax.annotate(
            "", xy=(loadings.loc[metrica, "PC1"] * escala, loadings.loc[metrica, "PC2"] * escala),
            xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="#C44E52", lw=1.5),
        )
        ax.text(loadings.loc[metrica, "PC1"] * escala * 1.1,
                loadings.loc[metrica, "PC2"] * escala * 1.1,
                _etiqueta_metrica(metrica), fontsize=8, color="#C44E52")
    ax.set_xlabel(f"PC1 ({varianza[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({varianza[1] * 100:.1f}%)")
    ax.axhline(0, color="grey", lw=0.7)
    ax.axvline(0, color="grey", lw=0.7)
    save_figure_pub(fig, "multivariado_pca_biplot", titulo="PCA - Biplot PC1-PC2 (scores y loadings)")

    # 2. Clustering jerarquico (Ward) ------------------------------------------
    distancias = ssd.pdist(matriz_z.values, metric="euclidean")
    enlace = sch.linkage(distancias, method="ward")
    cofenetico = float(sch.cophenet(enlace, distancias)[0])

    fig, ax = plt.subplots(figsize=(11, 6))
    sch.dendrogram(enlace, labels=matriz_z.index, ax=ax, color_threshold=0,
                   leaf_font_size=8)
    ax.set_ylabel("Distancia euclidiana")
    ax.set_xlabel("Aislado")
    save_figure_pub(fig, "multivariado_dendrograma",
                    titulo=f"Clustering jerárquico Ward (coef. cofenético={cofenetico:.3f})")

    # 3. KMeans ----------------------------------------------------------------
    k_opt, tabla_kmeans = _kmeans_optimo(matriz_z)
    guardar_tabla(tabla_kmeans, "susceptibilidad_kmeans_metricas", index=False)
    kmeans = KMeans(n_clusters=k_opt, random_state=SEMILLA_ALEATORIA, n_init=20)
    cluster = kmeans.fit_predict(matriz_z)

    fig, ejes = plt.subplots(1, 2, figsize=(13, 5))
    ejes[0].plot(range(1, len(tabla_kmeans) + 1), tabla_kmeans["inercia"], marker="o",
                 color="#4C72B0")
    ejes[0].set_xlabel("k")
    ejes[0].set_ylabel("Inercia")
    ejes[0].set_title("Metodo del codo")
    ejes[1].plot(range(2, len(tabla_kmeans) + 1), tabla_kmeans["silhouette"].dropna(),
                 marker="o", color="#C44E52")
    ejes[1].axvline(k_opt, color="grey", ls="--", lw=0.9)
    ejes[1].set_xlabel("k")
    ejes[1].set_ylabel("Silhouette")
    ejes[1].set_title(f"Silhouette (k optimo = {k_opt})")
    save_figure_pub(fig, "multivariado_kmeans_optimo", titulo="Selección del número de clusters (KMeans)")

    # 4. Tabla integradora -------------------------------------------------------
    tabla_final = pd.DataFrame({
        "aislamiento": matriz_z.index,
        "score_susceptibilidad": score.values,
        "categoria_susceptibilidad": categoria.values,
        "cluster_kmeans": cluster,
    })
    for metrica in VARIABLES_MULTIVARIADAS:
        tabla_final[metrica] = matriz[metrica].values
    tabla_final = tabla_final.sort_values("score_susceptibilidad", ascending=False)
    guardar_tabla(tabla_final, "susceptibilidad_clusters", index=False)

    # 5. Heatmap con dendrograma y categoria -------------------------------------
    orden_categoria = pd.Series(COLOR_CATEGORIAS).loc[
        [c for c in CATEGORIAS if c in categoria.unique()]
    ].to_dict()
    filas_color = categoria.map(COLOR_CATEGORIAS).to_frame("categoria")
    g = sns.clustermap(
        matriz_z, method="ward", metric="euclidean", cmap="RdBu_r", center=0,
        row_colors=filas_color, figsize=(9, 10),
        yticklabels=True, xticklabels=[_etiqueta_metrica(c) for c in matriz_z.columns],
    )
    g.ax_heatmap.set_xlabel("Metrica estandarizada")
    handles = [plt.Line2D([0], [0], color=v, lw=6, label=k) for k, v in orden_categoria.items()]
    g.ax_heatmap.legend(handles=handles, loc="lower right", title="Susceptibilidad relativa")
    save_figure_pub(g.fig, "multivariado_heatmap",
                    titulo="Heatmap de metricas estandarizadas por aislado (Ward)")

    # 6. Scatter PC1-PC2 coloreado por cluster --------------------------------------
    fig, ax = plt.subplots(figsize=(9, 7))
    paleta_cluster = sns.color_palette("Set2", n_colors=k_opt)
    for c in range(k_opt):
        sel = cluster == c
        ax.scatter(scores[sel, 0], scores[sel, 1], s=55, alpha=0.85,
                   color=paleta_cluster[c], label=f"Cluster {c + 1}", zorder=3)
    for i, aislado in enumerate(matriz_z.index):
        ax.annotate(aislado, (scores[i, 0], scores[i, 1]),
                    fontsize=7, alpha=0.8, xytext=(4, 4), textcoords="offset points")
    ax.axhline(0, color="grey", lw=0.7)
    ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel(f"PC1 ({varianza[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({varianza[1] * 100:.1f}%)")
    ax.legend(title=f"KMeans (k={k_opt})")
    save_figure_pub(fig, "multivariado_clusters_pca",
                    titulo=f"Clusters KMeans (k={k_opt}) en el plano PC1-PC2")

    # 7. Cruzar categorias y clusters ----------------------------------------------
    cruce = pd.crosstab(
        tabla_final["cluster_kmeans"], tabla_final["categoria_susceptibilidad"]
    )
    guardar_tabla(cruce, "susceptibilidad_cruce_cluster_categoria", index=True)

    print("=" * 72)
    print("FASE 10 - Análisis multivariado")
    print("=" * 72)
    print(f"Varianza explicada PC1-PC2: {varianza[0]*100:.1f}% / {varianza[1]*100:.1f}%")
    print(f"Coeficiente cofenetico (Ward): {cofenetico:.4f}")
    print(f"k optimo KMeans (silhouette): {k_opt}")
    print("Cruce cluster x categoria:")
    print(cruce.to_string())

    return {
        "matriz": matriz,
        "matriz_z": matriz_z,
        "score_susceptibilidad": score,
        "categoria": categoria,
        "pca": pca,
        "varianza": varianza,
        "loadings": loadings,
        "scores": scores_df,
        "cofenetico": cofenetico,
        "k_optimo": k_opt,
        "cluster": cluster,
        "tabla_final": tabla_final,
        "tabla_kmeans": tabla_kmeans,
        "cruce": cruce,
    }


def _etiqueta_metrica(metrica: str) -> str:
    """Convierte 'inhib_micelial_maceracion' en una etiqueta legible."""
    tipo, metodo = metrica.replace("inhib_", "").rsplit("_", 1)
    tipo_label = "INH micelial" if tipo == "micelial" else "INH conidias"
    return f"{tipo_label} - {METODO_LABEL[metodo]}"
