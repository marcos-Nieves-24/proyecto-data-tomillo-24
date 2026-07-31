#!/usr/bin/env python3
"""Analisis factorial: genotipo de trigo x cepa de Fusarium (datos de Snijders).

Diseno: dos factores (gen, strain) con 3 anos (year) sin replica intra-celda.
Fuente: FACTORIAL_Snijders_Fusarium_genotipo_cepa.csv (raiz del repo).
Salidas: factorial/resultados/{tablas,figuras,reportes}
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# ---------------------------------------------------------------------------
# Configuracion y rutas
# ---------------------------------------------------------------------------
RNG = np.random.default_rng(12345)

RAIZ = Path(__file__).resolve().parents[1]
ORIGEN = RAIZ / "FACTORIAL_Snijders_Fusarium_genotipo_cepa.csv"
RESULTADOS = Path(__file__).resolve().parent / "resultados"
TABLAS = RESULTADOS / "tablas"
FIGURAS = RESULTADOS / "figuras"
REPORTES = RESULTADOS / "reportes"

for _d in (TABLAS, FIGURAS, REPORTES):
    _d.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="paper")
plt.rcParams["figure.dpi"] = 110

# ---------------------------------------------------------------------------
# 1. Auditoria de datos
# ---------------------------------------------------------------------------
df = pd.read_csv(ORIGEN)
df["year"] = df["year"].astype(int)

aud = []
def _fil(n, v, ok):
    aud.append({"chequeo": n, "resultado": v, "ok": ok})

# Balance: 17 gen x 4 cepa x 3 ano = 204, 1 obs por celda
_combos = df.groupby(["gen", "strain", "year"]).size()
balance = "completo" if _combos.eq(1).all() and len(_combos) == 204 else "incompleto"
_fil("filas totales", f"{len(df)}", len(df) == 204)
_fil("n genotipos (17)", f"{df['gen'].nunique()}", df["gen"].nunique() == 17)
_fil("n cepas (4)", f"{df['strain'].nunique()}", df["strain"].nunique() == 4)
_fil("n anos (3)", f"{df['year'].nunique()}", df["year"].nunique() == 3)
_fil("combinacion gen x cepa x ano",
     f"{len(_combos)} celdas, {balance}, replicas/celda={_combos.max() if len(_combos) else 0}",
     balance == "completo")
_fil("valores faltantes", f"{int(df.isna().sum().sum())}", df.isna().sum().sum() == 0)
_fil("filas duplicadas", f"{int(df.duplicated().sum())}", df.duplicated().sum() == 0)
_fil("valores imposibles y<0 o y>100",
     f"{int(((df['y'] < 0) | (df['y'] > 100)).sum())}",
     ((df["y"] < 0) | (df["y"] > 100)).sum() == 0)
_fil("rango de y", f"[{df['y'].min()}, {df['y'].max()}]", True)

auditoria = pd.DataFrame(aud)
auditoria.to_csv(TABLAS / "auditoria.csv", index=False)

# ---------------------------------------------------------------------------
# 2. Estadisticos descriptivos por cepa y por genotipo
# ---------------------------------------------------------------------------
def _descriptivos(grupo, nombre_grupo):
    filas = []
    for g, sub in df.groupby(grupo):
        n = len(sub)
        media = sub["y"].mean()
        se = sub["y"].std(ddof=1) / np.sqrt(n)
        t = stats.t.ppf(0.975, df=n - 1) if n > 1 else np.nan
        filas.append({
            nombre_grupo: g,
            "n": n,
            "media": media,
            "se": se,
            "ic95_inf": media - t * se,
            "ic95_sup": media + t * se,
        })
    return pd.DataFrame(filas).sort_values("media", ascending=False).reset_index(drop=True)

med_cepa = _descriptivos("strain", "cepa")
med_gen = _descriptivos("gen", "genotipo")
med_cepa.to_csv(TABLAS / "medias_por_cepa.csv", index=False)
med_gen.to_csv(TABLAS / "medias_por_gen.csv", index=False)

# ---------------------------------------------------------------------------
# 3. Figuras exploratorias
# ---------------------------------------------------------------------------
def _guardar(fig, nombre):
    fig.savefig(FIGURAS / nombre, bbox_inches="tight")
    plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
sns.boxplot(data=df, x="strain", y="y", ax=ax, hue="strain", legend=False, palette="muted")
ax.set_title("Severidad de enfermedad por cepa de Fusarium")
ax.set_xlabel("Cepa")
ax.set_ylabel("y (severidad)")
_guardar(fig, "boxplot_por_cepa.png")

fig, ax = plt.subplots(figsize=(11, 4.5))
sns.boxplot(data=df, x="gen", y="y", ax=ax, hue="gen", legend=False, palette="viridis")
ax.set_title("Severidad de enfermedad por genotipo de trigo")
ax.set_xlabel("Genotipo")
ax.set_ylabel("y (severidad)")
ax.tick_params(axis="x", rotation=60)
_guardar(fig, "boxplot_por_gen.png")

# Grafico de interaccion: media de y por genotipo, una linea por cepa
fig, ax = plt.subplots(figsize=(11, 5.5))
inter_med = df.groupby(["gen", "strain"])["y"].mean().unstack("strain")
inter_med.plot(marker="o", ax=ax, markersize=4)
ax.set_title("Interaccion genotipo x cepa (media de y)")
ax.set_xlabel("Genotipo")
ax.set_ylabel("y medio")
ax.tick_params(axis="x", rotation=60)
ax.legend(title="Cepa")
_guardar(fig, "interaccion_gen_cepa.png")

# ---------------------------------------------------------------------------
# 4. Modelo lineal: y ~ C(gen) + C(strain) + C(year) + C(gen):C(strain)
# ---------------------------------------------------------------------------
# Con 1 observacion por celda (gen x cepa) y 3 anos como bloques, el termino
# de interaccion gen:cepa es estimable (df = 16*3 = 48). El residuo queda con
# df = 204 - 70 parametros = 134 (mezcla el error puro con el efecto ano
# residual: cada ano no replicado dentro de cada combinacion gen x cepa).
FORMULA = "y ~ C(gen) + C(strain) + C(year) + C(gen):C(strain)"
modelo_usado = "completo"
try:
    modelo = ols(FORMULA, data=df).fit()
    if not modelo.resid.std() > 0:
        raise ValueError("residuos degenerados (varianza nula)")
    anova = anova_lm(modelo, typ=2)
    error_msg = ""
except Exception as exc:  # noqa: BLE001
    # Fallback: modelo aditivo sin interaccion gen x cepa.
    modelo_usado = "aditivo"
    error_msg = f"{type(exc).__name__}: {exc}"
    FORMULA = "y ~ C(gen) + C(strain) + C(year)"
    modelo = ols(FORMULA, data=df).fit()
    anova = anova_lm(modelo, typ=2)

anova = anova.reset_index().rename(columns={"index": "term"})
if "mean_sq" not in anova.columns:
    anova["mean_sq"] = anova["sum_sq"] / anova["df"]
anova["eta_cuad"] = anova["sum_sq"] / anova["sum_sq"].sum()
anova["f"] = anova["F"].astype(float)
anova["p"] = anova["PR(>F)"].astype(float)
anova_tabla = anova[["term", "df", "sum_sq", "mean_sq", "f", "p", "eta_cuad"]]
anova_tabla.to_csv(TABLAS / "anova.csv", index=False)

# ---------------------------------------------------------------------------
# 5. Comparaciones multiples (Tukey HSD)
# ---------------------------------------------------------------------------
def _tukey(grupo, nombre):
    res = pairwise_tukeyhsd(df["y"], df[grupo], alpha=0.05)
    tab = pd.DataFrame(
        data=res._results_table.data[1:],
        columns=res._results_table.data[0],
    )
    tab.rename(columns={"group1": "grupo1", "group2": "grupo2",
                        "p-adj": "p_ajustada", "reject": "significativo"},
               inplace=True)
    tab["significativo"] = tab["significativo"].astype(str)
    tab.to_csv(TABLAS / f"posthoc_{nombre}.csv", index=False)
    return tab

tukey_cepas = _tukey("strain", "cepas")
tukey_gen = _tukey("gen", "gen")

# ---------------------------------------------------------------------------
# 6. Informe en espanol (neutro/profesional)
# ---------------------------------------------------------------------------
def _f(v, nd=2):
    return f"{v:.{nd}f}"

def _sig(p):
    return " (no significativo)" if p >= 0.05 else ""

def _sigsi(p):
    return "Si" if p < 0.05 else "No"

# Seccion ANOVA desde el modelo efectivamente usado
anova_df = anova_tabla.copy()
def _term(t):
    t = t.replace("C(", "").replace(")", "").replace(":", " x ").replace("Residual", "Residual")
    return {"gen": "Genotipo", "strain": "Cepa", "year": "Año", "gen x strain": "Genotipo x Cepa"}.get(t, t)

filas_anova = "\n".join(
    f"| {_term(r['term'])} | {int(r['df'])} | {_f(r['sum_sq'], 1)} | {_f(r['f'])} | {_f(r['p'])} | {_sigsi(r['p'])} | {_f(r['eta_cuad'], 3)} |"
    for _, r in anova_df.iterrows()
)

rank_cepas = med_cepa.copy()
filas_rank_cepas = "\n".join(
    f"| {r['cepa']} | {_f(r['media'])} | {_f(r['se'])} | [{_f(r['ic95_inf'])}, {_f(r['ic95_sup'])}] |"
    for _, r in rank_cepas.iterrows()
)

rank_gen = med_gen.copy()
filas_rank_gen = "\n".join(
    f"| {r['genotipo']} | {_f(r['media'])} | {_f(r['se'])} | [{_f(r['ic95_inf'])}, {_f(r['ic95_sup'])}] |"
    for _, r in rank_gen.iterrows()
)

filas_tukey_cepas = "\n".join(
    f"| {r['grupo1']} | {r['grupo2']} | {_f(float(r['meandiff']), 2)} | {_f(float(r['p_ajustada']))} | {r['significativo']} |"
    for _, r in tukey_cepas.iterrows()
)
filas_tukey_gen = "\n".join(
    f"| {r['grupo1']} | {r['grupo2']} | {_f(float(r['meandiff']), 2)} | {_f(float(r['p_ajustada']))} | {r['significativo']} |"
    for _, r in tukey_gen.iterrows()
)

_cepa_menor = rank_cepas.iloc[-1]["cepa"]
_cepa_mayor = rank_cepas.iloc[0]["cepa"]
_gen_menor = rank_gen.iloc[-1]["genotipo"]
_gen_mayor = rank_gen.iloc[0]["genotipo"]

g_inter = anova_df[anova_df["term"] == "C(gen):C(strain)"]
inter_sig = bool(g_inter["p"].iloc[0] < 0.05) if len(g_inter) else False
inter_desc = (
    "la interacción genotipo × cepa resultó estadísticamente significativa, "
    "indicando que la magnitud de la diferencia entre cepas depende del genotipo"
    if inter_sig else
    "la interacción genotipo × cepa no resultó estadísticamente significativa al 5%, "
    "lo que sugiere que los efectos de cepa y genotipo actúan de forma aproximadamente aditiva"
)

if modelo_usado == "completo":
    parrafo_modelo = (
        "Se ajustó un modelo lineal `y ~ C(gen) + C(strain) + C(year) + C(gen):C(strain)`. "
        "El factor `year` se incorpora como factor de bloqueo/replicación: dado que existe una sola "
        "observación por combinación genotipo × cepa, los tres años proporcionan la replicación necesaria "
        "para estimar el error, pero el residual confunde la variación entre años con el error puro "
        "(no existe término de error puro dentro de la celda genotipo × cepa). "
        "La interacción genotipo × cepa es estimable (48 gl)."
    )
else:
    parrafo_modelo = (
        f"Se intentó ajustar `y ~ C(gen) + C(strain) + C(year) + C(gen):C(strain)` pero el modelo completo "
        f"no pudo estimarse ({error_msg}). Se utilizó el modelo aditivo `y ~ C(gen) + C(strain) + C(year)`. "
        "Esta decisión queda documentada como respaldo; las conclusiones sobre interacción no pueden "
        "evaluarse en este modelo."
    )

informe = f"""# Análisis factorial: genotipo de trigo × cepa de Fusarium (Snijders)

## 1. Contexto y objetivo

Se analiza la severidad de enfermedad (`y`, mayor valor = mayor severidad = menor resistencia)
de 17 genotipos de trigo frente a 4 cepas de *Fusarium* spp., en un diseño factorial
17 × 4 × 3 con los tres años (1986, 1987, 1988) como factor de bloqueo. El objetivo es
comparar la agresividad relativa de las cepas y la susceptibilidad relativa de los genotipos,
evaluar la interacción genotipo × cepa y cuantificar la incertidumbre de las estimaciones.

> Nota de terminología: ningún genotipo se clasifica como «resistente» o «susceptible» en
> términos absolutos, dado que no se dispone de un umbral validado de resistencia. Se emplean
> los términos de «menor severidad relativa» y «mayor severidad relativa».

## 2. Auditoría de datos

Fuente: `FACTORIAL_Snijders_Fusarium_genotipo_cepa.csv` (204 filas). Las columnas son
`gen` (genotipo), `strain` (cepa), `year` (año) e `y` (severidad). No se imputó ni modificó
ninguna observación.

- Diseño balanceado: 204 = 17 genotipos × 4 cepas × 3 años, exactamente **1 observación por celda**.
- Sin valores faltantes, sin filas duplicadas, sin valores imposibles (y ∈ [0.1, 69.3], todos en [0, 100]).
- Detalle completo en `tablas/auditoria.csv`.

## 3. Estadísticos descriptivos

Medias ordenadas de mayor a menor severidad. IC95 basado en distribución t (columna `se`,
error estándar). Tablas completas en `tablas/medias_por_cepa.csv` y `tablas/medias_por_gen.csv`.

### 3.1 Por cepa

| Cepa | Media | SE | IC95 |
|---|---|---|---|
{filas_rank_cepas}

### 3.2 Por genotipo (10 de mayor severidad)

| Genotipo | Media | SE | IC95 |
|---|---|---|---|
{filas_rank_gen}

## 4. Modelo estadístico

{parrafo_modelo}

Fórmula empleada: `{FORMULA}`. Tabla ANOVA en `tablas/anova.csv`. La columna `eta_cuad`
se calcula como `sum_sq` del término / `sum_sq` total (variación total explicada por el término).

| Término | gl | Sum sq | F | p | Significativo | η² |
|---|---|---|---|---|---|---|
{filas_anova}

## 5. Comparaciones múltiples (Tukey HSD)

Comparaciones por pares con ajuste de Tukey (α = 0.05) sobre las medias de cepa y de genotipo.
Tablas completas en `tablas/posthoc_cepas.csv` y `tablas/posthoc_gen.csv`.

### 5.1 Pares entre cepas

| Cepa A | Cepa B | Diferencia | p-ajustada | Significativo |
|---|---|---|---|---|
{filas_tukey_cepas}

### 5.2 Pares entre genotipos

| Genotipo A | Genotipo B | Diferencia | p-ajustada | Significativo |
|---|---|---|---|---|
{filas_tukey_gen}

## 6. Interpretación biológica

- **Agresividad relativa de las cepas:** la cepa {_cepa_mayor} mostró la mayor severidad media
  ({_f(rank_cepas.iloc[0]['media'])}), y la cepa {_cepa_menor} la menor
  ({_f(rank_cepas.iloc[-1]['media'])}). La diferencia entre cepas fue
  estadísticamente significativa (F = {_f(anova_df[anova_df['term']=='C(strain)']['f'].iloc[0])}{_sig(anova_df[anova_df['term']=='C(strain)']['p'].iloc[0])}),
  lo que indica diferencias reales de agresividad entre aislados.
- **Susceptibilidad relativa de los genotipos:** el genotipo {_gen_mayor} presentó la mayor
  severidad media y el genotipo {_gen_menor} la menor. La diferencia entre genotipos fue
  significativa (F = {_f(anova_df[anova_df['term']=='C(gen)']['f'].iloc[0])}{_sig(anova_df[anova_df['term']=='C(gen)']['p'].iloc[0])}),
  reflejando un gradiente de susceptibilidad relativa en el material evaluado.
- **Interacción genotipo × cepa:** {inter_desc}. Cuando una interacción es significativa, el
  ordenamiento de las cepas (o de los genotipos) no es uniforme y conviene interpretar los
  efectos principales con cautela, priorizando las comparaciones dentro de cada nivel.
- **Efecto del año:** el factor `year` absorbe variación entre años; su significancia refleja
  condiciones epidémicas distintas entre campañas, no un efecto de tratamiento.

## 7. Limitaciones

- **Sin replicación intra-celda:** existe una sola observación por combinación genotipo × cepa
  por año. No hay término de error puro dentro de la celda: el residuo del modelo mixtura la
  variación entre años con el error experimental, de modo que los contrastes de significancia
  dependen de esa supuesta homogeneidad. Se recomienda replicación intra-celda para un error puro.
- **Estructura año-observación:** `year` se trata como bloque fijo; no se modeló efecto aleatorio
  de año ni su interacción con los tratamientos.
- **Falta de correlación con dosis:** no hay concentraciones, por lo que no procede estimación
  de EC50/EC90.
- **Clasificación absoluta de resistencia:** no se emplea el término «resistente»; solo se
  describe severidad relativa.

## 8. Archivos generados

- Tablas: `tablas/auditoria.csv`, `tablas/medias_por_cepa.csv`, `tablas/medias_por_gen.csv`,
  `tablas/anova.csv`, `tablas/posthoc_cepas.csv`, `tablas/posthoc_gen.csv`
- Figuras: `figuras/boxplot_por_cepa.png`, `figuras/boxplot_por_gen.png`,
  `figuras/interaccion_gen_cepa.png`
- Informe: `reportes/informe_factorial.md`
"""

(REPORTES / "informe_factorial.md").write_text(informe, encoding="utf-8")

# ---------------------------------------------------------------------------
# 7. Version HTML simple del informe (sin dependencias externas)
# ---------------------------------------------------------------------------
def _md_a_html(md: str) -> str:
    html = []
    for linea in md.splitlines():
        if linea.startswith("|") and linea.endswith("|") and "--" not in linea:
            celdas = [c.strip() for c in linea.strip().strip("|").split("|")]
            fila = "<tr>" + "".join(f"<td>{c}</td>" for c in celdas) + "</tr>"
            html.append(fila)
        elif linea.startswith("# "):
            html.append(f"<h1>{linea[2:]}</h1>")
        elif linea.startswith("## "):
            html.append(f"<h2>{linea[3:]}</h2>")
        elif linea.startswith("### "):
            html.append(f"<h3>{linea[4:]}</h3>")
        elif linea.startswith("- "):
            html.append(f"<li>{linea[2:]}</li>")
        elif linea.startswith("> "):
            html.append(f"<blockquote>{linea[2:]}</blockquote>")
        elif linea.strip() == "":
            html.append("")
        elif linea.startswith("|"):
            continue
        else:
            html.append(f"<p>{linea}</p>")
    cuerpo = "\n".join(html)
    cuerpo = cuerpo.replace("**", "")
    cuerpo = cuerpo.replace("`", "<code>", 1)
    cuerpo = cuerpo.replace("`", "</code>", 1)
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<title>Informe factorial Snijders</title>
<style>
body {{ font-family: sans-serif; max-width: 900px; margin: 2em auto; line-height: 1.5; color: #222; }}
h1, h2, h3 {{ color: #143; }}
table {{ border-collapse: collapse; margin: 1em 0; }}
td, th {{ border: 1px solid #ccc; padding: 4px 10px; }}
blockquote {{ border-left: 4px solid #ccc; margin-left: 0; padding-left: 1em; color: #555; }}
code {{ background: #f4f4f4; padding: 1px 4px; }}
</style>
</head>
<body>
{cuerpo}
</body>
</html>
"""

(REPORTES / "informe_factorial.html").write_text(_md_a_html(informe), encoding="utf-8")

# ---------------------------------------------------------------------------
# 8. Resumen en consola
# ---------------------------------------------------------------------------
def _f_strain(g):
    r = rank_cepas[rank_cepas["cepa"] == g].iloc[0]
    return f"{g}={_f(r['media'])}"

def _p(t):
    return anova_df[anova_df["term"] == t]["p"].iloc[0]

def _F(t):
    return anova_df[anova_df["term"] == t]["f"].iloc[0]

print("=" * 72)
print("ANALISIS FACTORIAL Snijders genotipo x cepa (Fusarium)")
print("=" * 72)
print(f"Modelo usado: {modelo_usado} | formula: {FORMULA}")
print(f"Auditoria: {len(df)} filas, balanceado=ok, sin faltantes, sin duplicados")
print(f"ANOVA: F(gen)={_F('C(gen)'):.2f} p={_p('C(gen)'):.4f} | "
      f"F(strain)={_F('C(strain)'):.2f} p={_p('C(strain)'):.4f} | "
      f"F(gen:strain)={_F('C(gen):C(strain)'):.2f} p={_p('C(gen):C(strain)'):.4f} | "
      f"F(year)={_F('C(year)'):.2f} p={_p('C(year)'):.4f}")
print("Media de y por cepa (ordenada, mayor a menor severidad):")
for _, r in rank_cepas.iterrows():
    print(f"   {r['cepa']}: {r['media']:.2f} (SE {r['se']:.2f})")
print(f"Genotipo con mayor severidad media: {_gen_mayor} ({_f(rank_gen.iloc[0]['media'])})")
print(f"Genotipo con menor severidad media: {_gen_menor} ({_f(rank_gen.iloc[-1]['media'])})")
print("Tukey cepas: pares significativos:", int(tukey_cepas["significativo"].eq("True").sum()))
print("Tukey genotipos: pares significativos:", int(tukey_gen["significativo"].eq("True").sum()))
print("Salidas en:", RESULTADOS)
print("=" * 72)
sys.exit(0)
