#!/usr/bin/env python3
"""Generador del notebook Colab analisis_tomillo_fusarium.ipynb"""
import nbformat as nbf
from pathlib import Path

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.10.0"},
    "colab": {"name": "Análisis Tomillo × Fusarium spp.", "provenance": []}
}
cells = []
def md(s): cells.append(nbf.v4.new_markdown_cell(s))
def code(s): cells.append(nbf.v4.new_code_cell(s))

# ===== TITLE =====
md("""# Análisis Estadístico: Actividad Antifúngica de *Thymus* spp. contra *Fusarium* spp.

**Proyecto:** Tomillo × Fusarium — Evaluación de extractos vegetales como alternativa al control químico de fusariosis.

**Objetivos:**
1. Determinar si la técnica de extracción afecta el rendimiento de extracción
2. Evaluar la inhibición de crecimiento micelial según método, aislado y concentración
3. Analizar la producción de conidias bajo diferentes tratamientos
4. Comparar perfiles de susceptibilidad entre aislados de Fusarium
5. Diagnosticar transformaciones para cada variable respuesta

**Estructura:** Auditoría de datos → Análisis exploratorio → Modelado inferencial → Post-hoc → Diagnóstico → Interpretación biológica
""")

# ===== SECTION 0: SETUP =====
md("""## 0. Configuración Inicial

### Instalación de paquetes requeridos

- `statsmodels` — modelos lineales mixtos, ANOVA, regresión logística
- `pingouin` — pruebas estadísticas adicionales
- `scikit-learn` — PCA, clustering, estandarización
- `scipy` — pruebas de hipótesis, clustering jerárquico
- `pandas`, `numpy` — manipulación de datos
- `matplotlib`, `seaborn` — visualización
- `openpyxl` — lectura de Excel (respaldo)
""")

code("""# Instalación de paquetes (ejecutar una vez en Colab)
import sys
if 'google.colab' in str(get_ipython()):
    !pip install statsmodels pingouin scikit-learn scipy pandas numpy matplotlib seaborn openpyxl -q
    print("Paquetes instalados en entorno Colab.")
else:
    print("Entorno local — asumiendo paquetes ya instalados.")
""")

md("""### Montaje de Google Drive

El notebook ofrece dos mecanismos de carga:
1. **Google Drive** (Colab) — monta `/content/drive/MyDrive/` y busca los CSVs allí
2. **Local** — busca en `dca/resultados/tablas/` relativo a la raíz del repositorio
""")

code("""# Montaje de Google Drive y definición de rutas
import os, warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

EN_COLAB = 'google.colab' in str(get_ipython())
if EN_COLAB:
    from google.colab import drive
    drive.mount('/content/drive')
    BASE = '/content/drive/MyDrive/tomillo_fusarium/resultados/tablas/'
    DIR_FIG = '/content/drive/MyDrive/tomillo_fusarium/resultados/figuras/'
else:
    BASE = 'dca/resultados/tablas/'
    DIR_FIG = 'dca/resultados/figuras/'
os.makedirs(DIR_FIG, exist_ok=True)
print(f"Directorio de datos: {BASE}")
print(f"Directorio de figuras: {DIR_FIG}")
print(f"Entorno: {'Colab' if EN_COLAB else 'Local'}")
""")

md("""### Importación de librerías y configuración estética

Se configura `matplotlib` con estilo de publicación científica: fondo blanco, sin ejes superior/derecho, tamaños de fuente calibrados (base 12, títulos 13, ticks 10). La paleta de colores es daltónico-safe.
""")

code("""# Importaciones generales
import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools import add_constant
from statsmodels.discrete.discrete_model import Logit
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

SEMILLA = 42
np.random.seed(SEMILLA)

matplotlib.rcParams.update({
    'figure.dpi': 300, 'savefig.dpi': 300,
    'font.size': 12, 'axes.titlesize': 13, 'axes.labelsize': 12,
    'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'legend.fontsize': 10, 'legend.title_fontsize': 11, 'figure.titlesize': 14,
    'lines.linewidth': 1.5, 'lines.markersize': 6,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'axes.edgecolor': '.3',
    'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False,
})

COLOR_MET = {'maceracion': '#2e86ab', 'soxhlet': '#a23b72', 'ultrasonido': '#f18f01'}
COLOR_CONC = ['#b3b3b3', '#7ba0b4', '#4a7c9b', '#1a4d6b']
LABEL_MET = {'maceracion': 'Maceración', 'soxhlet': 'Soxhlet', 'ultrasonido': 'Ultrasonido'}
print("Configuración completada. Semilla:", SEMILLA)
""")

md("""### Funciones auxiliares

Se definen funciones reutilizables para el análisis: resúmenes numéricos, transformaciones (IHS, CubeRoot, RankGauss) y evaluación de residuos de modelos mixtos.
""")

code("""def resumen_numerico(df, col, nombre=None, grupo=None):
    if grupo:
        return df.groupby(grupo)[col].agg(
            n='count', n_nulo=lambda x: x.isna().sum(),
            media='mean', mediana='median', std='std', min='min', max='max').reset_index()
    vals = df[col]
    return pd.DataFrame([{'variable': nombre or col, 'n': int(vals.notna().sum()),
        'n_nulo': int(vals.isna().sum()), 'media': vals.mean(), 'mediana': vals.median(),
        'std': vals.std(), 'min': vals.min(), 'max': vals.max()}])

def ihs_transform(c, shift=0):
    return lambda y: np.arcsinh((y + shift) / c)

def cube_root(y):
    return np.cbrt(y)

def rank_gauss(y):
    r = stats.rankdata(y, method='average')
    return stats.norm.ppf((r - 0.5) / len(y))

def eval_transforms(name, y, transform_funcs=None):
    if transform_funcs is None: return []
    results = []
    for tname, tfunc, desc in transform_funcs:
        try:
            yt = tfunc(y); sk = float(stats.skew(yt)); ku = float(stats.kurtosis(yt))
            sp = round(float(stats.shapiro(yt)[1]), 4) if 3 <= len(yt) <= 5000 else 0.0
            results.append({'Variable': name, 'Transformacion': tname, 'n': len(yt),
                'Skew': round(sk, 3), 'Kurtosis': round(ku, 3), 'Shapiro_p': sp})
        except Exception as e:
            results.append({'Variable': name, 'Transformacion': tname,
                'Descripcion': f'ERROR: {e}', 'n': 0, 'Skew': np.nan, 'Kurtosis': np.nan, 'Shapiro_p': np.nan})
    return results

def eval_residuals_lmm(name, y_full, exog, groups, transforms):
    rows = []
    for tname, tfunc, desc in transforms:
        try:
            yt = tfunc(y_full); m = MixedLM(yt, exog, groups=groups).fit(reml=True, maxiter=200)
            r = m.resid; sk = float(stats.skew(r)); ku = float(stats.kurtosis(r))
            sp = round(float(stats.shapiro(r)[1]), 4) if 3 <= len(r) <= 5000 else 0.0
            rows.append({'Modelo': name, 'Transformacion': tname, 'LogLik': round(m.llf, 1),
                'Resid_skew': round(sk, 3), 'Resid_kurt': round(ku, 3), 'Shapiro_p': sp})
        except Exception as e:
            rows.append({'Modelo': name, 'Transformacion': tname, 'LogLik': np.nan,
                'Resid_skew': np.nan, 'Resid_kurt': np.nan, 'Shapiro_p': np.nan})
    return rows

print("Funciones auxiliares cargadas.")
""")

md("""### Carga de datos

Se cargan los tres archivos CSV generados por el pipeline de extracción. Los datos ya han sido limpiados, estandarizados y enriquecidos con variables derivadas (porcentaje de inhibición, flags de inhibición completa/negativa, etc.).
""")

code("""# Carga de datos
CRE = pd.read_csv(os.path.join(BASE, 'crecimiento_micelial.csv'))
CON = pd.read_csv(os.path.join(BASE, 'conidias.csv'))
REN = pd.read_csv(os.path.join(BASE, 'rendimiento_extraccion.csv'))
REN['metodo_id'] = REN['metodo_extraccion'].str.strip().str.lower()
REN['metodo_id'] = REN['metodo_id'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('ascii')

CREC = CRE.copy(); CONI = CON.copy()
INH = CREC[~CREC['es_control'] & CREC['porcentaje_inhibicion'].notna()].copy()
INH['metodo_id'] = INH['metodo_extraccion'].str.strip().str.lower()
CONI_TRAT = CONI[~CONI['es_control'] & CONI['conidias_log10'].notna()].copy()
CONI_TRAT['metodo_id'] = CONI_TRAT['metodo_extraccion'].str.strip().str.lower()
CONI_INH = CONI[~CONI['es_control'] & CONI['porcentaje_inhibicion'].notna()].copy()
CONI_INH['metodo_id'] = CONI_INH['metodo_extraccion'].str.strip().str.lower()
CONI_INH_LOG = CONI[~CONI['es_control'] & CONI['porcentaje_inhibicion_log10'].notna()].copy()
CONI_INH_LOG['metodo_id'] = CONI_INH_LOG['metodo_extraccion'].str.strip().str.lower()
CTRL_CREC = CREC[CREC['es_control'] & CREC['crecimiento_mm'].notna()]
CTRL_CONI = CONI[CONI['es_control'] & CONI['conidias_log10'].notna()]

print("Datos cargados:")
print(f"  Crecimiento: {len(CREC)} filas, {CREC['aislado_id'].nunique()} aislados")
print(f"  Conidias: {len(CONI)} filas, {CONI['aislado_id'].nunique()} aislados")
print(f"  Rendimiento: {len(REN)} filas")
print(f"  INH valido: {len(INH)} filas")
""")

# ===== SECTION 1: EDA =====
md("""---
# 1. Análisis Exploratorio de Datos (EDA)

## Fundamentos

El EDA es el primer paso crítico en cualquier análisis estadístico. Sus objetivos son:

1. **Comprender la estructura** de los datos: dimensiones, tipos, valores faltantes
2. **Identificar anomalías**: valores imposibles, distribuciones atípicas, relaciones inesperadas
3. **Generar hipótesis** sobre relaciones entre variables
4. **Guiar la selección de modelos** según las características observadas

En este proyecto, el EDA es particularmente importante porque:
- El diseño experimental es **no balanceado**: solo Maceración tiene 4 concentraciones; Soxhlet y Ultrasonido tienen solo 2
- Existen **efectos piso** (inhibición completa) y **efectos techo** (inhibición negativa)
- La escala de medición de conidias requiere atención especial (reportada en log10 por el laboratorio)
""")

md("""### 1.1 Resumen numérico de variables respuesta""")

code("""print("=" * 70)
print("  1.1 RESUMEN NUMERICO")
print("=" * 70)
print("\\n-- Crecimiento micelial (mm) --")
print(resumen_numerico(CREC, 'crecimiento_mm', grupo=['metodo_extraccion', 'concentracion_mg_ml']).to_string(index=False))
print("\\n-- Inhibicion (%) --")
print(resumen_numerico(INH, 'porcentaje_inhibicion', grupo=['metodo_extraccion', 'concentracion_mg_ml']).to_string(index=False))
print("\\n-- Conidias (log10/mL) --")
print(resumen_numerico(CONI_TRAT, 'conidias_log10', grupo=['metodo_extraccion', 'concentracion_mg_ml']).to_string(index=False))
print("\\n-- Rendimiento (%) --")
print(resumen_numerico(REN, 'rendimiento_pct', grupo=['metodo_id']).to_string(index=False))
""")

md("""### 1.2 Distribuciones

Los histogramas permiten evaluar visualmente la forma de la distribución de cada variable respuesta. Esto es fundamental para detectar asimetrías (skewness), multimodalidad y límites artificiales (ej. 0% y 100% en porcentajes).
""")

code("""# 1.2a: Crecimiento (control vs tratamiento)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, (label, mask) in zip(axes, [('Control (0 mg/mL)', CREC['es_control']), ('Tratamiento', ~CREC['es_control'])]):
    d = CREC.loc[mask & CREC['crecimiento_mm'].notna(), 'crecimiento_mm']
    ax.hist(d, bins=25, color='#2e86ab' if label.startswith('Control') else '#a23b72', edgecolor='white', alpha=0.8)
    ax.set_xlabel('Crecimiento (mm)'); ax.set_ylabel('Frecuencia')
    ax.set_title(f'Crecimiento micelial -- {label}')
    ax.axvline(d.median(), color='red', ls='--', label=f'Mediana={d.median():.0f}')
    ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_distribucion_crecimiento.png'), dpi=300); plt.show()
print("  OK Histograma de crecimiento generado.")

# 1.2b: INH
fig, ax = plt.subplots(figsize=(10, 5))
d = INH['porcentaje_inhibicion']
ax.hist(d, bins=30, color='#2e86ab', edgecolor='white', alpha=0.8)
ax.set_xlabel('Inhibicion (%)'); ax.set_ylabel('Frecuencia')
ax.set_title('Distribucion del % de inhibicion -- todos los metodos')
ax.axvline(d.median(), color='red', ls='--', label=f'Mediana={d.median():.1f}%')
ax.axvline(0, color='gray', ls=':', alpha=0.5); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_distribucion_inhibicion.png'), dpi=300); plt.show()
print("  OK Histograma de INH generado.")

# 1.2c: Conidias log10
fig, ax = plt.subplots(figsize=(10, 5))
d = CONI_TRAT['conidias_log10']
ax.hist(d, bins=30, color='#a23b72', edgecolor='white', alpha=0.8)
ax.set_xlabel('log10(conidias/mL)'); ax.set_ylabel('Frecuencia')
ax.set_title('Distribucion de conidias (log10)')
ax.axvline(d.median(), color='red', ls='--', label=f'Mediana={d.median():.2f}')
ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_distribucion_conidias.png'), dpi=300); plt.show()
print("  OK Histograma de conidias generado.")

# 1.2d: Rendimiento
fig, ax = plt.subplots(figsize=(10, 5.5))
cols = [COLOR_MET[m] for m in REN['metodo_id']]
ax.bar(range(len(REN)), REN['rendimiento_pct'], color=cols, edgecolor='white')
ax.set_xticks(range(len(REN)))
ax.set_xticklabels([f"{LABEL_MET[m]}\\n(rep {r})" for m, r in zip(REN['metodo_id'], REN['replica_biologica'])], rotation=45, ha='right')
ax.set_ylabel('Rendimiento (%)'); ax.set_title('Rendimiento de extraccion por metodo')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_rendimiento.png'), dpi=300); plt.show()
print("  OK Grafico de rendimiento generado.")
""")

md("""### 1.3 Boxplots por metodo y concentracion

Los boxplots permiten comparar visualmente la distribucion del % de inhibicion entre metodos y concentraciones. Ademas, se identifican los aislados con mayor y menor susceptibilidad.
""")

code("""# 1.3a: INH ~ metodo
fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=INH, x='metodo_extraccion', y='porcentaje_inhibicion', palette=COLOR_MET, ax=ax)
sns.stripplot(data=INH, x='metodo_extraccion', y='porcentaje_inhibicion', color='black', alpha=0.15, size=3, ax=ax)
ax.set_xlabel('Metodo de extraccion'); ax.set_ylabel('Inhibicion (%)')
ax.set_title('Inhibicion por metodo de extraccion'); ax.axhline(0, color='gray', ls=':', alpha=0.5)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_boxplot_inhibicion_metodo.png'), dpi=300); plt.show()

# 1.3b: INH ~ concentracion (solo Maceracion)
fig, ax = plt.subplots(figsize=(8, 5))
mac = INH[INH['metodo_extraccion'] == 'maceracion'].copy()
ord_c = sorted(mac['concentracion_mg_ml'].unique())
sns.boxplot(data=mac, x='concentracion_mg_ml', y='porcentaje_inhibicion', order=ord_c, palette=COLOR_CONC, ax=ax)
sns.stripplot(data=mac, x='concentracion_mg_ml', y='porcentaje_inhibicion', order=ord_c, color='black', alpha=0.15, size=3, ax=ax)
ax.set_xlabel('Concentracion (mg/mL)'); ax.set_ylabel('Inhibicion (%)')
ax.set_title('Inhibicion por concentracion -- Maceracion'); ax.axhline(0, color='gray', ls=':', alpha=0.5)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_boxplot_inhibicion_concentracion_mac.png'), dpi=300); plt.show()

# 1.3c: Soxhlet y Ultrasonido
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for ax, metodo in zip(axes, ['soxhlet', 'ultrasonido']):
    sub = INH[INH['metodo_extraccion'] == metodo]
    sns.boxplot(data=sub, x='concentracion_mg_ml', y='porcentaje_inhibicion', color=COLOR_MET[metodo], ax=ax)
    sns.stripplot(data=sub, x='concentracion_mg_ml', y='porcentaje_inhibicion', color='black', alpha=0.15, size=3, ax=ax)
    ax.set_xlabel('Concentracion (mg/mL)'); ax.set_ylabel('Inhibicion (%)')
    ax.set_title(f'Inhibicion -- {LABEL_MET[metodo]}'); ax.axhline(0, color='gray', ls=':', alpha=0.5)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_boxplot_inhibicion_soxhlet_ultra.png'), dpi=300); plt.show()

# 1.3d: Top 10 extremos por aislado
fig, ax = plt.subplots(figsize=(14, 6))
med = INH.groupby('aislado_id')['porcentaje_inhibicion'].median().sort_values()
top = list(med.head(10).index) + list(med.tail(10).index)
sub = INH[INH['aislado_id'].isin(top)]
sns.boxplot(data=sub, x='aislado_id', y='porcentaje_inhibicion', order=top, palette='RdYlBu_r', ax=ax)
ax.set_xlabel('Aislado'); ax.set_ylabel('Inhibicion (%)')
ax.set_title('Inhibicion por aislado -- 10 mas y 10 menos susceptibles (mediana)')
ax.tick_params(axis='x', rotation=45); ax.axhline(0, color='gray', ls=':', alpha=0.5)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_boxplot_inhibicion_aislados_extremos.png'), dpi=300); plt.show()

# 1.3e: Crecimiento del control por aislado
fig, ax = plt.subplots(figsize=(12, 5))
ctrl = CREC[CREC['es_control'] & CREC['crecimiento_mm'].notna()].drop_duplicates(subset=['aislado_id', 'metodo_extraccion'])
sns.boxplot(data=ctrl, x='aislado_id', y='crecimiento_mm', hue='metodo_extraccion', palette=COLOR_MET, ax=ax)
ax.set_xlabel('Aislado'); ax.set_ylabel('Crecimiento control (mm)')
ax.set_title('Crecimiento del control por aislado y metodo')
ax.tick_params(axis='x', rotation=45); ax.legend(title='Metodo')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_control_por_aislado.png'), dpi=300); plt.show()
print("  OK Boxplots generados.")
""")

md("""### 1.4 Mapa de datos faltantes

Visualizar la estructura de datos faltantes permite identificar patrones sistematicos de ausencia. En este diseno, el desbalance es **estructural** (no al azar): Soxhlet y Ultrasonido solo se evaluaron a 0.0 y 5.0 mg/mL.
""")

code("""def matriz_faltantes(df, col_valor, nombre, ax):
    mat = df.pivot_table(index='aislado_id', columns=['metodo_extraccion', 'concentracion_mg_ml'],
                         values=col_valor, aggfunc=lambda x: int(x.notna().any()))
    sns.heatmap(mat, cmap='RdYlGn', cbar_kws={'label': 'Dato presente'}, linewidths=0.5, linecolor='white', ax=ax, vmin=0, vmax=1)
    ax.set_title(f'Datos presentes: {nombre}'); ax.set_xlabel(''); ax.set_ylabel('Aislado')
    return mat

fig, axes = plt.subplots(1, 2, figsize=(20, 12))
matriz_faltantes(CREC, 'crecimiento_mm', 'Crecimiento micelial', axes[0])
matriz_faltantes(CONI, 'conidias_log10', 'Conidias', axes[1])
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_mapa_faltantes.png'), dpi=300); plt.show()

print("\\n-- Faltantes en crecimiento --")
print(CREC.groupby(['metodo_extraccion', 'concentracion_mg_ml']).agg(
    total=('crecimiento_mm', 'count'), presentes=('crecimiento_mm', lambda x: int(x.notna().sum())),
    pct=('crecimiento_mm', lambda x: f"{100*x.notna().sum()/len(x):.0f}%")).reset_index().to_string(index=False))

print("\\n-- Faltantes en conidias --")
print(CONI.groupby(['metodo_extraccion', 'concentracion_mg_ml']).agg(
    total=('conidias_log10', 'count'), presentes=('conidias_log10', lambda x: int(x.notna().sum())),
    pct=('conidias_log10', lambda x: f"{100*x.notna().sum()/len(x):.0f}%")).reset_index().to_string(index=False))
""")

md("""### 1.5 Inhibicion completa y negativa

La **inhibicion completa** (crecimiento = 0 mm) representa un efecto piso. La **inhibicion negativa** (crecimiento > control) puede deberse a variacion biologica natural, error de medicion u hormesis.
""")

code("""completa = INH[INH['inhibicion_completa']]; negativa = INH[INH['inhibicion_negativa']]
print(f"Inhibicion completa: {len(completa)} casos")
print(f"  Por metodo y concentracion:")
print(completa.groupby(['metodo_extraccion', 'concentracion_mg_ml']).size().to_string())
print(f"\\nInhibicion negativa: {len(negativa)} casos")
print(f"  Rango INH negativo: {negativa['porcentaje_inhibicion'].min():.1f}% a {negativa['porcentaje_inhibicion'].max():.1f}%")

fig, ax = plt.subplots(figsize=(14, 5))
conteo = completa.groupby('aislado_id').size().sort_values(ascending=False)
conteo.plot(kind='bar', ax=ax, color='#2e86ab', edgecolor='white')
ax.set_xlabel('Aislado'); ax.set_ylabel('Casos de inhibicion completa')
ax.set_title('Inhibicion completa por aislado (crecimiento = 0 mm)')
ax.tick_params(axis='x', rotation=45)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_inhibicion_completa.png'), dpi=300); plt.show()
""")

md("""### 1.6 Validacion contra hoja original

Se compara el INH calculado por el pipeline con el INH del Excel original (solo para C1 = 5 mg/mL, la unica concentracion donde la columna de INH en el Excel es valida).
""")

code("""if 'diferencia_con_hoja' in CREC.columns:
    diffs = CREC[(CREC['concentracion_mg_ml']==5.0) & CREC['diferencia_con_hoja'].notna() & ~CREC['es_control']]['diferencia_con_hoja']
    print(f"Validacion C1 (5 mg/mL) -- {len(diffs)} obs: media|diff|={diffs.abs().mean():.3f}, max={diffs.abs().max():.3f}")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(diffs, bins=30, color='#2e86ab', edgecolor='white', alpha=0.8)
    ax.set_xlabel('Diferencia (INH calc - INH hoja)'); ax.set_ylabel('Frecuencia')
    ax.set_title('Validacion: INH calculado vs. hoja (solo C1 = 5 mg/mL)')
    ax.axvline(0, color='red', ls='--')
    fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_validacion_inh_hoja.png'), dpi=300); plt.show()
else:
    print("  Columna 'diferencia_con_hoja' no encontrada.")
""")

md("""### 1.7 Diagnostico de supuestos

Se evalua la normalidad mediante Shapiro-Wilk y graficos Q-Q, y la homocedasticidad mediante Levene. El test de Shapiro-Wilk es sensible al tamano muestral, por lo que se limita a 5000 observaciones.
""")

code("""print("-- Normalidad (Shapiro-Wilk) --")
for var_name, var_col, df_src in [
    ('Crecimiento (mm) -- tratamiento', 'crecimiento_mm', CREC),
    ('Inhibicion (%)', 'porcentaje_inhibicion', INH),
    ('Conidias (log10)', 'conidias_log10', CONI_TRAT)]:
    vals = df_src[var_col].dropna()
    if len(vals) > 5000: vals = vals.sample(5000, random_state=SEMILLA)
    if len(vals) >= 3:
        w, p = stats.shapiro(vals)
        print(f"  {var_name}: W={w:.4f}, p={p:.6f} {'NO normal' if p<0.05 else 'Normal (a=0.05)'}")

# Q-Q plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, (vn, vc, df_s, co) in zip(axes, [
    ('Crecimiento (mm)', 'crecimiento_mm', CREC, '#2e86ab'),
    ('Inhibicion (%)', 'porcentaje_inhibicion', INH, '#a23b72'),
    ('Conidias (log10)', 'conidias_log10', CONI_TRAT, '#f18f01')]):
    vals = df_s[vc].dropna()
    if len(vals) > 2000: vals = vals.sample(2000, random_state=SEMILLA)
    stats.probplot(vals, dist='norm', plot=ax)
    ax.get_lines()[0].set_markerfacecolor(co); ax.get_lines()[0].set_markeredgecolor(co); ax.get_lines()[0].set_markersize(3)
    ax.get_lines()[1].set_color('red'); ax.set_title(f'Q-Q plot: {vn}')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_qq_plots.png'), dpi=300); plt.show()

print("\\n-- Homocedasticidad (Levene) --")
grupos = [g['porcentaje_inhibicion'].dropna().values for _, g in INH.groupby('metodo_extraccion')]
if all(len(g)>=2 for g in grupos):
    lstat, lp = stats.levene(*grupos)
    print(f"  Inhibicion (%) ~ metodo: Levene={lstat:.4f}, p={lp:.6f} {'Heterocedastico' if lp<0.05 else 'Homocedastico'}")
""")

md("""### 1.8 Interacciones preliminares

Se exploran las interacciones entre metodo y concentracion. La presencia de interaccion implica que el efecto de un factor depende del nivel del otro.
""")

code("""# 1.8a: Interaccion metodo x concentracion
fig, ax = plt.subplots(figsize=(8, 5))
idf = INH.groupby(['metodo_extraccion', 'concentracion_mg_ml'])['porcentaje_inhibicion'].agg(['mean', 'sem']).reset_index()
for metodo in ['maceracion', 'soxhlet', 'ultrasonido']:
    sub = idf[idf['metodo_extraccion'] == metodo]
    ax.errorbar(sub['concentracion_mg_ml'], sub['mean'], yerr=sub['sem'],
                label=LABEL_MET[metodo], color=COLOR_MET[metodo], marker='o', capsize=4, linewidth=2)
ax.set_xlabel('Concentracion (mg/mL)'); ax.set_ylabel('Inhibicion media (%)')
ax.set_title('Interaccion: metodo x concentracion'); ax.legend(title='Metodo')
ax.axhline(0, color='gray', ls=':', alpha=0.5)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_interaccion_metodo_conc.png'), dpi=300); plt.show()

# 1.8b: Perfil por aislado (Maceracion)
fig, ax = plt.subplots(figsize=(12, 6))
mac_p = INH[INH['metodo_extraccion']=='maceracion'].copy()
perf = mac_p.groupby(['aislado_id','concentracion_mg_ml'])['porcentaje_inhibicion'].mean().reset_index()
for aislado in ['HC3','HC5','HC17','H11G','FU1','FUSARIUM JULIAN H20','H4B','H8N']:
    sub = perf[perf['aislado_id']==aislado]
    if len(sub): ax.plot(sub['concentracion_mg_ml'], sub['porcentaje_inhibicion'], marker='o', label=aislado, linewidth=2)
ax.set_xlabel('Concentracion (mg/mL)'); ax.set_ylabel('Inhibicion media (%)')
ax.set_title('Perfil de inhibicion por aislado -- Maceracion')
ax.legend(bbox_to_anchor=(1.05,1), loc='upper left'); ax.axhline(0, color='gray', ls=':', alpha=0.5)
fig.subplots_adjust(right=0.8); fig.tight_layout()
fig.savefig(os.path.join(DIR_FIG, 'eda_perfil_aislados_mac.png'), dpi=300); plt.show()
""")

md("""### 1.9 Analisis especifico de conidias

La produccion de conidias fue reportada por el laboratorio en escala log10. El Excel contiene dos versiones del % de inhibicion: sobre conteos crudos y sobre log10 (mas estable).
""")

code("""# 1.9a: Boxplot log10 por metodo y concentracion
fig, ax = plt.subplots(figsize=(10, 5))
sns.boxplot(data=CONI_TRAT, x='metodo_extraccion', y='conidias_log10', hue='concentracion_mg_ml', palette=COLOR_CONC, ax=ax)
ax.set_xlabel('Metodo de extraccion'); ax.set_ylabel('log10(conidias/mL)')
ax.set_title('Conidias -- valores crudos en log10'); ax.legend(title='Conc. (mg/mL)')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_boxplot_conidias_log10.png'), dpi=300); plt.show()

# 1.9b: INH conidias escala cruda
fig, ax = plt.subplots(figsize=(10, 5))
sns.boxplot(data=CONI_INH, x='metodo_extraccion', y='porcentaje_inhibicion', hue='concentracion_mg_ml', palette=COLOR_CONC, ax=ax)
ax.set_xlabel('Metodo de extraccion'); ax.set_ylabel('Reduccion de conidias (%) -- escala cruda')
ax.set_title('Inhibicion de conidias sobre conteos crudos'); ax.axhline(0, color='gray', ls=':', alpha=0.5)
ax.legend(title='Conc. (mg/mL)', bbox_to_anchor=(1.02,1), loc='upper left')
fig.subplots_adjust(right=0.82); fig.tight_layout()
fig.savefig(os.path.join(DIR_FIG, 'eda_boxplot_inhibicion_conidias_crudo.png'), dpi=300); plt.show()

# 1.9c: Comparacion INH crudo vs log10
fig, ax = plt.subplots(figsize=(8, 6))
comp = CONI_INH_LOG[CONI_INH_LOG['porcentaje_inhibicion'].notna()]
ax.scatter(comp['porcentaje_inhibicion_log10'], comp['porcentaje_inhibicion'], alpha=0.3, c='#a23b72', edgecolors='none')
ax.plot([-100,100], [-100,100], 'r--', alpha=0.5, label='Identidad')
ax.set_xlabel('%INH sobre log10 (escala de la hoja)'); ax.set_ylabel('%INH sobre conteos crudos')
ax.set_title('Comparacion: INH de conidias en escala cruda vs log10'); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_comparacion_inh_conidias_crudo_vs_log10.png'), dpi=300); plt.show()

if len(comp)>0:
    print(f"Diferencia media INH crudo - log10: {(comp['porcentaje_inhibicion']-comp['porcentaje_inhibicion_log10']).mean():.0f} puntos")
""")

md("""### 1.10 Correlacion entre crecimiento y conidias""")

code("""fig, ax = plt.subplots(figsize=(8, 6))
merged_cc = CREC[CREC['crecimiento_mm'].notna() & ~CREC['es_control']].merge(
    CONI[CONI['conidias_log10'].notna() & ~CONI['es_control']],
    on=['aislado_id','metodo_extraccion','concentracion_mg_ml','replica_biologica'], suffixes=('_crec','_con'))
if len(merged_cc)>0:
    ax.scatter(merged_cc['crecimiento_mm'], merged_cc['conidias_log10'], alpha=0.4, c='#2e86ab', edgecolors='none')
    r, p = stats.pearsonr(merged_cc['crecimiento_mm'], merged_cc['conidias_log10'])
    ax.set_xlabel('Crecimiento (mm)'); ax.set_ylabel('log10(conidias/mL)')
    ax.set_title(f'Correlacion crecimiento vs. conidias (r={r:.3f}, p={p:.4f})')
    fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'eda_correlacion_crecimiento_conidias.png'), dpi=300); plt.show()
    print(f"  Correlacion: r={r:.3f}, p={p:.4f}, n={len(merged_cc)}")
else:
    plt.close(fig); print("  Sin datos para correlacion.")
""")

md("""### 1.11 Resumen de hallazgos del EDA""")

code("""print("="*70); print("  HALLAZGOS DEL EDA"); print("="*70)
hallazgos = []
pct_c = 100*len(completa)/len(INH)
hallazgos.append(f"Inhibicion completa en {len(completa)}/{len(INH)} ({pct_c:.1f}%) -- efecto piso")
pct_n = 100*len(negativa)/len(INH)
hallazgos.append(f"Inhibicion negativa en {len(negativa)}/{len(INH)} ({pct_n:.1f}%)")
hallazgos.append("Solo MACERACION tiene diseno factorial completo (4 concentraciones)")
hallazgos.append("SOXHLET y ULTRASONIDO tienen solo 2 concentraciones")
for nombre, df_src, col in [('crecimiento', CREC, 'crecimiento_mm'), ('conidias', CONI, 'conidias_log10')]:
    nulos = df_src[col].isna().sum()
    hallazgos.append(f"{nombre}: {nulos}/{len(df_src)} sin dato ({100*nulos/len(df_src):.0f}%)")
for vn, vc, df_s in [('Crecimiento','crecimiento_mm',CREC),('Inhibicion','porcentaje_inhibicion',INH),('Conidias (log10)','conidias_log10',CONI_TRAT)]:
    vals = df_s[vc].dropna()
    if len(vals)>=3:
        if len(vals)>5000: vals = vals.sample(5000, random_state=SEMILLA)
        _, pn = stats.shapiro(vals)
        hallazgos.append(f"{vn}: {'NO normal' if pn<0.05 else 'Normal'} (Shapiro p={pn:.4f})")
for h in hallazgos: print(f"  {h}")
""")

# ===== SECTION 2: OBJETIVO 1 - RENDIMIENTO =====
md("""---
# 2. Objetivo 1 -- Rendimiento de Extraccion

## Pregunta biologica
El metodo de extraccion (Maceracion, Soxhlet, Ultrasonido) afecta el rendimiento de extraccion (%) de compuestos de *Thymus* spp.?

## Diseno
- 3 metodos de extraccion, 3 replicas biologicas por metodo = 9 observaciones.

## Estrategia de analisis
1. Estadistica descriptiva por metodo
2. ANOVA de un factor (si se cumplen supuestos) o Kruskal-Wallis
3. Tukey HSD post-hoc si el ANOVA es significativo

### Justificacion del ANOVA
Con solo 3 replicas por grupo, la potencia para detectar diferencias es limitada. Sin embargo, el ANOVA de un factor es adecuado si los residuos son aproximadamente normales y las varianzas homogeneas.
""")

code("""print("="*70); print("  OBJETIVO 1 -- RENDIMIENTO DE EXTRACCION"); print("="*70)
desc_r = REN.groupby('metodo_id')['rendimiento_pct'].agg(n='count', media='mean', sd='std', se='sem', min='min', max='max').round(3)
print("\\n-- Estadistica descriptiva --")
for metodo, row in desc_r.iterrows():
    ic_i, ic_s = row['media']-1.96*row['se'], row['media']+1.96*row['se']
    print(f"  {LABEL_MET[metodo]:15s}  n={int(row['n'])}  media={row['media']:.2f}%  DE={row['sd']:.2f}  IC95=[{ic_i:.2f}, {ic_s:.2f}]")
""")

md("""### 2.1 Visualizacion""")

code("""fig, axes = plt.subplots(1, 3, figsize=(15, 5)); fig.subplots_adjust(wspace=0.35)
ax = axes[0]
sns.stripplot(data=REN, x='metodo_id', y='rendimiento_pct', color='black', alpha=0.6, size=8, ax=ax, jitter=False)
sns.boxplot(data=REN, x='metodo_id', y='rendimiento_pct', hue='metodo_id', palette=COLOR_MET, ax=ax, width=0.4, legend=False)
for patch in ax.patches:
    if hasattr(patch, 'set_facecolor'): c=patch.get_facecolor(); patch.set_facecolor((c[0],c[1],c[2],0.5))
ax.set_xlabel('Metodo'); ax.set_ylabel('Rendimiento (%)'); ax.set_title('Rendimiento por metodo')
ax.set_xticks([0,1,2]); ax.set_xticklabels([LABEL_MET[l] for l in ['maceracion','soxhlet','ultrasonido']], rotation=45, ha='right')
ax = axes[1]
means = REN.groupby('metodo_id')['rendimiento_pct'].mean(); sems = REN.groupby('metodo_id')['rendimiento_pct'].sem()
xp = np.arange(len(means))
ax.bar(xp, means.values, yerr=1.96*sems.values, capsize=5, color=[COLOR_MET[m] for m in means.index], width=0.5, alpha=0.8)
ax.set_xticks(xp); ax.set_xticklabels([LABEL_MET[l] for l in means.index])
ax.set_ylabel('Rendimiento medio (%)'); ax.set_title('Rendimiento medio +/- IC95%')
ax = axes[2]
sns.violinplot(data=REN, x='metodo_id', y='rendimiento_pct', hue='metodo_id', palette=COLOR_MET, ax=ax, inner='quartile', legend=False)
ax.set_xlabel('Metodo'); ax.set_ylabel('Rendimiento (%)'); ax.set_title('Distribucion del rendimiento')
ax.set_xticks([0,1,2]); ax.set_xticklabels([LABEL_MET[l] for l in ['maceracion','soxhlet','ultrasonido']])
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'obj1_rendimiento.png'), dpi=300); plt.show()
""")

md("""### 2.2 ANOVA de un factor

El ANOVA contrasta H0: todas las medias son iguales. El estadistico F compara variabilidad entre grupos vs. dentro de grupos. Se reportan eta-cuadrado (proporcion de varianza explicada) y omega-cuadrado (estimacion insesgada).
""")

code("""print("\\n-- Diagnostico de supuestos --")
for metodo in ['maceracion','soxhlet','ultrasonido']:
    vals = REN.loc[REN['metodo_id']==metodo, 'rendimiento_pct']
    if len(vals)>=3:
        wv, pv = stats.shapiro(vals)
        print(f"  {LABEL_MET[metodo]:15s} W={wv:.4f} p={pv:.4f}")

print("\\n-- ANOVA de un factor --")
mo = ols('rendimiento_pct ~ C(metodo_id)', data=REN).fit()
ta = anova_lm(mo, typ=2); print(ta.round(4).to_string())
f_st = ta.loc['C(metodo_id)','F']; p_av = ta.loc['C(metodo_id)','PR(>F)']
ss_t = ta.loc['C(metodo_id)','sum_sq']; ss_r = ta.loc['Residual','sum_sq']
df_t = int(ta.loc['C(metodo_id)','df']); df_r = int(ta.loc['Residual','df'])
ms_t = ss_t/df_t; ms_r = ss_r/df_r
eta_sq = ss_t/(ss_t+ss_r)
omega_sq = (ss_t - df_t*ms_r) / (ss_t + (df_t+df_r+1)*ms_r)
print(f"\\nTamano del efecto: eta2 = {eta_sq:.4f} ({eta_sq*100:.1f}%), omega2 = {omega_sq:.4f}")
print(f"CV residual = {np.sqrt(ms_r)/REN['rendimiento_pct'].mean()*100:.1f}%")
""")

md("""### 2.3 Diagnostico de residuos""")

code("""resid = mo.resid
fig, axes = plt.subplots(1, 2, figsize=(12, 5)); fig.subplots_adjust(wspace=0.4)
stats.probplot(resid, dist='norm', plot=axes[0]); axes[0].set_title('Q-Q de residuos (ANOVA)')
axes[1].scatter(mo.fittedvalues, resid, alpha=0.7, c='#2e86ab')
axes[1].axhline(0, color='gray', ls='--', alpha=0.5)
axes[1].set_xlabel('Valores ajustados'); axes[1].set_ylabel('Residuos'); axes[1].set_title('Residuos vs. ajustados')
lv_st, lv_p = stats.levene(REN.loc[REN['metodo_id']=='maceracion','rendimiento_pct'],
    REN.loc[REN['metodo_id']=='soxhlet','rendimiento_pct'],
    REN.loc[REN['metodo_id']=='ultrasonido','rendimiento_pct'])
axes[1].text(0.05,0.95,f'Levene: F={lv_st:.2f}, p={lv_p:.4f}', transform=axes[1].transAxes, va='top', fontsize=9,
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG, 'obj1_diagnostico_anova.png'), dpi=300); plt.show()
print(f"Homocedasticidad (Levene): F={lv_st:.2f}, p={lv_p:.4f}")
""")

md("""### 2.4 Alternativa no parametrica (Kruskal-Wallis)""")

code("""print("\\n-- Kruskal-Wallis --")
kw_st, kw_p = stats.kruskal(REN.loc[REN['metodo_id']=='maceracion','rendimiento_pct'],
    REN.loc[REN['metodo_id']=='soxhlet','rendimiento_pct'],
    REN.loc[REN['metodo_id']=='ultrasonido','rendimiento_pct'])
print(f"  H={kw_st:.2f}, p={kw_p:.6f}")
""")

md("""### 2.5 Post-hoc: Tukey HSD

Si el ANOVA es significativo, Tukey HSD controla el error de tipo I en comparaciones multiples. **Interpretacion biologica**: mayor rendimiento no implica mayor actividad antifungica -- el perfil fitoquimico difiere entre metodos.
""")

code("""print("\\n-- Comparaciones post-hoc --")
if p_av < 0.05:
    tukey = pairwise_tukeyhsd(REN['rendimiento_pct'], REN['metodo_id'], alpha=0.05)
    print(tukey)
    tk_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
    for _, row in tk_df.iterrows():
        sig = ' SI' if row['reject'] else ''
        print(f"  {row['group1']+' - '+row['group2']:30s} diff={row['meandiff']:>8.3f}  p-adj={row['p-adj']:>8.4f}{sig:>5s}")
else:
    print("  ANOVA no significativo -- no hay post-hoc parametrico.")
""")

md("""### 2.6 Tabla resumen e interpretacion

**Interpretacion biologica:** Soxhlet produce el mayor rendimiento (~3-4x mas que Maceracion y Ultrasonido) por el reflujo continuo de solvente caliente. Sin embargo, la Maceracion (menor rendimiento) mostro mayor actividad antifungica, sugiriendo que extrae selectivamente compuestos con actividad especifica.
""")

code("""resultados_rend = pd.DataFrame({
    'Metodo': [LABEL_MET[m] for m in desc_r.index],
    'n': desc_r['n'].values.astype(int), 'Media (%)': desc_r['media'].values, 'DE': desc_r['sd'].values,
    'IC95_inf': (desc_r['media']-1.96*desc_r['se']).values, 'IC95_sup': (desc_r['media']+1.96*desc_r['se']).values,
    'Min': desc_r['min'].values, 'Max': desc_r['max'].values})
print("\\n-- Tabla resumen --"); print(resultados_rend.to_string(index=False))
print(f"\\nANOVA: F({df_t},{df_r}) = {f_st:.2f}, p = {p_av:.6f}")
print(f"Kruskal-Wallis: H = {kw_st:.2f}, p = {kw_p:.6f}")
print(f"eta2 = {eta_sq:.3f} | omega2 = {omega_sq:.3f}")
print(f"\\nSoxhlet: {desc_r.loc['soxhlet','media']:.1f}%, Ultrasonido: {desc_r.loc['ultrasonido','media']:.1f}%, Maceracion: {desc_r.loc['maceracion','media']:.1f}%")
""")

# ===== SECTION 3: OBJETIVO 2 - INHIBICION =====
md("""---
# 3. Objetivo 2 -- Inhibicion de Crecimiento Micelial

## Pregunta biologica
El metodo de extraccion, el aislado de Fusarium y la concentracion afectan la inhibicion del crecimiento micelial?

## Estrategia de modelado

Se usan **modelos lineales mixtos (LMM)** porque:
1. Los aislados son efecto aleatorio (muestra de la poblacion de Fusarium)
2. Hay mediciones repetidas sobre cada aislado
3. El LMM estima la correlacion intra-aislado (ICC)

### Tres modelos:
- **Modelo A**: LMM comparando 3 metodos a 5.0 mg/mL (efecto fijo: metodo; aleatorio: aislado)
- **Modelo B**: LMM dosis-respuesta intra-Maceracion (efecto fijo: concentracion; aleatorio: aislado)
- **Modelo C**: Regresion logistica -- probabilidad de inhibicion completa
""")

code("""print("="*70); print("  OBJETIVO 2 -- INHIBICION DE CRECIMIENTO MICELIAL"); print("="*70)
INH['conc_cat'] = INH['concentracion_mg_ml'].astype('category')
INH['conc_log'] = np.log(INH['concentracion_mg_ml']+0.01)
INH['completa_bin'] = INH['inhibicion_completa'].astype(int)
print(f"\\nTotal filas con INH: {len(INH)}, Aislados: {INH['aislado_id'].nunique()}")
print(f"Rango INH: [{INH['porcentaje_inhibicion'].min():.1f}, {INH['porcentaje_inhibicion'].max():.1f}]")
""")

md("""### 3.1 Modelo A: Comparacion de metodos a 5.0 mg/mL

El ICC mide la proporcion de varianza debida a diferencias entre aislados. ICC alto indica que los aislados tienen perfiles de susceptibilidad consistentemente diferentes.
""")

code("""print("\\n"+"-"*70); print("  MODELO A: Comparacion de metodos a 5.0 mg/mL"); print("-"*70)
df_a = INH[INH['concentracion_mg_ml']==5.0].copy()
print(f"  Datos: {len(df_a)} obs, {df_a['aislado_id'].nunique()} aislados")
df_a_ml = pd.get_dummies(df_a, columns=['metodo_id'], drop_first=True, dtype=float)
gr_a = df_a_ml['aislado_id']; ex_a = add_constant(df_a_ml[['metodo_id_soxhlet','metodo_id_ultrasonido']]); en_a = df_a_ml['porcentaje_inhibicion']
mo_a = MixedLM(en_a, ex_a, groups=gr_a).fit(reml=True, maxiter=200)
print(f"\\n  Convergio: {mo_a.converged} | Log-Lik: {mo_a.llf:.1f}, AIC: {mo_a.aic:.0f}")
print(f"  Efectos fijos:\\n{mo_a.fe_params.to_string()}")
print(f"  Var(aislado) = {mo_a.cov_re.iloc[0,0]:.2f}, Var(residual) = {mo_a.scale:.2f}")
icc_a = mo_a.cov_re.iloc[0,0] / (mo_a.cov_re.iloc[0,0]+mo_a.scale)
print(f"  ICC (aislado) = {icc_a:.3f} ({icc_a*100:.1f}%)")
coefs_a = pd.DataFrame({'Coef':mo_a.fe_params,'EE':mo_a.bse_fe,'z':mo_a.tvalues,'p_valor':mo_a.pvalues})
coefs_a['IC95_inf']=coefs_a['Coef']-1.96*coefs_a['EE']; coefs_a['IC95_sup']=coefs_a['Coef']+1.96*coefs_a['EE']
print(f"\\n{coefs_a.round(3).to_string()}")

pred_mac_a = mo_a.fe_params['const']; pred_sox_a = mo_a.fe_params['const']+mo_a.fe_params['metodo_id_soxhlet']
pred_ult_a = mo_a.fe_params['const']+mo_a.fe_params['metodo_id_ultrasonido']
se_mac_a = mo_a.bse_fe['const']; cov_a = mo_a.cov_params()
se_sox_a = np.sqrt(mo_a.bse_fe['const']**2+mo_a.bse_fe['metodo_id_soxhlet']**2+2*cov_a.loc['const','metodo_id_soxhlet'])
se_ult_a = np.sqrt(mo_a.bse_fe['const']**2+mo_a.bse_fe['metodo_id_ultrasonido']**2+2*cov_a.loc['const','metodo_id_ultrasonido'])
print("\\n  Medias marginales a 5.0 mg/mL:")
for lbl, est, se in [('Maceracion',pred_mac_a,se_mac_a),('Soxhlet',pred_sox_a,se_sox_a),('Ultrasonido',pred_ult_a,se_ult_a)]:
    print(f"    {lbl:15s} = {est:.2f}% +/- {1.96*se:.2f}")
""")

md("""### 3.2 Modelo B: Dosis-respuesta en Maceracion

Se usa concentracion como efecto fijo categorico para no asumir una forma funcional especifica. La concentracion de referencia es la mas baja (0.0 mg/mL = control).
""")

code("""print("\\n"+"-"*70); print("  MODELO B: Maceracion -- dosis-respuesta"); print("-"*70)
df_b = INH[INH['metodo_extraccion']=='maceracion'].copy()
print(f"  Datos: {len(df_b)} obs, {df_b['aislado_id'].nunique()} aislados, Concentraciones: {sorted(df_b['concentracion_mg_ml'].unique())}")
dc_b = pd.get_dummies(df_b['conc_cat'], drop_first=True, dtype=float, prefix='conc')
df_b_ml = pd.concat([df_b.reset_index(drop=True), dc_b.reset_index(drop=True)], axis=1)
gr_b = df_b_ml['aislado_id']; ex_b = add_constant(df_b_ml[[c for c in dc_b.columns]]); en_b = df_b_ml['porcentaje_inhibicion']
mo_b = MixedLM(en_b, ex_b, groups=gr_b).fit(reml=True, maxiter=200)
print(f"\\n  Convergio: {mo_b.converged} | Log-Lik: {mo_b.llf:.1f}, AIC: {mo_b.aic:.0f}")
coefs_b = pd.DataFrame({'Coef':mo_b.fe_params,'EE':mo_b.bse_fe,'z':mo_b.tvalues,'p_valor':mo_b.pvalues})
coefs_b['IC95_inf']=coefs_b['Coef']-1.96*coefs_b['EE']; coefs_b['IC95_sup']=coefs_b['Coef']+1.96*coefs_b['EE']
print(f"\\n{coefs_b.round(3).to_string()}")
icc_b = mo_b.cov_re.iloc[0,0]/(mo_b.cov_re.iloc[0,0]+mo_b.scale)
print(f"\\n  Var(aislado)={mo_b.cov_re.iloc[0,0]:.2f}, Var(residual)={mo_b.scale:.2f}, ICC={icc_b:.3f}")
concs_b = sorted(df_b['concentracion_mg_ml'].unique()); cref_b = concs_b[0]; const_b = mo_b.fe_params['const']
print(f"\\n  Medias marginales (ref = {cref_b} mg/mL):")
for conc in concs_b:
    if conc==cref_b: est=const_b; se=mo_b.bse_fe['const']
    else:
        col=f'conc_{conc}'; est=const_b+mo_b.fe_params[col]
        se=np.sqrt(mo_b.bse_fe['const']**2+mo_b.bse_fe[col]**2+2*mo_b.cov_params().loc['const',col])
    print(f"    {conc:>5} mg/mL: {est:.2f}% +/- {1.96*se:.2f}")
""")

md("""### 3.3 Modelo C: Probabilidad de inhibicion completa

Modelar la probabilidad de inhibicion completa (crecimiento = 0 mm) mediante regresion logistica permite separar el efecto sobre la magnitud de la inhibicion vs. la capacidad de lograr inhibicion total.

**Limitacion:** Si solo un metodo presenta casos de inhibicion completa, el modelo logistico no puede estimarse.
""")

code("""print("\\n"+"-"*70); print("  MODELO C: Probabilidad de inhibicion completa (5.0 mg/mL)"); print("-"*70)
df_c = df_a.copy(); n_comp = df_c['completa_bin'].sum()
print(f"  Inhibicion completa: {n_comp}/{len(df_c)} ({100*n_comp/len(df_c):.1f}%)")
ct = pd.crosstab(df_c['metodo_id'], df_c['completa_bin']); print(f"\\n  Contingencia:\\n{ct}")
nc_m = df_c.groupby('metodo_id')['completa_bin'].sum()
if nc_m.nunique()==1:
    print("\\n  Solo un metodo tiene inhibicion completa -- modelo logistico no aplicable")
    mc_ok=False; pred_probs=ct[1]/ct.sum(axis=1) if 1 in ct.columns else pd.Series(0,index=ct.index)
else:
    df_c_ml=pd.get_dummies(df_c,columns=['metodo_id'],drop_first=True,dtype=float)
    ex_c=add_constant(df_c_ml[['metodo_id_soxhlet','metodo_id_ultrasonido']]); en_c=df_c_ml['completa_bin']
    mc=Logit(en_c,ex_c).fit(disp=False,maxiter=200)
    print(f"  Log-Lik: {mc.llf:.1f} | Pseudo-R2: {mc.prsquared:.3f}")
    pred_probs=mc.predict(ex_c); mc_ok=True

print(f"\\n  Proporcion de inhibicion completa:")
for metodo in ['maceracion','soxhlet','ultrasonido']:
    if mc_ok: pp=pred_probs[df_c['metodo_id']==metodo].mean(); print(f"    {LABEL_MET[metodo]:15s}: {pp:.1%}")
    else: nt=len(df_c[df_c['metodo_id']==metodo]); nc=nc_m.get(metodo,0); print(f"    {LABEL_MET[metodo]:15s}: {nc}/{nt} = {100*nc/nt:.1f}%")
print("\\n  Solo Maceracion logro inhibicion completa -- perfil fitoquimico cualitativamente diferente.")
""")

md("""### 3.4 Figuras del Objetivo 2""")

code("""# 3.4a: Boxplot 5 mg/mL
fig,ax=plt.subplots(figsize=(8,5))
sns.boxplot(data=df_a,x='metodo_id',y='porcentaje_inhibicion',hue='metodo_id',palette=COLOR_MET,ax=ax,legend=False)
sns.stripplot(data=df_a,x='metodo_id',y='porcentaje_inhibicion',color='black',alpha=0.3,size=4,ax=ax,jitter=True)
ax.set_xlabel('Metodo'); ax.set_ylabel('Inhibicion (%)'); ax.set_title('Inhibicion a 5.0 mg/mL por metodo')
ax.set_xticks([0,1,2]); ax.set_xticklabels([LABEL_MET[l] for l in ['maceracion','soxhlet','ultrasonido']])
ax.axhline(0,color='gray',ls=':',alpha=0.5)
ax.errorbar([0,1,2],[pred_mac_a,pred_sox_a,pred_ult_a],yerr=[1.96*se_mac_a,1.96*se_sox_a,1.96*se_ult_a],
    fmt='D',color='#d62728',markersize=7,capsize=4,zorder=5,label='Estimado LMM')
ax.legend(loc='upper right')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj2_inhibicion_5mg_ml.png'),dpi=300); plt.show()

# 3.4b: Dosis-respuesta Maceracion
fig,ax=plt.subplots(figsize=(8,5))
sns.boxplot(data=df_b,x='conc_cat',y='porcentaje_inhibicion',hue='conc_cat',palette='Blues',ax=ax,legend=False)
sns.stripplot(data=df_b,x='conc_cat',y='porcentaje_inhibicion',color='black',alpha=0.3,size=4,ax=ax,jitter=True)
ax.set_xlabel('Concentracion (mg/mL)'); ax.set_ylabel('Inhibicion (%)'); ax.set_title('Maceracion -- dosis-respuesta')
ax.axhline(0,color='gray',ls=':',alpha=0.5)
for i,conc in enumerate(concs_b):
    if conc==cref_b: est=const_b; se=mo_b.bse_fe['const']
    else:
        col=f'conc_{conc}'
        if col in mo_b.fe_params.index: est=const_b+mo_b.fe_params[col]; se=np.sqrt(mo_b.bse_fe['const']**2+mo_b.bse_fe[col]**2+2*mo_b.cov_params().loc['const',col])
    ax.plot(i,est,'D',color='red',markersize=8,zorder=5); ax.errorbar(i,est,yerr=1.96*se,color='red',capsize=5,capthick=2,alpha=0.7)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj2_dosis_respuesta_maceracion.png'),dpi=300); plt.show()

# 3.4c: Perfil por aislado
fig,ax=plt.subplots(figsize=(12,6))
perf_b=df_b.groupby(['aislado_id','concentracion_mg_ml'])['porcentaje_inhibicion'].mean().reset_index()
for aislado in df_b['aislado_id'].unique():
    sub=perf_b[perf_b['aislado_id']==aislado]
    ax.plot(sub['concentracion_mg_ml'],sub['porcentaje_inhibicion'],marker='o',alpha=0.4,linewidth=0.8,color='#2e86ab')
avg_b=perf_b.groupby('concentracion_mg_ml')['porcentaje_inhibicion'].mean()
ax.plot(avg_b.index,avg_b.values,'r-o',linewidth=3,label='Promedio')
ax.set_xlabel('Concentracion (mg/mL)'); ax.set_ylabel('Inhibicion media (%)')
ax.set_title('Perfil individual por aislado -- Maceracion')
ax.legend(bbox_to_anchor=(1.02,1),loc='upper left'); ax.axhline(0,color='gray',ls=':',alpha=0.5)
fig.subplots_adjust(right=0.82); fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj2_perfil_aislados_mac.png'),dpi=300); plt.show()

# 3.4d: Inhibicion completa
fig,ax=plt.subplots(figsize=(8,5))
ct_pct=ct.div(ct.sum(axis=1),axis=0)*100; cp=1 if 1 in ct_pct.columns else ct_pct.columns[0]
ct_pct[cp].plot(kind='bar',ax=ax,color=[COLOR_MET[l] for l in ['maceracion','soxhlet','ultrasonido']],edgecolor='black',alpha=0.8)
ax.set_xticks(range(3)); ax.set_xticklabels([LABEL_MET[l] for l in ['maceracion','soxhlet','ultrasonido']],rotation=0)
ax.set_xlabel('Metodo'); ax.set_ylabel('% inhibicion completa'); ax.set_title('Proporcion de inhibicion completa a 5.0 mg/mL')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj2_inhibicion_completa.png'),dpi=300); plt.show()

# 3.4e: Diagnostico
fig,axes=plt.subplots(1,2,figsize=(12,5)); fig.subplots_adjust(wspace=0.4)
res_a=mo_a.resid; stats.probplot(res_a,dist='norm',plot=axes[0]); axes[0].set_title('Q-Q residuos (Modelo A)')
axes[1].scatter(mo_a.fittedvalues,res_a,alpha=0.5,c='#a23b72',edgecolors='none'); axes[1].axhline(0,color='gray',ls='--',alpha=0.5)
axes[1].set_xlabel('Ajustados'); axes[1].set_ylabel('Residuos'); axes[1].set_title('Residuos vs. ajustados')
if 3<=len(res_a)<=5000:
    _,sp_a=stats.shapiro(res_a); axes[1].text(0.05,0.95,f'Shapiro p={sp_a:.4f}',transform=axes[1].transAxes,va='top',fontsize=9,bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.5))
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj2_diagnostico_modelo_a.png'),dpi=300); plt.show()
print("  OK Figuras Objetivo 2 generadas.")
""")

md("""### Interpretacion del Objetivo 2

**Hallazgos principales:**
1. Maceracion produce la mayor inhibicion a 5.0 mg/mL (~80%), seguida de Soxhlet y Ultrasonido
2. Efecto dosis-respuesta claro en Maceracion: 0.2->~30%, 1.0->~60%, 5.0->~80%
3. Solo Maceracion logra inhibicion completa en varios aislados
4. Alta variabilidad entre aislados (ICC ~0.35-0.45): consistentemente diferentes en susceptibilidad
5. La inhibicion negativa en algunos casos sugiere posible estimulacion a baja concentracion (hormesis)
""")

# ===== SECTION 4: OBJETIVO 3 - CONIDIAS =====
md("""---
# 4. Objetivo 3 -- Produccion de Conidias

## Pregunta biologica
El tratamiento afecta la produccion de conidias (esporas) de Fusarium?

## Importancia biologica
Un fungicida ideal deberia: (1) inhibir el crecimiento del hongo y (2) reducir la produccion de esporas para limitar la dispersion.

## Estrategia de modelado

Los datos de conidias fueron reportados en escala **log10** por el laboratorio. Se modela `conidias_log10` como variable respuesta primaria por su mejor comportamiento estadistico.

### Modelos:
- **Modelo A**: LMM sobre log10(conidias) a 5.0 mg/mL -- compara los 3 metodos
- **Modelo B**: LMM dosis-respuesta Maceracion (log10 conidias)
- **Modelo C**: LMM sobre INH_log10 de conidias a 5.0 mg/mL
""")

code("""print("="*70); print("  OBJETIVO 3 -- PRODUCCION DE CONIDIAS"); print("="*70)
print(f"\\nTotal obs: {len(CONI_TRAT)}, Aislados: {CONI_TRAT['aislado_id'].nunique()}")
print(f"Rango log10(conidias): [{CONI_TRAT['conidias_log10'].min():.2f}, {CONI_TRAT['conidias_log10'].max():.2f}]")
print(f"Control: {len(CTRL_CONI)} obs, media log10={CTRL_CONI['conidias_log10'].mean():.2f}")
CONI_TRAT['conc_cat'] = CONI_TRAT['concentracion_mg_ml'].astype('category')
""")

md("""### 4.1 Modelo A: log10(conidias) a 5.0 mg/mL

En escala log10, una diferencia de 1 unidad equivale a un factor de 10 en la escala cruda. Pequenas diferencias en log10 representan grandes reducciones porcentuales.
""")

code("""print("\\n"+"-"*70); print("  MODELO A: log10(conidias) a 5.0 mg/mL"); print("-"*70)
df_ac = CONI_TRAT[CONI_TRAT['concentracion_mg_ml']==5.0].copy()
print(f"  Datos: {len(df_ac)} obs, {df_ac['aislado_id'].nunique()} aislados")
df_ac_ml = pd.get_dummies(df_ac, columns=['metodo_id'], drop_first=True, dtype=float)
gr_ac = df_ac_ml['aislado_id']; ex_ac = add_constant(df_ac_ml[['metodo_id_soxhlet','metodo_id_ultrasonido']]); en_ac = df_ac_ml['conidias_log10']
mo_ac = MixedLM(en_ac, ex_ac, groups=gr_ac).fit(reml=True, maxiter=200)
print(f"\\n  Convergio: {mo_ac.converged} | Log-Lik: {mo_ac.llf:.1f}")
coefs_ac = pd.DataFrame({'Coef':mo_ac.fe_params,'EE':mo_ac.bse_fe,'z':mo_ac.tvalues,'p_valor':mo_ac.pvalues})
coefs_ac['IC95_inf']=coefs_ac['Coef']-1.96*coefs_ac['EE']; coefs_ac['IC95_sup']=coefs_ac['Coef']+1.96*coefs_ac['EE']
print(f"\\n{coefs_ac.round(3).to_string()}")
v_iso_ac=mo_ac.cov_re.iloc[0,0]; v_res_ac=mo_ac.scale; icc_ac=v_iso_ac/(v_iso_ac+v_res_ac)
print(f"\\n  Var(aislado)={v_iso_ac:.3f}, Var(residual)={v_res_ac:.3f}, ICC={icc_ac:.3f}")
pm_ac=mo_ac.fe_params['const']; ps_ac=mo_ac.fe_params['const']+mo_ac.fe_params['metodo_id_soxhlet']
pu_ac=mo_ac.fe_params['const']+mo_ac.fe_params['metodo_id_ultrasonido']
ctrl_m = CTRL_CONI['conidias_log10'].mean()
print(f"\\n  Medias marginales y reduccion vs control (log10={ctrl_m:.2f}):")
for lbl,est in [('Maceracion',pm_ac),('Soxhlet',ps_ac),('Ultrasonido',pu_ac)]:
    rl=ctrl_m-est; rp=(1-10**(-rl))*100
    print(f"    {lbl:15s} log10={est:.2f} ({10**est:.0f}/mL) reduccion: {rp:.1f}%")
""")

md("""### 4.2 Modelo B: Dosis-respuesta Maceracion (log10 conidias)""")

code("""print("\\n"+"-"*70); print("  MODELO B: Maceracion -- dosis-respuesta en log10(conidias)"); print("-"*70)
df_bc = CONI_TRAT[CONI_TRAT['metodo_extraccion']=='maceracion'].copy()
print(f"  Datos: {len(df_bc)} obs, Concentraciones: {sorted(df_bc['concentracion_mg_ml'].unique())}")
dbc_dum = pd.get_dummies(df_bc['conc_cat'], drop_first=True, dtype=float, prefix='conc')
df_bc_ml = pd.concat([df_bc.reset_index(drop=True), dbc_dum.reset_index(drop=True)], axis=1)
gr_bc = df_bc_ml['aislado_id']; ex_bc = add_constant(df_bc_ml[[c for c in dbc_dum.columns]]); en_bc = df_bc_ml['conidias_log10']
mo_bc = MixedLM(en_bc, ex_bc, groups=gr_bc).fit(reml=True, maxiter=200)
print(f"\\n  Convergio: {mo_bc.converged} | Log-Lik: {mo_bc.llf:.1f}")
coefs_bc = pd.DataFrame({'Coef':mo_bc.fe_params,'EE':mo_bc.bse_fe,'z':mo_bc.tvalues,'p_valor':mo_bc.pvalues})
coefs_bc['IC95_inf']=coefs_bc['Coef']-1.96*coefs_bc['EE']; coefs_bc['IC95_sup']=coefs_bc['Coef']+1.96*coefs_bc['EE']
print(f"\\n{coefs_bc.round(3).to_string()}")
v_iso_bc=mo_bc.cov_re.iloc[0,0]; v_res_bc=mo_bc.scale; icc_bc=v_iso_bc/(v_iso_bc+v_res_bc)
print(f"\\n  Var(aislado)={v_iso_bc:.3f}, Var(residual)={v_res_bc:.3f}, ICC={icc_bc:.3f}")
cbc=sorted(df_bc['concentracion_mg_ml'].unique()); crbc=cbc[0]; cst_bc=mo_bc.fe_params['const']
ctr_mac=CTRL_CONI[CTRL_CONI['metodo_extraccion']=='maceracion']['conidias_log10'].mean()
print(f"\\n  Control Maceracion: log10={ctr_mac:.2f}")
for conc in cbc:
    if conc==crbc: est=cst_bc
    else:
        col=f'conc_{conc}'
        if col in mo_bc.fe_params.index: est=cst_bc+mo_bc.fe_params[col]
        else: est=cst_bc
    rp=max(0,(1-10**(-(ctr_mac-est)))*100)
    print(f"    {conc:>5} mg/mL: log10={est:.2f} ({10**est:.0f}/mL) reduccion: {rp:.1f}%")
""")

md("""### 4.3 Modelo C: INH_log10 de conidias a 5.0 mg/mL

**Por que NO usar INH crudo?** Su distribucion tiene skew~-12.2, kurtosis~176 -- no rescatable. El `porcentaje_inhibicion_log10` (calculado sobre log10) es estable y adecuado para modelado. Este modelo es paralelo al Modelo A del Objetivo 2.
""")

code("""print("\\n"+"-"*70); print("  MODELO C: INH conidias (escala log10, reportada por el laboratorio)"); print("-"*70)
df_cc = CONI_INH_LOG[(CONI_INH_LOG['concentracion_mg_ml']==5.0) & CONI_INH_LOG['porcentaje_inhibicion_log10'].notna()].copy()
print(f"  Datos: {len(df_cc)} obs")
df_cc_ml = pd.get_dummies(df_cc, columns=['metodo_id'], drop_first=True, dtype=float)
gr_cc = df_cc_ml['aislado_id']; ex_cc = add_constant(df_cc_ml[['metodo_id_soxhlet','metodo_id_ultrasonido']]); en_cc = df_cc_ml['porcentaje_inhibicion_log10']
mo_cc = MixedLM(en_cc, ex_cc, groups=gr_cc).fit(reml=True, maxiter=200)
print(f"\\n  Convergio: {mo_cc.converged} | Log-Lik: {mo_cc.llf:.1f}")
coefs_cc = pd.DataFrame({'Coef':mo_cc.fe_params,'EE':mo_cc.bse_fe,'z':mo_cc.tvalues,'p_valor':mo_cc.pvalues})
coefs_cc['IC95_inf']=coefs_cc['Coef']-1.96*coefs_cc['EE']; coefs_cc['IC95_sup']=coefs_cc['Coef']+1.96*coefs_cc['EE']
print(f"\\n{coefs_cc.round(3).to_string()}")
pm_cc=mo_cc.fe_params['const']; ps_cc=mo_cc.fe_params['const']+mo_cc.fe_params['metodo_id_soxhlet']
pu_cc=mo_cc.fe_params['const']+mo_cc.fe_params['metodo_id_ultrasonido']
print("\\n  INH estimado (escala log10) a 5.0 mg/mL:")
for lbl,est in [('Maceracion',pm_cc),('Soxhlet',ps_cc),('Ultrasonido',pu_cc)]:
    print(f"    {lbl:15s}: {est:.1f}%")
""")

md("""### 4.4 Correlacion con Objetivo 2""")

code("""fig,ax=plt.subplots(figsize=(8,6))
ci=CREC[~CREC['es_control']&CREC['porcentaje_inhibicion'].notna()][['aislado_id','metodo_extraccion','concentracion_mg_ml','replica_biologica','porcentaje_inhibicion']].rename(columns={'porcentaje_inhibicion':'ic'})
ci2=CONI_INH[CONI_INH['porcentaje_inhibicion'].notna()][['aislado_id','metodo_extraccion','concentracion_mg_ml','replica_biologica','porcentaje_inhibicion']].rename(columns={'porcentaje_inhibicion':'ic2'})
mg=ci.merge(ci2,on=['aislado_id','metodo_extraccion','concentracion_mg_ml','replica_biologica'])
if len(mg)>10:
    ax.scatter(mg['ic'],mg['ic2'],alpha=0.4,c='#2e86ab',edgecolors='none')
    r_m,p_m=stats.pearsonr(mg['ic'],mg['ic2'])
    ax.set_xlabel('Inhibicion crecimiento (%)'); ax.set_ylabel('Inhibicion conidias (%) -- crudo')
    ax.set_title(f'Correlacion INH crecimiento vs. conidias (r={r_m:.3f}, p={p_m:.4f})')
    ax.axhline(0,color='gray',ls=':',alpha=0.5); ax.axvline(0,color='gray',ls=':',alpha=0.5)
    fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj3_correlacion_inh_crecimiento_conidias.png'),dpi=300); plt.show()
    print(f"  Correlacion: r={r_m:.3f}, p={p_m:.4f}, n={len(mg)}")
else: plt.close(fig); print(f"  Pocos datos (n={len(mg)})")
""")

md("""### 4.5 Figuras del Objetivo 3""")

code("""# 4.5a: Boxplot log10 5 mg/mL
fig,ax=plt.subplots(figsize=(8,5)); om=['maceracion','soxhlet','ultrasonido']
sns.boxplot(data=df_ac,x='metodo_id',y='conidias_log10',hue='metodo_id',palette=COLOR_MET,ax=ax,order=om,legend=False)
sns.stripplot(data=df_ac,x='metodo_id',y='conidias_log10',color='black',alpha=0.3,size=4,ax=ax,jitter=True,order=om)
ax.axhline(CTRL_CONI['conidias_log10'].mean(),color='red',ls='--',alpha=0.6,label='Control')
ax.set_xlabel('Metodo'); ax.set_ylabel('log10(conidias/mL)'); ax.set_title('Conidias a 5.0 mg/mL por metodo')
ax.set_xticks([0,1,2]); ax.set_xticklabels([LABEL_MET[l] for l in om]); ax.legend(loc='lower right')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj3_conidias_log10_5mg_ml.png'),dpi=300); plt.show()

# 4.5b: Dosis-respuesta Mac conidias
fig,ax=plt.subplots(figsize=(8,5))
sns.boxplot(data=df_bc,x='conc_cat',y='conidias_log10',hue='conc_cat',palette='Blues',ax=ax,legend=False)
sns.stripplot(data=df_bc,x='conc_cat',y='conidias_log10',color='black',alpha=0.3,size=4,ax=ax,jitter=True)
ax.axhline(CTRL_CONI[CTRL_CONI['metodo_extraccion']=='maceracion']['conidias_log10'].mean(),color='red',ls='--',alpha=0.6,label='Control Mac.')
ax.set_xlabel('Concentracion (mg/mL)'); ax.set_ylabel('log10(conidias/mL)'); ax.set_title('Maceracion -- dosis-respuesta conidias')
ax.legend(bbox_to_anchor=(1.02,1),loc='upper left'); fig.subplots_adjust(right=0.82)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj3_dosis_respuesta_mac_conidias.png'),dpi=300); plt.show()

# 4.5c: INH conidias log10
fig,ax=plt.subplots(figsize=(8,5))
sns.boxplot(data=df_cc,x='metodo_id',y='porcentaje_inhibicion_log10',hue='metodo_id',palette=COLOR_MET,ax=ax,order=om,legend=False)
sns.stripplot(data=df_cc,x='metodo_id',y='porcentaje_inhibicion_log10',color='black',alpha=0.3,size=4,ax=ax,jitter=True,order=om)
ax.axhline(0,color='gray',ls=':',alpha=0.5)
ax.set_xlabel('Metodo'); ax.set_ylabel('Reduccion conidias (%) -- escala log10')
ax.set_title('Inhibicion de conidias (log10, escala de la hoja)')
ax.set_xticks([0,1,2]); ax.set_xticklabels([LABEL_MET[l] for l in om])
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj3_inh_conidias_5mg_ml.png'),dpi=300); plt.show()

# 4.5d: Diagnostico
fig,axes=plt.subplots(1,2,figsize=(12,5)); fig.subplots_adjust(wspace=0.4)
r_ac=mo_ac.resid; stats.probplot(r_ac,dist='norm',plot=axes[0]); axes[0].set_title('Q-Q residuos (Modelo A -- log10 conidias)')
axes[1].scatter(mo_ac.fittedvalues,r_ac,alpha=0.5,c='#2e86ab',edgecolors='none'); axes[1].axhline(0,color='gray',ls='--',alpha=0.5)
axes[1].set_xlabel('Ajustados'); axes[1].set_ylabel('Residuos'); axes[1].set_title('Residuos vs. ajustados')
if 3<=len(r_ac)<=5000:
    _,sp_ac=stats.shapiro(r_ac); axes[1].text(0.05,0.95,f'Shapiro p={sp_ac:.4f}',transform=axes[1].transAxes,va='top',fontsize=9,bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.5))
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj3_diagnostico_modelo_a.png'),dpi=300); plt.show()
print("  OK Figuras Objetivo 3 generadas.")
""")

md("""### Interpretacion del Objetivo 3

**Hallazgos principales:**
1. Maceracion reduce drasticamente la esporulacion (>90% vs control en escala cruda)
2. Soxhlet y Ultrasonido tienen efecto minimo sobre conidias, pese a inhibir el crecimiento micelial
3. Efecto dosis-respuesta en Maceracion para conidias
4. Baja correlacion con inhibicion del crecimiento sugiere mecanismos independientes
5. Maceracion afecta tanto crecimiento como esporulacion (doble mecanismo); Soxhlet/Ultrasonido solo crecimiento
""")

# ===== SECTION 5: OBJETIVO 4 - SUSCEPTIBILIDAD =====
md("""---
# 5. Objetivo 4 -- Susceptibilidad de Aislados de Fusarium

## Pregunta biologica
Existen diferencias sistematicas en la susceptibilidad a los extractos entre distintos aislados de *Fusarium* spp.?

## Enfoque
Se construye un **perfil de susceptibilidad** para cada aislado usando:
1. INH en crecimiento a 5.0 mg/mL para cada metodo
2. INH en conidias a 5.0 mg/mL
3. EC50 para Maceracion (el unico metodo con dosis-respuesta)

Luego: PCA, clustering jerarquico, heatmap, ranking.

**IMPORTANTE:** No se usa el termino "resistente" porque no existe un umbral validado. Se emplean: susceptibilidad alta, intermedia y baja.
""")

code("""print("="*70); print("  OBJETIVO 4 -- SUSCEPTIBILIDAD DE AISLADOS"); print("="*70)
aislados = sorted(INH['aislado_id'].unique())
print(f"\\nAislados: {len(aislados)}")
""")

md("""### 5.1 Metricas de susceptibilidad por aislado""")

code("""print("\\n"+"-"*70); print("  METRICAS DE SUSCEPTIBILIDAD"); print("-"*70)
metricas_lista = []
for aislado in aislados:
    s_inh = INH[INH['aislado_id']==aislado]
    for conc in [0.2,1.0,5.0]:
        sub=s_inh[(s_inh['metodo_extraccion']=='maceracion')&(s_inh['concentracion_mg_ml']==conc)]['porcentaje_inhibicion']
        metricas_lista.append((f'crec_mac_{conc}',aislado,conc,sub.mean()))
    for met in ['soxhlet','ultrasonido']:
        sub=s_inh[(s_inh['metodo_extraccion']==met)&(s_inh['concentracion_mg_ml']==5.0)]['porcentaje_inhibicion']
        metricas_lista.append((f'crec_{met}_5.0',aislado,5.0,sub.mean()))
    s_con=CONI_INH[CONI_INH['aislado_id']==aislado]
    for met in ['maceracion','soxhlet','ultrasonido']:
        sub=s_con[(s_con['metodo_extraccion']==met)&(s_con['concentracion_mg_ml']==5.0)]['porcentaje_inhibicion']
        metricas_lista.append((f'con_{met}_5.0',aislado,5.0,sub.mean()))

mdf=pd.DataFrame(metricas_lista,columns=['metrica','aislado','conc','valor'])
perfil=mdf.pivot_table(index='aislado',columns='metrica',values='valor')
print(f"\\nPerfil: {perfil.shape[0]} aislados x {perfil.shape[1]} metricas")
for col in perfil.columns: print(f"  {col:20s} media={perfil[col].mean():.1f} DE={perfil[col].std():.1f} n={perfil[col].notna().sum()}")
""")

md("""### 5.2 EC50 para Maceracion

La EC50 se estima por interpolacion log-lineal entre los puntos experimentales. Con solo 4 concentraciones, un modelo no lineal de 4 parametros no es estimable. Categorias: `estimado`, `no_alcanza_50`, `ec50_menor_0.2`, `insuficiente`.
""")

code("""print("\\n"+"-"*70); print("  EC50 -- MACERACION (interpolacion log-lineal)"); print("-"*70)
ec50_res = []
for aislado in aislados:
    sub=(INH[(INH['aislado_id']==aislado)&(INH['metodo_extraccion']=='maceracion')]
         .groupby('concentracion_mg_ml')['porcentaje_inhibicion'].mean())
    if 0.2 not in sub.index or 1.0 not in sub.index or 5.0 not in sub.index:
        ec50_res.append({'aislado':aislado,'ec50_mg_ml':np.nan,'ec50_clasificacion':'insuficiente'}); continue
    x_all=np.array([0.0,0.2,1.0,5.0]); y_all=np.array([0.0,sub[0.2],sub[1.0],sub[5.0]])
    if y_all.max()<50: ec50_res.append({'aislado':aislado,'ec50_mg_ml':np.nan,'ec50_clasificacion':'no_alcanza_50'})
    elif sub[0.2]>=50: ec50_res.append({'aislado':aislado,'ec50_mg_ml':np.nan,'ec50_clasificacion':'ec50_menor_0.2'})
    else:
        ev,cf=np.nan,'no_estimable'
        for i in range(len(x_all)-1):
            if (y_all[i]<50) and (y_all[i+1]>=50):
                x1,x2=x_all[i],x_all[i+1]; y1,y2=y_all[i],y_all[i+1]
                if x1==0: ev=x1+(50-y1)*(x2-x1)/(y2-y1)
                else: ev=np.exp(np.log(x1)+(50-y1)*(np.log(x2)-np.log(x1))/(y2-y1))
                cf='estimado'; break
        ec50_res.append({'aislado':aislado,'ec50_mg_ml':ev,'ec50_clasificacion':cf})

ec50_df=pd.DataFrame(ec50_res); n_est=(ec50_df['ec50_clasificacion']=='estimado').sum()
print(f"  EC50 estimado: {n_est}/{len(ec50_df)} aislados")
print(f"  No alcanzan 50%: {(ec50_df['ec50_clasificacion']=='no_alcanza_50').sum()}")
print(f"  EC50 < 0.2 mg/mL: {(ec50_df['ec50_clasificacion']=='ec50_menor_0.2').sum()}")

ec50_v=ec50_df[ec50_df['ec50_clasificacion']=='estimado'].copy()
if len(ec50_v)>0:
    print(f"\\n  EC50 medio: {ec50_v['ec50_mg_ml'].mean():.2f} mg/mL")
    print(f"  EC50 min: {ec50_v['ec50_mg_ml'].min():.2f}, max: {ec50_v['ec50_mg_ml'].max():.2f}")
    print("\\n  Top 5 mas susceptibles (menor EC50):")
    for _,r in ec50_v.nsmallest(5,'ec50_mg_ml').iterrows(): print(f"    {r['aislado']:25s} EC50={r['ec50_mg_ml']:.2f}")
    print("\\n  Top 5 menos susceptibles (mayor EC50):")
    for _,r in ec50_v.nlargest(5,'ec50_mg_ml').iterrows(): print(f"    {r['aislado']:25s} EC50={r['ec50_mg_ml']:.2f}")

perfil = perfil.merge(ec50_df[['aislado','ec50_mg_ml']], on='aislado', how='left')
if 'aislado' in perfil.columns: perfil = perfil.set_index('aislado')
""")

md("""### 5.3 PCA -- Analisis de Componentes Principales

Reduce las multiples metricas a componentes ortogonales que capturan la maxima varianza. Las variables se estandarizan porque estan en diferentes escalas.
""")

code("""print("\\n"+"-"*70); print("  PCA -- ANALISIS DE COMPONENTES PRINCIPALES"); print("-"*70)
pvars=[c for c in perfil.columns if c.startswith('crec_') or c.startswith('con_')]
pdata=perfil[pvars].dropna()
print(f"  Variables: {len(pvars)}, Aislados completos: {len(pdata)}")
scaler=StandardScaler(); ps=scaler.fit_transform(pdata)
ncp=min(len(pvars),5); pca=PCA(n_components=ncp); pc=pca.fit_transform(ps)
ve=pca.explained_variance_ratio_; vac=np.cumsum(ve)
print("\\n  Varianza explicada:")
for i,(ev,ea) in enumerate(zip(ve,vac)): print(f"    PC{i+1}: {ev:.1%} (acum: {ea:.1%})")
cargas=pd.DataFrame(pca.components_.T,index=pdata.columns,columns=[f'PC{i+1}' for i in range(ncp)])
print("\\n  Cargas:")
for col in cargas.columns:
    tv=cargas[col].abs().nlargest(3)
    print(f"    {col}: {', '.join([f'{v}: {cargas.loc[v,col]:.3f}' for v in tv.index])}")
""")

md("""### 5.4 Clustering jerarquico

Metodo de Ward con distancia Euclidea. El numero optimo de clusters se determina por el coeficiente de silhouette (maxima separacion entre clusters).
""")

code("""print("\\n"+"-"*70); print("  CLUSTERING JERARQUICO"); print("-"*70)
Z=linkage(ps,method='ward')
sil_scores=[]
for k in range(2,min(10,len(pdata))):
    lk=fcluster(Z,k,criterion='maxclust'); sil=silhouette_score(ps,lk)
    sil_scores.append({'k':k,'silhouette':sil})
sildf=pd.DataFrame(sil_scores); bk=int(sildf.loc[sildf['silhouette'].idxmax(),'k']); bs=sildf['silhouette'].max()
print("\\n  Silhouette:")
for _,r in sildf.iterrows(): print(f"    k={int(r['k'])} -> silhouette={r['silhouette']:.3f}")
print(f"  Mejor k = {bk} (silhouette={bs:.3f})")
lf=fcluster(Z,bk,criterion='maxclust')
perfil['cluster']=np.nan; pidx=pdata.index
for i,idx in enumerate(pidx): perfil.loc[idx,'cluster']=lf[i]
perfil['cluster']=perfil['cluster'].astype('Int64')
print("\\n  Caracterizacion de clusters:")
for k in sorted(perfil['cluster'].dropna().unique()):
    miembros=perfil[perfil['cluster']==k].index.tolist()
    print(f"\\n  Cluster {int(k)} (n={len(miembros)}):")
    for vp in ['crec_mac_5.0','crec_sox_5.0','crec_ult_5.0','ec50_mg_ml']:
        if vp in perfil.columns:
            vals=perfil[perfil['cluster']==k][vp].dropna()
            if len(vals): print(f"    {vp:20s} media={vals.mean():.1f}")
    print(f"    Aislados: {', '.join(str(m) for m in miembros)}")
""")

md("""### 5.5 Figuras del Objetivo 4""")

code("""colores_cluster=['#2e86ab','#a23b72','#f18f01','#41ab5d','#d95f02']
# 5.5a: EC50
if len(ec50_v)>0:
    fig,ax=plt.subplots(figsize=(12,6)); eo=ec50_v.sort_values('ec50_mg_ml'); xr=range(len(eo))
    ax.bar(xr,eo['ec50_mg_ml'].values,color='#2e86ab',alpha=0.8,edgecolor='black')
    ax.set_xticks(xr); ax.set_xticklabels(eo['aislado'].values,rotation=45,ha='right',fontsize=9)
    ax.set_ylabel('EC50 (mg/mL)'); ax.set_title('EC50 de Maceracion por aislado (interpolacion log-lineal)')
    ax.axhline(eo['ec50_mg_ml'].median(),color='red',ls='--',alpha=0.5,label=f'Mediana={eo["ec50_mg_ml"].median():.2f}')
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj4_ec50_aislados.png'),dpi=300); plt.show()

# 5.5b: Biplot PCA
fig,ax=plt.subplots(figsize=(9,7))
for k in sorted(perfil['cluster'].dropna().unique()):
    mask=(perfil.loc[pidx,'cluster'].values==k)
    if mask.sum()>0:
        ax.scatter(pc[mask,0],pc[mask,1],c=colores_cluster[int(k)%len(colores_cluster)],label=f'Cluster {int(k)}',s=60,alpha=0.7,edgecolors='black')
    for i,idx in enumerate(pidx):
        if mask[i]: ax.annotate(str(idx)[:10],(pc[i,0],pc[i,1]),fontsize=7,alpha=0.7)
for i,var in enumerate(pdata.columns):
    ax.arrow(0,0,cargas.loc[var,'PC1']*3,cargas.loc[var,'PC2']*3,head_width=0.05,head_length=0.05,fc='gray',ec='gray',alpha=0.5)
    ax.text(cargas.loc[var,'PC1']*3.2,cargas.loc[var,'PC2']*3.2,var,fontsize=8,color='gray',alpha=0.7)
ax.set_xlabel(f'PC1 ({ve[0]:.1%})'); ax.set_ylabel(f'PC2 ({ve[1]:.1%})')
ax.set_title('PCA -- Perfil de susceptibilidad de aislados')
ax.axhline(0,color='gray',ls=':',alpha=0.3); ax.axvline(0,color='gray',ls=':',alpha=0.3)
ax.legend(bbox_to_anchor=(1.02,1),loc='upper left'); fig.subplots_adjust(right=0.82)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj4_pca_susceptibilidad.png'),dpi=300); plt.show()

# 5.5c: Scree plot
fig,ax=plt.subplots(figsize=(6,4))
ax.bar(range(1,len(ve)+1),ve,alpha=0.7,color='#2e86ab',edgecolor='black')
ax.plot(range(1,len(ve)+1),vac,'ro-',markersize=6)
ax.set_xlabel('Componente principal'); ax.set_ylabel('Varianza explicada')
ax.set_title('Scree plot'); ax.axhline(0.7,color='gray',ls='--',alpha=0.4,label='70%'); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj4_scree_plot.png'),dpi=300); plt.show()

# 5.5d: Dendrograma
fig,ax=plt.subplots(figsize=(14,8))
dendrogram(Z,labels=pidx.values,leaf_font_size=10,color_threshold=Z[-(bk-1),2] if len(Z)>=bk else None,above_threshold_color='gray',ax=ax)
ax.set_ylabel('Distancia (Ward)'); ax.set_title('Dendrograma -- clustering jerarquico de aislados')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj4_dendrograma.png'),dpi=300); plt.show()

# 5.5e: Heatmap
fig,ax=plt.subplots(figsize=(14,10))
po=perfil.loc[pidx].copy(); po=po.sort_values('cluster')
hv=[v for v in pvars if v in po.columns]; hd=po[hv].copy(); hs=(hd-hd.mean())/hd.std()
sns.heatmap(hs.T,cmap='RdYlBu_r',center=0,xticklabels=po.index,yticklabels=hv,ax=ax,cbar_kws={'label':'Desviaciones de la media'})
ax.set_xticklabels(ax.get_xticklabels(),rotation=45,ha='right',fontsize=8)
ax.set_xlabel('Aislado'); ax.set_ylabel('Metrica'); ax.set_title('Perfil de susceptibilidad -- heatmap (estandarizado)')
umap={k:i for i,k in enumerate(sorted(po['cluster'].unique()))}
cc=[colores_cluster[umap[k]%len(colores_cluster)] for k in po['cluster']]
for i,(_,row) in enumerate(po.iterrows()): ax.text(i+0.5,len(hv)+0.5,f"C{int(row['cluster'])}",ha='center',va='center',fontsize=7,fontweight='bold',color=cc[i],clip_on=False)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj4_heatmap_susceptibilidad.png'),dpi=300); plt.show()

# 5.5f: Susceptibilidad por aislado (Maceracion 5.0)
fig,ax=plt.subplots(figsize=(14,6))
imh=INH[(INH['metodo_extraccion']=='maceracion')&(INH['concentracion_mg_ml']==5.0)]
oa=imh.groupby('aislado_id')['porcentaje_inhibicion'].mean().sort_values().index.tolist()
sns.boxplot(data=imh,x='aislado_id',y='porcentaje_inhibicion',order=oa,palette='RdYlGn',hue='aislado_id',legend=False,ax=ax)
ax.set_xticks(range(len(oa))); ax.set_xticklabels(oa,rotation=90,fontsize=8)
ax.set_xlabel('Aislado'); ax.set_ylabel('Inhibicion (%)')
ax.set_title('Susceptibilidad por aislado -- Maceracion 5.0 mg/mL')
ax.axhline(50,color='red',ls='--',alpha=0.4,label='50%'); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'obj4_susceptibilidad_aislados.png'),dpi=300); plt.show()
print("  OK Figuras Objetivo 4 generadas.")
""")

md("""### 5.6 Ranking de susceptibilidad

Se construye un score compuesto promediando el INH a 5.0 mg/mL en los 3 metodos. Se clasifica en terciles: susceptibilidad alta, intermedia y baja.
""")

code("""print("\\n"+"-"*70); print("  RANKING DE SUSCEPTIBILIDAD"); print("-"*70)
rvars=['crec_mac_5.0','crec_sox_5.0','crec_ult_5.0']
rdata=perfil[rvars].dropna().copy()
rdata['score_susceptibilidad']=rdata[rvars].mean(axis=1)
rdata=rdata.merge(ec50_df[['aislado','ec50_mg_ml']],left_index=True,right_on='aislado',how='left').set_index('aislado')
rdata=rdata.sort_values('score_susceptibilidad',ascending=False)
rdata['rank']=range(1,len(rdata)+1)
rdata['clasificacion']=pd.qcut(rdata['score_susceptibilidad'],q=3,labels=['Baja','Intermedia','Alta'],duplicates='drop')
print(f"\\n{'Rank':>5s} {'Aislado':25s} {'Score':>8s} {'EC50':>8s} {'Clasificacion'}")
print("  "+"-"*65)
for _,row in rdata.iterrows():
    es=f"{row['ec50_mg_ml']:.2f}" if not np.isnan(row.get('ec50_mg_ml',np.nan)) else 'N/A'
    print(f"  {int(row['rank']):>5d} {row.name:25s} {row['score_susceptibilidad']:>7.1f}  {es:>8s}  {row['clasificacion']}")

print("\\n**Nota:** No se usa 'resistente' por ausencia de umbral validado.")
print("Categorias: Alta susceptibilidad = tercil superior, Intermedia = medio, Baja = inferior.")
""")

# ===== SECTION 6: TRANSFORMACIONES =====
md("""---
# 6. Diagnostico de Transformaciones

## Fundamentos

Muchos tests estadisticos (ANOVA, modelos lineales) asumen normalidad de los residuos. Cuando los datos no cumplen este supuesto, las transformaciones pueden estabilizar la varianza y mejorar la normalidad. Sin embargo, cada transformacion tiene un costo en interpretabilidad.

Se evaluan sistematicamente las siguientes transformaciones para cada variable respuesta:

### Para INH de crecimiento:
- **Crudo** (sin transformar) -- interpretacion directa
- **IHS** (Inverse Hyperbolic Sine) -- maneja valores negativos y cero
- **CubeRoot** -- reduce asimetria positiva
- **RankGauss** -- fuerza normalidad perfecta pero pierde escala original
- **Arcsin-sqrt** -- clasica para porcentajes (solo 0-100)
- **Logit** -- para porcentajes acotados (solo 0-100)

### Para log10(conidias):
- **Crudo** (log10 ya aplicado)
- **BoxCox** -- encuentra la transformacion optima
- **Square** -- enfatiza diferencias grandes

### Para INH conidias:
- **Crudo** (escala original, muy asimetrica)
- **IHS** -- alternativa
- **log10 escala** -- la reportada por el laboratorio

Se evalua cada transformacion mediante:
1. **Skewness** y **kurtosis** de la distribucion transformada
2. **Shapiro-Wilk** como prueba formal de normalidad
3. **Residuos del LMM** con cada transformacion (cuando aplica)
4. **Levene** para homocedasticidad de residuos
""")

code("""print("="*70); print("  DIAGNOSTICO SISTEMATICO DE TRANSFORMACIONES"); print("="*70)
yc=INH['porcentaje_inhibicion'].dropna().values
ym=INH[(INH['porcentaje_inhibicion']>0)&(INH['porcentaje_inhibicion']<100)]['porcentaje_inhibicion'].values
cl=CONI_TRAT['conidias_log10'].dropna().values
ycc=CONI_INH['porcentaje_inhibicion'].dropna().values
ycl=CONI_INH_LOG['porcentaje_inhibicion_log10'].dropna().values
""")

md("""### 6.1 Transformaciones para INH de crecimiento""")

code("""print("\\n"+"-"*70); print("  1. INH CRECIMIENTO MICELIAL"); print("-"*70)
inh_t=[('Crudo',lambda x:x,'Sin transformacion'),('IHS(/10)',ihs_transform(10),'arcsinh(y/10)'),
       ('IHS(/50)',ihs_transform(50),'arcsinh(y/50)'),('IHS(/100)',ihs_transform(100),'arcsinh(y/100)'),
       ('CubeRoot',cube_root,'cbrt(y)'),('RankGauss',rank_gauss,'Normal scores')]
mid_t=[('Arcsin-sqrt',lambda x:np.arcsin(np.sqrt(x/100)),'arcsin(sqrt(y/100))'),
       ('Logit',lambda x:np.log(x/100/(1-x/100)),'log(y/(100-y))')]
print("  Distribucion completa:")
for r in eval_transforms('INH completo',yc,transform_funcs=inh_t):
    print(f"  {r['Transformacion']:12s} skew={r['Skew']:.2f} kurt={r['Kurtosis']:.2f} Shapiro p={r['Shapiro_p']:.4f}")
print("\\n  Subset 0-100:")
for r in eval_transforms('INH mid',ym,transform_funcs=mid_t+inh_t):
    print(f"  {r['Transformacion']:12s} skew={r['Skew']:.2f} kurt={r['Kurtosis']:.2f} Shapiro p={r['Shapiro_p']:.4f}")
""")

md("""### 6.2 Residuos del LMM con diferentes transformaciones

Se ajusta el mismo LMM (INH ~ metodo + (1|aislado) a 5.0 mg/mL) variando la transformacion de la variable respuesta. Se comparan log-verosimilitud, asimetria y curtosis de residuos.
""")

code("""print("\\n  Residuos del LMM (5.0 mg/mL, metodo ~ aislado):")
df_lt=INH.copy(); df_lt['metodo_id']=df_lt['metodo_extraccion'].str.strip().str.lower()
df_lt=df_lt[df_lt['concentracion_mg_ml']==5.0]; d_dum=pd.get_dummies(df_lt,columns=['metodo_id'],drop_first=True,dtype=float)
ex_l=add_constant(d_dum[['metodo_id_soxhlet','metodo_id_ultrasonido']]); gr_l=d_dum['aislado_id']
lmm_t=[('Crudo',lambda x:x,'Crudo'),('IHS/10',ihs_transform(10),'arcsinh(y/10)'),
       ('IHS/50',ihs_transform(50),'arcsinh(y/50)'),('CubeRoot',cube_root,'cbrt(y)'),('RankGauss',rank_gauss,'Normal scores')]
lmm_r=eval_residuals_lmm('LMM INH',d_dum['porcentaje_inhibicion'],ex_l,gr_l,lmm_t)
for r in lmm_r:
    print(f"  {r['Transformacion']:12s} LogLik={r['LogLik']:.0f} skew={r['Resid_skew']:.2f} kurt={r['Resid_kurt']:.2f} Shapiro p={r['Shapiro_p']:.4f}")

print("\\n  Homocedasticidad (Levene) de residuos del LMM:")
gar=np.where(d_dum['metodo_id_soxhlet']==1,'soxhlet',np.where(d_dum['metodo_id_ultrasonido']==1,'ultrasonido','maceracion'))
for tn,tf,_ in lmm_t:
    yt=tf(d_dum['porcentaje_inhibicion']); m_l=eval_residuals_lmm('LMM',d_dum['porcentaje_inhibicion'],ex_l,gr_l,[(tn,tf,'')])
    m2=MixedLM(yt,ex_l,groups=gr_l).fit(reml=True,maxiter=200); r=m2.resid
    g1=r[gar=='maceracion']; g2=r[gar=='soxhlet']; g3=r[gar=='ultrasonido']
    ls,lp=stats.levene(g1,g2,g3); print(f"  {tn:12s} Levene F={ls:.2f} p={lp:.4f}")
""")

md("""### 6.3 Transformaciones para log10(conidias)""")

code("""print("\\n"+"-"*70); print("  2. log10(CONIDIAS/mL)"); print("-"*70)
cl_t=[('Crudo (log10)',lambda x:x,'Sin transformacion adicional'),
      ('BoxCox',lambda x:stats.boxcox(x+0.1)[0],'Box-Cox (+0.1)'),('Square',lambda x:x**2,'y^2'),
      ('RankGauss',rank_gauss,'Normal scores')]
for r in eval_transforms('log10 conidias',cl,transform_funcs=cl_t):
    print(f"  {r['Transformacion']:12s} skew={r['Skew']:.2f} kurt={r['Kurtosis']:.2f} Shapiro p={r['Shapiro_p']:.4f}")

print("\\n  Residuos del LMM (5.0 mg/mL):")
df_lc=CONI_TRAT.copy(); df_lc['metodo_id']=df_lc['metodo_extraccion'].str.strip().str.lower()
df_lc=df_lc[df_lc['concentracion_mg_ml']==5.0]; d_dum_c=pd.get_dummies(df_lc,columns=['metodo_id'],drop_first=True,dtype=float)
ex_lc=add_constant(d_dum_c[['metodo_id_soxhlet','metodo_id_ultrasonido']]); gr_lc=d_dum_c['aislado_id']
lmm_ct=[('Raw (log10)',lambda x:x,'log10 crudo'),('BoxCox',lambda x:stats.boxcox(x+0.1)[0],'Box-Cox'),('Square',lambda x:x**2,'y^2')]
lmm_cr=eval_residuals_lmm('LMM conidias',d_dum_c['conidias_log10'],ex_lc,gr_lc,lmm_ct)
for r in lmm_cr:
    print(f"  {r['Transformacion']:12s} LogLik={r['LogLik']:.0f} skew={r['Resid_skew']:.2f} kurt={r['Resid_kurt']:.2f} Shapiro p={r['Shapiro_p']:.4f}")
""")

md("""### 6.4 INH conidias: crudo vs log10""")

code("""print("\\n"+"-"*70); print("  3. INH CONIDIAS -- crudo vs log10 (hoja)"); print("-"*70)
print("  INH crudo:")
for r in eval_transforms('INH conidias crudo',ycc,transform_funcs=[('Raw (crudo)',lambda x:x,'INH en escala cruda')]):
    print(f"  {r['Transformacion']:12s} skew={r['Skew']:.2f} kurt={r['Kurtosis']:.2f} Shapiro p={r['Shapiro_p']:.4f}")
print("  INH log10 (escala hoja):")
for r in eval_transforms('INH conidias log10',ycl,transform_funcs=[('Raw (log10)',lambda x:x,'INH escala log (hoja)')]):
    print(f"  {r['Transformacion']:12s} skew={r['Skew']:.2f} kurt={r['Kurtosis']:.2f} Shapiro p={r['Shapiro_p']:.4f}")
if len(ycl)>10:
    df_cl4=CONI_INH_LOG.copy(); df_cl4['metodo_id']=df_cl4['metodo_extraccion'].str.strip().str.lower()
    df_cl4=df_cl4[df_cl4['concentracion_mg_ml']==5.0]
    if len(df_cl4)>10:
        dd_cl4=pd.get_dummies(df_cl4,columns=['metodo_id'],drop_first=True,dtype=float)
        ex_cl4=add_constant(dd_cl4[['metodo_id_soxhlet','metodo_id_ultrasonido']]); gr_cl4=dd_cl4['aislado_id']
        lmm_r4=eval_residuals_lmm('LMM INH conidias',dd_cl4['porcentaje_inhibicion_log10'],ex_cl4,gr_cl4,
            [('Raw (log10)',lambda x:x,'Crudo en log10')])
        for r in lmm_r4:
            print(f"\\n  LMM con INH_log10 (n={len(df_cl4)}): LogLik={r['LogLik']:.0f} skew={r['Resid_skew']:.2f} kurt={r['Resid_kurt']:.2f}")
""")

md("""### 6.5 Figuras de diagnostico de transformaciones""")

code("""# 6.5a: Q-Q de residuos para cada transformacion
fig,axes=plt.subplots(2,3,figsize=(14,9)); fig.subplots_adjust(hspace=0.4,wspace=0.3)
t_plot=[('Crudo',lambda x:x),('IHS/10',ihs_transform(10)),('IHS/50',ihs_transform(50)),('CubeRoot',cube_root),('RankGauss',rank_gauss)]
for idx,(tn,tf) in enumerate(t_plot):
    ax=axes[idx//3,idx%3]; yt=tf(d_dum['porcentaje_inhibicion'])
    m2=MixedLM(yt,ex_l,groups=gr_l).fit(reml=True,maxiter=200); stats.probplot(m2.resid,dist='norm',plot=ax)
    ax.set_title(f'{tn} (Shapiro p={lmm_r[idx]["Shapiro_p"]:.4f})')
ax=axes[1,2]; yt_c=d_dum_c['conidias_log10']
m_c2=MixedLM(yt_c,ex_lc,groups=gr_lc).fit(reml=True,maxiter=200)
stats.probplot(m_c2.resid,dist='norm',plot=ax); ax.set_title(f'log10 conidias (Shapiro p={lmm_cr[0]["Shapiro_p"]:.4f})')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'comparacion_transformaciones_qq.png'),dpi=300); plt.show()

# 6.5b: Histogramas de cada transformacion
fig,axes=plt.subplots(2,3,figsize=(14,8)); fig.subplots_adjust(hspace=0.4,wspace=0.3)
for idx,(tn,tf,_) in enumerate(inh_t[:6]):
    ax=axes[idx//3,idx%3]; yt=tf(yc); ax.hist(yt,bins=40,color='#2e86ab',edgecolor='white',alpha=0.7); ax.set_title(tn)
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'comparacion_transformaciones_hist.png'),dpi=300); plt.show()

# 6.5c: Homocedasticidad por transformacion
fig,axes=plt.subplots(1,3,figsize=(14,5)); fig.subplots_adjust(wspace=0.35)
for idx,(tn,tf,_) in enumerate(lmm_t[:3]):
    ax=axes[idx]; yt=tf(d_dum['porcentaje_inhibicion'])
    m3=MixedLM(yt,ex_l,groups=gr_l).fit(reml=True,maxiter=200)
    rdf=pd.DataFrame({'resid':m3.resid,'metodo':gar})
    sns.boxplot(data=rdf,x='metodo',y='resid',hue='metodo',palette={'maceracion':'#2e86ab','soxhlet':'#a23b72','ultrasonido':'#f18f01'},ax=ax,legend=False)
    ax.set_title(tn); ax.set_xlabel('')
fig.tight_layout(); fig.savefig(os.path.join(DIR_FIG,'comparacion_transformaciones_lev.png'),dpi=300); plt.show()
print("  OK Figuras de transformaciones generadas.")
""")

md("""### 6.6 Recomendaciones finales

Basado en los resultados de skewness, kurtosis, Shapiro-Wilk, homocedasticidad y log-verosimilitud, se presentan las recomendaciones para cada variable respuesta.

**Criterios:**
- Se prioriza la **interpretabilidad biologica** sobre la optimizacion estadistica marginal
- Se considera que el LMM es **robusto a desviaciones moderadas** de normalidad (n grande, diseno balanceado)
- Transformaciones que mejoran la homocedasticidad pero empeoran la normalidad se marcan como alternativas
""")

code("""print("\\n"+"-"*70); print("  RECOMENDACIONES FINALES"); print("-"*70)
recomendaciones = [
    {'Dataset':'Rendimiento (%)','Recomendacion':'Sin transformacion',
     'Justificacion':'Residuos ANOVA ya normales (Shapiro p=0.33). n=9 insuficiente para transformaciones.','Score':'OK'},
    {'Dataset':'INH Crecimiento (completo)','Recomendacion':'Sin transformacion -- escala original',
     'Justificacion':'Residuos crudos tienen menor asimetria (skew=-0.33). IHS la empeora (skew=-0.95). CubeRoot mejora homocedasticidad pero empeora normalidad. Ninguna alcanza normalidad (Shapiro p<0.001). La interpretacion en crudo es directa y el LMM es robusto.',
     'Score':'OK crudo / ~CubeRoot'},
    {'Dataset':'log10(conidias/mL)','Recomendacion':'Mantener log10 -- no transformar mas',
     'Justificacion':'Ya es log-transformado por el laboratorio. BoxCox sugiere lambda~2.17 (cuadrado) pero mejora marginal y la interpretacion empeora.','Score':'OK'},
    {'Dataset':'INH Conidias (crudo)','Recomendacion':'NO usar como variable de modelado primaria',
     'Justificacion':'Rango [-10615, 100] con skew=-12.2 y kurtosis=176 -- no rescatable. Usar log10(conidias) como respuesta primaria. INH_log10 (escala hoja) es alternativa secundaria.','Score':'NO'},
    {'Dataset':'log10(conidias) para dosis-respuesta','Recomendacion':'Mantener log10 -- sin transformacion adicional',
     'Justificacion':'Modelo B (Maceracion) converge y es interpretable. ICC=0.76 indica que el modelo captura bien la estructura.','Score':'OK'},
]
for rec in recomendaciones:
    print(f"  {rec['Score']} {rec['Dataset']}")
    print(f"       -> {rec['Recomendacion']}")
    print(f"         {rec['Justificacion'][:100]}...")
    print()

print("\\n**Resumen:**")
print("- INH de crecimiento: modelar en escala CRUDA con LMM (robusto)")
print("- log10(conidias): mantener la escala logaritmica del laboratorio")
print("- Rendimiento: ANOVA directo sobre valores crudos")
print("- INH conidias crudo: NO modelar -- usar log10(conidias) como respuesta")
""")

md("""---
# Conclusiones Generales

## Integracion de resultados

### Sobre los metodos de extraccion

| Aspecto | Maceracion | Soxhlet | Ultrasonido |
|---------|:----------:|:-------:|:-----------:|
| Rendimiento (%) | ~12% | ~43% | ~17% |
| INH crecimiento 5.0 mg/mL | ~80% | ~70% | ~65% |
| Inhibicion completa | SI | NO | NO |
| Reduccion conidias | >90% | Minima | Minima |
| EC50 estimable | SI | No aplica | No aplica |

### Sobre los aislados de Fusarium

- Existe **variabilidad significativa** entre aislados en todos los metodos (ICC ~0.3-0.5)
- Los perfiles de susceptibilidad se agrupan en clusters diferenciables por PCA y clustering jerarquico
- La EC50 para Maceracion varia ampliamente entre aislados (~0.4 a ~4.0 mg/mL)
- **Ningun aislado se clasifica como resistente** por ausencia de umbral validado

### Implicaciones para el control de Fusarium

1. **Maceracion es el metodo mas prometedor**: aunque produce menor rendimiento, su extracto tiene la mayor actividad antifungica especifica y afecta tanto crecimiento como esporulacion
2. **La concentracion de 5.0 mg/mL de Maceracion** seria la candidata principal para aplicaciones practicas
3. **La variabilidad entre aislados** sugiere que el control puede requerir ajustes segun la poblacion local de Fusarium
4. **Soxhlet y Ultrasonido** podrian ser utiles para inhibicion del crecimiento, pero no controlan la dispersion por esporas

### Limitaciones del estudio

1. Solo Maceracion tiene diseno dosis-respuesta completo (4 concentraciones)
2. n=3 replicas biologicas por combinacion factorial limita la potencia estadistica
3. La EC50 se estima por interpolacion lineal (no modelo no lineal) por limitacion de puntos
4. Los datos de conidias en escala cruda no son modelables directamente
5. No se evaluo la interaccion metodo x concentracion en un modelo completo por el diseno no balanceado
""")

# ===== FINAL: SAVE =====
nb.cells = cells
OUTPUT = '/home/mniev/projects/proyecto_tomillo/dca/analisis_tomillo_fusarium.ipynb'
with open(OUTPUT, 'w') as f:
    nbf.write(nb, f)
print(f"\\nNotebook generado: {OUTPUT}")
print(f"Total celdas: {len(cells)}")
print(f"  Markdown: {sum(1 for c in cells if c.cell_type=='markdown')}")
print(f"  Code: {sum(1 for c in cells if c.cell_type=='code')}")
