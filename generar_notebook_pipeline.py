#!/usr/bin/env python3
"""Genera el notebook orquestador ``dca/analisis_dca.ipynb``.

Construye un notebook autocontenido y educativo con nbformat (kernel python3)
que orquesta las 12 fases de la pipeline reproducible. El código pesado vive
en ``pipeline/`` (módulos importables); aquí se definen las celdas markdown
(explicativas, en español profesional y neutro) y las celdas de código que
llaman a las funciones de cada fase.

Uso:
    python3 generar_notebook_pipeline.py
"""

from __future__ import annotations

from pathlib import Path

import nbformat

RAIZ = Path(__file__).resolve().parent
SALIDA = RAIZ / "dca" / "analisis_dca.ipynb"


# ---------------------------------------------------------------------------
# Celdas markdown (explicativas, en español profesional y neutro)
# ---------------------------------------------------------------------------


def _md_titulo() -> str:
    return """# Pipeline reproducible de análisis estadístico

## Actividad antifúngica de extractos de tomillo (*Thymus vulgaris*) frente a *Fusarium* spp.

Este notebook es el **orquestador educativo** de una pipeline estadística
automatizada, reproducible y reutilizable. Todo el código de cómputo pesado
vive en el paquete `pipeline/` (importable); aquí se explica cada fase, se
ejecuta y se interpretan los resultados.

### Contexto del proyecto

Se evaluó la actividad antifúngica de extractos de tomillo obtenidos por tres
técnicas de extracción (**maceración**, **Soxhlet** y **ultrasonido**) frente
a **31 aislados de *Fusarium* spp.**, todos ensayados a una única concentración
de **5 mg/mL**. Para cada combinación técnica × aislado se dispone de **3
réplicas biológicas** (unidad experimental: caja Petri), totalizando 279
observaciones de bioensayo, más 9 mediciones de rendimiento de extracción.

### Objetivos científicos

1. **Rendimiento de extracción**: determinar si la técnica afecta el rendimiento (%).
2. **Inhibición del crecimiento micelial**: comparar la actividad antifúngica
   de las técnicas y de los aislados.
3. **Inhibición de la producción de conidias**: evaluar el efecto sobre la
   esporulación (variable continua en escala log10).
4. **Susceptibilidad relativa de los aislados**: agrupar los aislados según su
   perfil de susceptibilidad (NUNCA se usa el término "resistente" sin un
   criterio validado).
5. **Ranking de técnicas**: integrar rendimiento y actividad en un score.
6. **Reproducibilidad**: que todo el análisis sea re-ejecutable con datos nuevos.

### Estructura de las 12 fases

| Fase | Módulo | Contenido |
|------|--------|-----------|
| 1-2 | `cargar_datos` | Carga de datos y auditoría de calidad |
| 3 | `limpiar` | Normalización del dataset maestro |
| 4 | `eda` | Exploración y descriptivos |
| 5 | `diseno` | Inferencia del diseño experimental |
| 6 | `supuestos` | Verificación de supuestos de los modelos |
| 7 | `modelos` | Análisis inferencial y selección automática |
| 8 | `comparaciones` | Comparaciones múltiples y letras CLD |
| 9 | `visualizar` | Figuras de resultados |
| 10 | `multivariado` | PCA, clustering y categorías biológicas |
| 11 | `ranking` | Ranking de técnicas |
| 12 | `informe` | Informe final (Markdown y HTML) |

> **Contrato científico**: se respetan las reglas del proyecto: no se eliminan
> atípicos automáticamente, no se llama "resistente" a ningún aislado, se
> reportan efecto e IC además del p-valor, y se distingue significancia
> estadística de significancia biológica.
"""


def _md_setup() -> str:
    return """## Configuración del entorno

Esta celda prepara el entorno: determina la raíz del proyecto, agrega el
paquete `pipeline/` al `sys.path`, fija la semilla aleatoria global (42) para
garantizar la reproducibilidad de los procedimientos estocásticos (KMeans,
PCA, etc.) y verifica las versiones de las librerías utilizadas.

**Por qué una semilla fija**: el análisis de clusters y cualquier procedimiento
con inicialización aleatoria producen resultados ligeramente distintos entre
ejecuciones. Fijar `random_state`/`seed` hace que el análisis sea determinista.
"""


def _md_fase_1() -> str:
    return """## Fase 1: Carga de datos

**Qué se hace**: se leen las dos fuentes de datos de la pipeline:
- `dca/resultados/database/consolidado_tidy.xlsx` (hoja "Consolidado", 279 filas):
  bioensayo con crecimiento micelial (mm), % de inhibición micelial, conidias
  (log10/mL) e inhibición de conidias para cada técnica × aislado × réplica.
- `dca/resultados/database/rendimiento_extraccion.csv` (9 filas): rendimiento de
  extracción por técnica y réplica biológica.

**Por qué**: se usa el dataset consolidado y limpio (no el Excel crudo) como
fuente única de entrada; esto garantiza trazabilidad.

**Cómo interpretar**: se verifica que las dimensiones coincidan con lo esperado
(279 y 9 filas) y se inspeccionan columnas, tipos y primeras filas antes de
cualquier transformación.
"""


def _md_fase_2() -> str:
    return """## Fase 2: Auditoría de calidad de datos

**Qué se hace**: se audita el dataset columna a columna: valores faltantes,
filas duplicadas exactas, columnas constantes, tipos, valores inconsistentes
(método en {Maceración, Soxhlet, Ultrasonido}, réplica en {R1,R2,R3}, %INH
micelial en [0,100], crecimiento en [0,100], conidias >= 0) y atípicos por
rango intercuartílico (1.5xIQR).

**Por qué**: la auditoría previa es obligatoria para detectar duplicados,
faltantes y valores imposibles antes del análisis (regla de integridad del
proyecto).

**Supuestos**: se asume que los límites de cada variable son los documentados
en el diccionario de datos; el %INH de conidias puede ser negativo por
diseño (el extracto puede inducir mayor esporulación) y por eso no se
restringe su rango.

**Cómo interpretar**: los atípicos se **flaguean pero nunca se eliminan**;
su existencia se documenta y su influencia se evalúa en los análisis.
Valores "inconsistentes" > 0 indicarían errores de registro que requieren
revisión con el investigador.
"""


def _md_fase_3() -> str:
    return """## Fase 3: Limpieza y normalización del dataset maestro

**Qué se hace**: se renombran las columnas al esquema canónico en inglés
(`metodo_extraccion`, `aislamiento`, `replica`, `crecimiento_micelial_mm`,
`porcentaje_inhibicion_micelial`, `conidias_log10_ml`,
`porcentaje_inhibicion_conidias`), se normaliza el método a minúsculas sin
acentos y la réplica a entero 1-3. Se **valida** que el dataset cumpla: 279
filas, sin valores nulos y diseño balanceado (31 aislados × 3 métodos × 3
réplicas). Luego se incorporan los **controles C4** del Excel crudo del
laboratorio (columnas `control_crecimiento_mm` y `control_conidias_log10`,
una por aislado y compartidas por sus 3 réplicas), dejando el **dataset
maestro** con 9 columnas. Se guarda en CSV y XLSX (hojas Bioensayo y
Rendimiento) junto con el diccionario de datos.

**Por qué**: un esquema canónico estable facilita reutilizar la pipeline con
nuevos archivos y evita errores de codificación (acentos, mayúsculas); los
controles C4 explícitos permiten validar la fórmula del %INH (integridad de
datos) y sirven de línea de base.

**Supuestos**: el dataset consolidado ya contiene los %INH calculados contra el
control C4 de cada aislado; la concentración es constante (5 mg/mL) y no se
modela. Si el Excel crudo no está disponible, los controles no se incorporan
y la validación del %INH se saltea con un aviso (sin romper el flujo).

**Cómo interpretar**: si la validación lanzara un error, habría que resolver el
problema de datos antes de continuar; un diseño desbalanceado cambiaría la
selección de modelos.
"""


def _md_validacion_inh() -> str:
    return """## Fase 3.5: Validación de la fórmula del %INH (integridad de datos)

**Qué se hace**: se reconstruye el %INH con la fórmula
`%INH = (1 - C1/C4) × 100` usando las columnas de control C4 incorporadas en la
Fase 3 y se compara contra los valores reportados por el laboratorio
(`porcentaje_inhibicion_micelial` y `porcentaje_inhibicion_conidias`), fila por
fila. Se reporta la máxima diferencia absoluta, el número de discrepancias y el
estado.

**Por qué**: es una verificación de integridad de datos: confirma que los %INH
reportados son consistentes con la fórmula documentada y con los controles.
**No reemplaza** las respuestas reportadas: la inferencia principal usa el
%INH reportado por el investigador (instrucción explícita).

**Supuestos**: para las conidias, la fórmula se aplica sobre la escala log10
directamente, porque el laboratorio reportó la reducción de conidias en log10
(no en conteos crudos). Si un control fuera 0, el %INH verificado se define
como 100 cuando C1=0 y como indefinido (NaN) cuando C1>0; en los datos reales
el control nunca es 0.

**Cómo interpretar**: si el estado es 'ok' (sin discrepancias por encima de la
tolerancia 1e-6), la fórmula del laboratorio es consistente con los datos; si
hubiera discrepancias, habría que investigar si se deben a escala (log10 vs
crudo) o a redondeo antes de reportar.
"""


def _md_fase_4() -> str:
    return """## Fase 4: Análisis exploratorio (EDA)

**Qué se hace**: se calculan descriptivos por método (n, media, DE, error
estándar, IC95%, mínimo y máximo) para las cuatro variables de respuesta, y se
generan figuras exploratorias: histogramas, boxplots, violin plots, densidades
con rug, QQ-plots, matriz de correlación (Pearson con p) y un scatter entre
crecimiento y conidias.

**Por qué**: el EDA revela la forma de las distribuciones (asimetrías, efecto
techo en %INH micelial), la variabilidad entre métodos y las relaciones entre
respuestas, orientando la elección del modelo.

**Supuestos**: ninguno inferencial; es una fase descriptiva.

**Cómo interpretar**:
- Si el %INH micelial se concentra en 100 (efecto techo), la distribución es
  asimétrica y los métodos no paramétricos o modelos mixtos ganan relevancia.
- Si el %INH de conidias toma valores negativos, el extracto indujo mayor
  esporulación en esas celdas (la escala es log10).
- Correlaciones altas entre respuestas sugieren que la inhibición micelial y
  la de conidias covarían (a revisar en el multivariado).
"""


def _md_fase_5() -> str:
    return """## Fase 5: Diseño experimental

**Qué se hace**: se infiere y documenta el diseño: DCA factorial
técnica × aislado, 3 réplicas biológicas, unidad experimental (caja Petri),
concentración única (5 mg/mL) y balanceo.

**Por qué**: conocer el diseño determina la estructura del modelo (factores
fijos/aleatorios, interacciones) y las limitaciones de la inferencia.

**Supuestos y caveats**:
- El aislado se trata como factor **fijo** en el ANOVA factorial y como
  **aleatorio** en el análisis de sensibilidad LMM.
- Cada %INH se calculó contra **un único control C4 compartido** por las tres
  réplicas del aislado: las réplicas de %INH no son totalmente independientes
  (pseudorreplicación del control). El crecimiento crudo en mm no presenta
  este problema.

**Cómo interpretar**: la tabla resume el diseño inferido; el texto resalta las
limitaciones que condicionarán las conclusiones (especialmente para %INH).
"""


def _md_fase_6() -> str:
    return """## Fase 6: Verificación de supuestos

**Qué se hace**: para cada variable de respuesta se ajusta el modelo OLS
factorial `variable ~ método * aislado` y se evalúa:
- **Normalidad** de residuos (Shapiro-Wilk),
- **Homocedasticidad** entre métodos (Levene y Bartlett),
- **Independencia** (estadístico de Durbin-Watson),
y se generan figuras de diagnóstico (histograma, QQ-plot, residuos vs
ajustados, residuos por método).

**Por qué**: el ANOVA factorial exige residuos normales, homocedásticos e
independientes; si no se cumplen, la inferencia debe apoyarse en la vía no
paramétrica (fase 7) o en modelos mixtos.

**Supuestos**: p > 0.05 se interpreta como ausencia de evidencia en contra del
supuesto (no como prueba de su cumplimiento).

**Cómo interpretar**:
- Shapiro p < 0.05 con %INH micelial refleja el efecto techo (muchos valores
  en 100).
- Bartlett puede no ser calculable si un grupo es constante (varianza nula),
  por ejemplo en crecimiento 0 (inhibición completa); se reporta como tal.
- Durbin-Watson cercano a 2 sugiere independencia de los residuos.
"""


def _md_fase_7() -> str:
    return """## Fase 7: Análisis inferencial y selección automática de modelos

**Qué se hace**:
1. **Rendimiento**: ANOVA de una vía (OLS) `rendimiento_pct ~ método` con
   tabla tipo II, eta², omega² y supuestos; si fallan, se complementa con
   Kruskal-Wallis.
2. **Factorial por variable**: se ajusta `variable ~ método * aislado` con
   ANOVA tipo II y tamaños de efecto (eta² y omega² parciales). La selección
   de la vía de inferencia es **automática y justificada**:
   - Si la variable es un **conteo entero >= 0 sobredisperso** -> rama GLM
     Poisson/Binomial negativa (función `glm_conteos`, documentada).
   - Si es continua y los supuestos se cumplen -> tabla F del ANOVA.
   - Si los supuestos fallan -> vía no paramétrica: Kruskal-Wallis por método
     + Scheirer-Ray-Hare (ANOVA tipo II sobre rangos) para la interacción.
3. **Sensibilidad LMM**: `variable ~ método + (1|aislamiento)` (modelo mixto
   con aislado aleatorio) para evaluar si la conclusión sobre el método es
   robusta; se reporta el ICC (proporción de varianza atribuible al aislado).
4. **Conidias**: diagnóstico de que `conidias_log10_ml` es continua (no
   entera), lo que **desactiva la rama Poisson/NB** (los conteos crudos no
   están disponibles y la escala ya es log10); se usa modelo lineal sobre
   log10.

**Por qué**: la filosofía del proyecto es NO elegir un test solo por producir
significancia; la elección se justifica con diagnóstico de datos.

**Supuestos**: OLS requiere residuos normales, homocedásticos e independientes;
el LMM requiere normalidad de los efectos aleatorios (aproximada).

**Cómo interpretar**:
- Reportar **F, p, eta², omega² e IC**; un p significativo con eta² pequeño
  no es biológicamente relevante.
- Un ICC alto indica que el aislado explica gran parte de la variación.
- La comparación entre el factorial y el LMM permite evaluar robustez.
"""


def _md_fase_8() -> str:
    return """## Fase 8: Comparaciones múltiples

**Qué se hace**: se comparan los métodos por pares según el modelo seleccionado
en la fase 7:
- Si el modelo fue paramétrico (ANOVA) -> **Tukey HSD** entre métodos.
- Si fue no paramétrico -> **test de Dunn** (manual, sobre rangos) con
  corrección **FDR** (Benjamini-Hochberg) y **Wilcoxon/Mann-Whitney** con FDR
  como robustez.
Se generan **letras compactas (CLD)**: métodos que comparten al menos una letra
NO difieren significativamente (p ajustada >= 0.05). Se guardan tablas y la
figura `posthoc_<variable>_letras`.

**Por qué**: comparar todas las técnicas por pares sin corregir infla el error
tipo I; las letras CLD facilitan la lectura de los grupos homogéneos.

**Supuestos**: el post-hoc debe corresponder al modelo principal (regla del
proyecto). Cuando la interacción método × aislado es significativa, la
comparación marginal entre métodos debe interpretarse con cautela: describe el
efecto promedio sobre los aislados.

**Cómo interpretar**: en cada variable, los métodos que comparten letra son
estadísticamente equivalentes al nivel 0.05 (ajustado); el orden de las medias
aporta la dirección del efecto.
"""


def _md_fase_9() -> str:
    return """## Fase 9: Visualización de resultados

**Qué se hace**: se generan figuras de resultados (prefijo `resultados_`) para
cada variable: medias por método con SD / SE / IC95%, gráfico de interacción
método × aislado (medias por celda), efectos principales (método y aislado) y
comparación con letras CLD. Se guarda la tabla `medias_<variable>.csv`.

**Por qué**: las figuras de calidad de publicación permiten comunicar los
patrones de forma directa y verificable.

**Cómo interpretar**:
- Barras separadas en la figura CLD = grupos distintos.
- Un gráfico de interacción con líneas cruzadas sugiere que el efecto del
  método depende del aislado (respalda el término de interacción del ANOVA).
- Los IC95% que no se solapan indican diferencias descriptivas entre medias.
"""


def _md_fase_10() -> str:
    return """## Fase 10: Análisis multivariado de susceptibilidad

**Qué se hace**:
1. Se construye la **matriz por aislado** (31 x 6): media de %INH micelial y
   %INH de conidias para cada técnica; se estandariza (z-score).
2. **PCA**: varianza explicada, scree plot y biplot PC1-PC2.
3. **Clustering jerárquico** de Ward con dendrograma y coeficiente cofenético.
4. **KMeans** con codo (inercia) y silhouette; k óptimo automático.
5. **Categorías biológicas**: score compuesto de susceptibilidad (promedio de
   los z de inhibición micelial y de conidias) -> terciles -> etiquetas
   "Alta / Moderada / Baja susceptibilidad relativa". Se cruzan con los
   clusters de KMeans.
6. **Heatmap** de las 6 métricas estandarizadas con dendrograma y anotaciones
   de categoría; scatter de clusters en PC1-PC2.

**Por qué**: la susceptibilidad es un perfil multidimensional; PCA y clustering
la resumen objetivamente.

**Supuestos**: estandarización previa (las escalas de % difieren en
variabilidad); el número de clusters se elige por criterio objetivo
(silhouette).

**Cómo interpretar**:
- Mayor score compuesto = mayor susceptibilidad relativa (mayor inhibición).
- **No se usa el término "resistente"**: sin un umbral biológico validado, se
  habla de susceptibilidad relativa.
- Clusters con pocos aislados deben interpretarse con cautela (tamaño de
  muestra pequeño).
"""


def _md_fase_11() -> str:
    return """## Fase 11: Ranking de técnicas de extracción

**Qué se hace**: por técnica se calculan tres métricas (rendimiento medio,
%INH micelial medio, %INH de conidias medio), se normalizan min-max a 0-1
(1 = mejor) y se promedian en un **score compuesto**; se ordenan las técnicas
y se genera un gráfico radar de tres ejes.

**Por qué**: integra rendimiento y actividad antifúngica en una única medida
de desempeño global, útil para la toma de decisiones.

**Supuestos**: las tres métricas se ponderan por igual (promedio simple); el
score es una herramienta descriptiva, no inferencial.

**Cómo interpretar**: la técnica con mayor score compuesto combina buen
rendimiento y alta actividad; verificar que las diferencias no contradigan las
comparaciones estadísticas de las fases 7-8 (una técnica puede liderar el
ranking pero no diferir significativamente de la segunda).
"""


def _md_fase_12() -> str:
    return """## Fase 12: Informe final

**Qué se hace**: se genera `dca/resultados/reportes/informe_final.md` (español
profesional y neutro) con las secciones: resumen ejecutivo, calidad de datos,
diseño experimental, supuestos, análisis seleccionados y su justificación,
comparaciones múltiples, análisis multivariado, ranking, interpretación
biológica, conclusiones y limitaciones. Luego se convierte a
`informe_final.html` con un CSS simple y legible. Además se exporta un libro
Excel de resumen (`dca/resultados/excel/resumen_analisis.xlsx`) con las hojas
Descriptivos, Rendimiento (ANOVA), Factorial, Posthoc, Ranking y
Susceptibilidad.

**Por qué**: consolida todos los resultados en un documento accionable y
compartible, siguiendo la filosofía de reportar efecto, IC y p-valor.
"""


def _md_conclusiones() -> str:
    return """## Conclusiones generales

1. **Rendimiento**: la técnica afectó significativamente el rendimiento de
   extracción; los resultados se reportan con tamaño de efecto (eta², omega²).
2. **Actividad antifúngica**: a 5 mg/mL los extractos inhibieron el crecimiento
   micelial con un marcado efecto techo en maceración; las comparaciones entre
   técnicas y aislados se apoyaron en la vía seleccionada por los supuestos.
3. **Conidias**: la variable es continua en escala log10; se modeló con un
   modelo lineal sobre log10 (la rama Poisson/NB no aplica y quedó documentada).
4. **Susceptibilidad**: los aislados se clasificaron en Alta / Moderada / Baja
   susceptibilidad relativa mediante score compuesto y clustering.
5. **Ranking**: la tabla de la fase 11 resume el desempeño global de las
   técnicas.

## Cómo re-ejecutar con un archivo nuevo

La pipeline está diseñada para reutilizarse sin modificar la lógica:

1. **Reemplazar las fuentes** editando las rutas en `pipeline/config.py`
   (`EXCEL_TIDY` y `CSV_RENDIMIENTO`) o, mejor, pasando rutas alternativas a
   `cargar_datos()` (acepta paths parametrizables).
2. **Respetar el esquema de columnas** del consolidado (7 columnas con los
   mismos nombres) o ajustar `RENOMBRES_BIO` en `pipeline/limpiar.py`. El
   master tendrá 7 columnas de respuesta más las 2 de control C4 cuando el
   Excel crudo del laboratorio esté disponible (9 columnas en total).
3. **Volver a generar y ejecutar**: `python3 generar_notebook_pipeline.py` y
   ejecutar el notebook completo.
4. **Revisar la validación** de la fase 3: si las dimensiones o el balanceo
   cambian, la validación lanzará un error y habrá que decidir cómo proceder
   (el pipeline nunca imputa datos ni elimina atípicos automáticamente).
5. **Controles y validación**: si el Excel crudo
   (`datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx`) no está presente,
   los controles C4 se omiten y la validación del %INH se saltea con un aviso;
   la inferencia no se ve afectada.

Los resultados se escriben siempre en `dca/resultados/` (tablas, figuras,
reportes y Excel), sin tocar los datos fuente.
"""


# ---------------------------------------------------------------------------
# Celdas de código
# ---------------------------------------------------------------------------


def _code_setup() -> str:
    return '''
import os
import sys
import warnings
from pathlib import Path

# Determinar la raíz del proyecto (directorio que contiene pipeline/)
RAIZ = Path(os.getcwd()).resolve()
for candidato in (RAIZ, RAIZ.parent, RAIZ.parent.parent):
    if (candidato / "pipeline" / "config.py").exists():
        RAIZ = candidato
        break
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))
os.chdir(RAIZ)
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from IPython.display import display, Markdown

from pipeline.config import (
    fijar_semilla, METODOS, METODO_LABEL, VARIABLES_RESPUESTA, VARIABLE_LABEL,
)
from pipeline import (
    cargar_datos, limpiar, eda, diseno, supuestos, modelos, comparaciones,
    visualizar, multivariado, ranking, informe,
)

fijar_semilla(42)

print("Versión de librerías:")
import importlib
for nombre in ("pandas", "numpy", "scipy", "statsmodels", "sklearn", "pingouin", "matplotlib", "seaborn"):
    try:
        mod = importlib.import_module(nombre)
        print(f"  {nombre} {mod.__version__}")
    except Exception as exc:  # pragma: no cover
        print(f"  {nombre}: {exc}")

# Diccionario acumulador de resultados para el informe final
resultados = {}
print("\\nRaíz del proyecto:", RAIZ)
'''


def _code_fase_1() -> str:
    return '''
datos = cargar_datos.cargar_datos()
resultados["datos"] = datos
'''


def _code_fase_2() -> str:
    return '''
tabla_auditoria, resumen_auditoria = cargar_datos.auditoria_calidad(datos["bio"])
resultados["auditoria"] = {"tabla": tabla_auditoria, "resumen": resumen_auditoria}
display(tabla_auditoria)
'''


def _code_fase_3() -> str:
    return '''
df_bio, df_rend = limpiar.normalizar_master(datos["bio"], datos["rend"])
limpiar.guardar_master(df_bio, df_rend)
resultados["master"] = {"bio": df_bio, "rend": df_rend}
display(df_bio.head(8))
'''


def _code_validacion_inh() -> str:
    return '''
validacion_inh = limpiar.validar_inh(df_bio)
resultados["validacion_inh"] = validacion_inh
'''


def _code_fase_4() -> str:
    return '''
tabla_desc = eda.resumen_descriptivo(df_bio)
resultados["eda"] = {"descriptivos": tabla_desc}
display(tabla_desc)

figuras_eda = eda.figuras_eda(df_bio)
resultados["eda"]["figuras"] = figuras_eda
'''


def _code_fase_5() -> str:
    return '''
diseno_info = diseno.inferir_diseno(df_bio)
resultados["diseno"] = diseno_info
display(diseno_info["detalle"])
'''


def _code_fase_6() -> str:
    return '''
supuestos_result = {}
for variable in VARIABLES_RESPUESTA:
    supuestos_result[variable] = supuestos.verificar_supuestos(df_bio, variable)
    display(supuestos_result[variable]["tabla_supuestos"])
resultados["supuestos"] = supuestos_result
'''


def _code_fase_7() -> str:
    return '''
# 7.1 Rendimiento de extracción
an_rendimiento = modelos.analisis_rendimiento(df_rend)
resultados["rendimiento"] = an_rendimiento
display(an_rendimiento["tabla_anova"])

# 7.2 Factorial por variable de respuesta
modelos_res = {}
lmm_res = {}
for variable in VARIABLES_RESPUESTA:
    m = modelos.analisis_factorial(df_bio, variable)
    modelos_res[variable] = m
    display(m["tabla_anova"])

# 7.3 Sensibilidad con modelo mixto (aislado aleatorio)
for variable in VARIABLES_RESPUESTA:
    lmm_res[variable] = modelos.analisis_sensibilidad_lmm(df_bio, variable)
    display(lmm_res[variable]["tabla_efectos_fijos"])

# 7.4 Diagnóstico de conidias
conidias_info = modelos.analisis_conidias(df_bio)
display(conidias_info["diagnostico"])

resultados["modelos"] = modelos_res
resultados["lmm"] = lmm_res
resultados["conidias"] = conidias_info
'''


def _code_fase_8() -> str:
    return '''
posthoc_result = {}
for variable in VARIABLES_RESPUESTA:
    tipo_modelo = modelos_res[variable]["tipo_modelo"]
    posthoc_result[variable] = comparaciones.comparaciones_posthoc(df_bio, variable, tipo_modelo)
    display(posthoc_result[variable]["tabla_pares"])
resultados["posthoc"] = posthoc_result
'''


def _code_fase_9() -> str:
    return '''
figuras_resultado = visualizar.figuras_resultados(df_bio, posthoc_result)
resultados["figuras_resultados"] = figuras_resultado
'''


def _code_fase_10() -> str:
    return '''
multivariado_info = multivariado.analisis_multivariado(df_bio)
resultados["multivariado"] = multivariado_info
display(multivariado_info["tabla_final"].head(15))
display(multivariado_info["cruce"])
'''


def _code_fase_11() -> str:
    return '''
ranking_info = ranking.ranking_tecnicas(df_bio, df_rend)
resultados["ranking"] = ranking_info
display(ranking_info["tabla"])
'''


def _code_fase_12() -> str:
    return '''
informe.generar_informe(resultados)

# Libro Excel de resumen con las hojas solicitadas
from pipeline.config import exportar_excel, DIR_EXCEL

hojas_resumen = {
    "Descriptivos": tabla_desc,
    "Rendimiento_ANOVA": an_rendimiento["tabla_anova"],
    "Factorial_INH_micelial": modelos_res["porcentaje_inhibicion_micelial"]["tabla_anova"],
    "Posthoc_INH_micelial": posthoc_result["porcentaje_inhibicion_micelial"]["tabla_pares"],
    "Ranking": ranking_info["tabla"],
    "Susceptibilidad": multivariado_info["tabla_final"],
}
ruta_excel = exportar_excel(hojas_resumen, DIR_EXCEL / "resumen_analisis.xlsx")
print("Libro Excel de resumen:", ruta_excel)
'''


def _md_pie() -> str:
    return """### Nota final sobre interpretación

- **Significancia estadística** (p < 0.05) no implica **relevancia biológica**;
  siempre se evalúan junto con el tamaño de efecto y los IC95%.
- Las limitaciones principales (control compartido, %INH de conidias en escala
  log10 y ausencia de dosis-respuesta) están documentadas en la sección 10 del
  informe final.
- Todos los archivos generados quedan en `dca/resultados/`; el notebook solo
  orquesta la ejecución.
"""


# ---------------------------------------------------------------------------
# Construcción del notebook
# ---------------------------------------------------------------------------


def construir_notebook() -> nbformat.NotebookNode:
    """Construye el notebook con nbformat y lo devuelve."""
    nb = nbformat.v4.new_notebook()
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    }
    nb.cells = [
        nbformat.v4.new_markdown_cell(_md_titulo()),
        nbformat.v4.new_markdown_cell(_md_setup()),
        nbformat.v4.new_code_cell(_code_setup()),
        nbformat.v4.new_markdown_cell(_md_fase_1()),
        nbformat.v4.new_code_cell(_code_fase_1()),
        nbformat.v4.new_markdown_cell(_md_fase_2()),
        nbformat.v4.new_code_cell(_code_fase_2()),
        nbformat.v4.new_markdown_cell(_md_fase_3()),
        nbformat.v4.new_code_cell(_code_fase_3()),
        nbformat.v4.new_markdown_cell(_md_validacion_inh()),
        nbformat.v4.new_code_cell(_code_validacion_inh()),
        nbformat.v4.new_markdown_cell(_md_fase_4()),
        nbformat.v4.new_code_cell(_code_fase_4()),
        nbformat.v4.new_markdown_cell(_md_fase_5()),
        nbformat.v4.new_code_cell(_code_fase_5()),
        nbformat.v4.new_markdown_cell(_md_fase_6()),
        nbformat.v4.new_code_cell(_code_fase_6()),
        nbformat.v4.new_markdown_cell(_md_fase_7()),
        nbformat.v4.new_code_cell(_code_fase_7()),
        nbformat.v4.new_markdown_cell(_md_fase_8()),
        nbformat.v4.new_code_cell(_code_fase_8()),
        nbformat.v4.new_markdown_cell(_md_fase_9()),
        nbformat.v4.new_code_cell(_code_fase_9()),
        nbformat.v4.new_markdown_cell(_md_fase_10()),
        nbformat.v4.new_code_cell(_code_fase_10()),
        nbformat.v4.new_markdown_cell(_md_fase_11()),
        nbformat.v4.new_code_cell(_code_fase_11()),
        nbformat.v4.new_markdown_cell(_md_fase_12()),
        nbformat.v4.new_code_cell(_code_fase_12()),
        nbformat.v4.new_markdown_cell(_md_conclusiones()),
        nbformat.v4.new_markdown_cell(_md_pie()),
    ]
    return nb


def main() -> None:
    """Construye y guarda el notebook, e imprime el número de celdas."""
    nb = construir_notebook()
    nbformat.write(nb, SALIDA)
    n_codigo = sum(1 for c in nb.cells if c.cell_type == "code")
    n_markdown = sum(1 for c in nb.cells if c.cell_type == "markdown")
    print(f"Notebook generado: {SALIDA}")
    print(f"Número de celdas: {len(nb.cells)} ({n_markdown} markdown, {n_codigo} codigo)")


if __name__ == "__main__":
    main()
