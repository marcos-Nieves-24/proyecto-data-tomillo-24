# Análisis factorial: genotipo de trigo × cepa de Fusarium (Snijders)

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
| F39 | 24.30 | 2.83 | [18.61, 30.00] |
| F436 | 10.41 | 1.37 | [7.65, 13.17] |
| F329 | 6.36 | 0.95 | [4.46, 8.27] |
| F348 | 5.72 | 0.96 | [3.79, 7.66] |

### 3.2 Por genotipo (10 de mayor severidad)

| Genotipo | Media | SE | IC95 |
|---|---|---|---|
| SVP72005-20-3-1 | 26.06 | 5.41 | [14.16, 37.96] |
| SVP75059-32 | 24.83 | 4.94 | [13.95, 35.72] |
| Nautica | 22.99 | 5.03 | [11.93, 34.06] |
| SVP75059-46 | 18.98 | 5.29 | [7.33, 30.64] |
| SVP73012-1-2-3 | 18.56 | 5.14 | [7.25, 29.87] |
| SVP73016-2-4 | 16.09 | 4.34 | [6.53, 25.65] |
| SVP73030-8-1-I | 13.64 | 4.93 | [2.79, 24.49] |
| SVP75059-28 | 9.08 | 4.36 | [-0.52, 18.68] |
| SVP77079-15 | 9.03 | 3.13 | [2.15, 15.90] |
| SVP72003-4-2-4 | 7.63 | 2.03 | [3.16, 12.11] |
| SVP77078-30 | 6.81 | 1.73 | [2.99, 10.62] |
| Saiga | 6.73 | 1.32 | [3.82, 9.65] |
| SVP77076-1 | 4.23 | 0.99 | [2.06, 6.41] |
| SVP77076-38 | 3.97 | 1.50 | [0.67, 7.28] |
| Arina | 3.94 | 1.00 | [1.74, 6.15] |
| SVP77076-4 | 3.69 | 1.20 | [1.06, 6.32] |
| SVP72017-17-5-10 | 2.62 | 0.55 | [1.40, 3.84] |

## 4. Modelo estadístico

Se ajustó un modelo lineal `y ~ C(gen) + C(strain) + C(year) + C(gen):C(strain)`. El factor `year` se incorpora como factor de bloqueo/replicación: dado que existe una sola observación por combinación genotipo × cepa, los tres años proporcionan la replicación necesaria para estimar el error, pero el residual confunde la variación entre años con el error puro (no existe término de error puro dentro de la celda genotipo × cepa). La interacción genotipo × cepa es estimable (48 gl).

Fórmula empleada: `y ~ C(gen) + C(strain) + C(year) + C(gen):C(strain)`. Tabla ANOVA en `tablas/anova.csv`. La columna `eta_cuad`
se calcula como `sum_sq` del término / `sum_sq` total (variación total explicada por el término).

| Término | gl | Sum sq | F | p | Significativo | η² |
|---|---|---|---|---|---|---|
| Genotipo | 16 | 12368.2 | 8.75 | 0.00 | Si | 0.298 |
| Cepa | 3 | 11461.1 | 43.22 | 0.00 | Si | 0.277 |
| Año | 2 | 454.5 | 2.57 | 0.08 | No | 0.011 |
| Genotipo x Cepa | 48 | 5308.2 | 1.25 | 0.16 | No | 0.128 |
| Residual | 134 | 11843.7 | nan | nan | No | 0.286 |

## 5. Comparaciones múltiples (Tukey HSD)

Comparaciones por pares con ajuste de Tukey (α = 0.05) sobre las medias de cepa y de genotipo.
Tablas completas en `tablas/posthoc_cepas.csv` y `tablas/posthoc_gen.csv`.

### 5.1 Pares entre cepas

| Cepa A | Cepa B | Diferencia | p-ajustada | Significativo |
|---|---|---|---|---|
| F329 | F348 | -0.64 | 0.99 | False |
| F329 | F39 | 17.94 | 0.00 | True |
| F329 | F436 | 4.05 | 0.34 | False |
| F348 | F39 | 18.58 | 0.00 | True |
| F348 | F436 | 4.69 | 0.22 | False |
| F39 | F436 | -13.89 | 0.00 | True |

### 5.2 Pares entre genotipos

| Genotipo A | Genotipo B | Diferencia | p-ajustada | Significativo |
|---|---|---|---|---|
| Arina | Nautica | 19.05 | 0.02 | True |
| Arina | SVP72003-4-2-4 | 3.69 | 1.00 | False |
| Arina | SVP72005-20-3-1 | 22.12 | 0.00 | True |
| Arina | SVP72017-17-5-10 | -1.32 | 1.00 | False |
| Arina | SVP73012-1-2-3 | 14.62 | 0.26 | False |
| Arina | SVP73016-2-4 | 12.15 | 0.59 | False |
| Arina | SVP73030-8-1-I | 9.70 | 0.89 | False |
| Arina | SVP75059-28 | 5.14 | 1.00 | False |
| Arina | SVP75059-32 | 20.89 | 0.01 | True |
| Arina | SVP75059-46 | 15.04 | 0.22 | False |
| Arina | SVP77076-1 | 0.29 | 1.00 | False |
| Arina | SVP77076-38 | 0.03 | 1.00 | False |
| Arina | SVP77076-4 | -0.25 | 1.00 | False |
| Arina | SVP77078-30 | 2.87 | 1.00 | False |
| Arina | SVP77079-15 | 5.08 | 1.00 | False |
| Arina | Saiga | 2.79 | 1.00 | False |
| Nautica | SVP72003-4-2-4 | -15.36 | 0.19 | False |
| Nautica | SVP72005-20-3-1 | 3.07 | 1.00 | False |
| Nautica | SVP72017-17-5-10 | -20.38 | 0.01 | True |
| Nautica | SVP73012-1-2-3 | -4.43 | 1.00 | False |
| Nautica | SVP73016-2-4 | -6.90 | 1.00 | False |
| Nautica | SVP73030-8-1-I | -9.35 | 0.92 | False |
| Nautica | SVP75059-28 | -13.91 | 0.34 | False |
| Nautica | SVP75059-32 | 1.84 | 1.00 | False |
| Nautica | SVP75059-46 | -4.01 | 1.00 | False |
| Nautica | SVP77076-1 | -18.76 | 0.03 | True |
| Nautica | SVP77076-38 | -19.02 | 0.02 | True |
| Nautica | SVP77076-4 | -19.30 | 0.02 | True |
| Nautica | SVP77078-30 | -16.18 | 0.13 | False |
| Nautica | SVP77079-15 | -13.97 | 0.34 | False |
| Nautica | Saiga | -16.26 | 0.12 | False |
| SVP72003-4-2-4 | SVP72005-20-3-1 | 18.43 | 0.04 | True |
| SVP72003-4-2-4 | SVP72017-17-5-10 | -5.02 | 1.00 | False |
| SVP72003-4-2-4 | SVP73012-1-2-3 | 10.93 | 0.76 | False |
| SVP72003-4-2-4 | SVP73016-2-4 | 8.46 | 0.96 | False |
| SVP72003-4-2-4 | SVP73030-8-1-I | 6.01 | 1.00 | False |
| SVP72003-4-2-4 | SVP75059-28 | 1.45 | 1.00 | False |
| SVP72003-4-2-4 | SVP75059-32 | 17.20 | 0.07 | False |
| SVP72003-4-2-4 | SVP75059-46 | 11.35 | 0.70 | False |
| SVP72003-4-2-4 | SVP77076-1 | -3.40 | 1.00 | False |
| SVP72003-4-2-4 | SVP77076-38 | -3.66 | 1.00 | False |
| SVP72003-4-2-4 | SVP77076-4 | -3.94 | 1.00 | False |
| SVP72003-4-2-4 | SVP77078-30 | -0.82 | 1.00 | False |
| SVP72003-4-2-4 | SVP77079-15 | 1.39 | 1.00 | False |
| SVP72003-4-2-4 | Saiga | -0.90 | 1.00 | False |
| SVP72005-20-3-1 | SVP72017-17-5-10 | -23.44 | 0.00 | True |
| SVP72005-20-3-1 | SVP73012-1-2-3 | -7.50 | 0.99 | False |
| SVP72005-20-3-1 | SVP73016-2-4 | -9.97 | 0.87 | False |
| SVP72005-20-3-1 | SVP73030-8-1-I | -12.42 | 0.55 | False |
| SVP72005-20-3-1 | SVP75059-28 | -16.98 | 0.08 | False |
| SVP72005-20-3-1 | SVP75059-32 | -1.23 | 1.00 | False |
| SVP72005-20-3-1 | SVP75059-46 | -7.08 | 0.99 | False |
| SVP72005-20-3-1 | SVP77076-1 | -21.82 | 0.00 | True |
| SVP72005-20-3-1 | SVP77076-38 | -22.08 | 0.00 | True |
| SVP72005-20-3-1 | SVP77076-4 | -22.37 | 0.00 | True |
| SVP72005-20-3-1 | SVP77078-30 | -19.25 | 0.02 | True |
| SVP72005-20-3-1 | SVP77079-15 | -17.03 | 0.08 | False |
| SVP72005-20-3-1 | Saiga | -19.32 | 0.02 | True |
| SVP72017-17-5-10 | SVP73012-1-2-3 | 15.94 | 0.14 | False |
| SVP72017-17-5-10 | SVP73016-2-4 | 13.47 | 0.40 | False |
| SVP72017-17-5-10 | SVP73030-8-1-I | 11.03 | 0.75 | False |
| SVP72017-17-5-10 | SVP75059-28 | 6.47 | 1.00 | False |
| SVP72017-17-5-10 | SVP75059-32 | 22.22 | 0.00 | True |
| SVP72017-17-5-10 | SVP75059-46 | 16.37 | 0.11 | False |
| SVP72017-17-5-10 | SVP77076-1 | 1.62 | 1.00 | False |
| SVP72017-17-5-10 | SVP77076-38 | 1.36 | 1.00 | False |
| SVP72017-17-5-10 | SVP77076-4 | 1.07 | 1.00 | False |
| SVP72017-17-5-10 | SVP77078-30 | 4.19 | 1.00 | False |
| SVP72017-17-5-10 | SVP77079-15 | 6.41 | 1.00 | False |
| SVP72017-17-5-10 | Saiga | 4.12 | 1.00 | False |
| SVP73012-1-2-3 | SVP73016-2-4 | -2.47 | 1.00 | False |
| SVP73012-1-2-3 | SVP73030-8-1-I | -4.92 | 1.00 | False |
| SVP73012-1-2-3 | SVP75059-28 | -9.47 | 0.91 | False |
| SVP73012-1-2-3 | SVP75059-32 | 6.28 | 1.00 | False |
| SVP73012-1-2-3 | SVP75059-46 | 0.42 | 1.00 | False |
| SVP73012-1-2-3 | SVP77076-1 | -14.32 | 0.29 | False |
| SVP73012-1-2-3 | SVP77076-38 | -14.58 | 0.26 | False |
| SVP73012-1-2-3 | SVP77076-4 | -14.87 | 0.23 | False |
| SVP73012-1-2-3 | SVP77078-30 | -11.75 | 0.65 | False |
| SVP73012-1-2-3 | SVP77079-15 | -9.53 | 0.90 | False |
| SVP73012-1-2-3 | Saiga | -11.82 | 0.64 | False |
| SVP73016-2-4 | SVP73030-8-1-I | -2.45 | 1.00 | False |
| SVP73016-2-4 | SVP75059-28 | -7.01 | 0.99 | False |
| SVP73016-2-4 | SVP75059-32 | 8.74 | 0.95 | False |
| SVP73016-2-4 | SVP75059-46 | 2.89 | 1.00 | False |
| SVP73016-2-4 | SVP77076-1 | -11.86 | 0.63 | False |
| SVP73016-2-4 | SVP77076-38 | -12.12 | 0.59 | False |
| SVP73016-2-4 | SVP77076-4 | -12.40 | 0.55 | False |
| SVP73016-2-4 | SVP77078-30 | -9.28 | 0.92 | False |
| SVP73016-2-4 | SVP77079-15 | -7.07 | 0.99 | False |
| SVP73016-2-4 | Saiga | -9.36 | 0.92 | False |
| SVP73030-8-1-I | SVP75059-28 | -4.56 | 1.00 | False |
| SVP73030-8-1-I | SVP75059-32 | 11.19 | 0.73 | False |
| SVP73030-8-1-I | SVP75059-46 | 5.34 | 1.00 | False |
| SVP73030-8-1-I | SVP77076-1 | -9.41 | 0.91 | False |
| SVP73030-8-1-I | SVP77076-38 | -9.67 | 0.89 | False |
| SVP73030-8-1-I | SVP77076-4 | -9.95 | 0.87 | False |
| SVP73030-8-1-I | SVP77078-30 | -6.83 | 1.00 | False |
| SVP73030-8-1-I | SVP77079-15 | -4.62 | 1.00 | False |
| SVP73030-8-1-I | Saiga | -6.91 | 1.00 | False |
| SVP75059-28 | SVP75059-32 | 15.75 | 0.16 | False |
| SVP75059-28 | SVP75059-46 | 9.90 | 0.87 | False |
| SVP75059-28 | SVP77076-1 | -4.85 | 1.00 | False |
| SVP75059-28 | SVP77076-38 | -5.11 | 1.00 | False |
| SVP75059-28 | SVP77076-4 | -5.39 | 1.00 | False |
| SVP75059-28 | SVP77078-30 | -2.27 | 1.00 | False |
| SVP75059-28 | SVP77079-15 | -0.06 | 1.00 | False |
| SVP75059-28 | Saiga | -2.35 | 1.00 | False |
| SVP75059-32 | SVP75059-46 | -5.85 | 1.00 | False |
| SVP75059-32 | SVP77076-1 | -20.60 | 0.01 | True |
| SVP75059-32 | SVP77076-38 | -20.86 | 0.01 | True |
| SVP75059-32 | SVP77076-4 | -21.14 | 0.01 | True |
| SVP75059-32 | SVP77078-30 | -18.02 | 0.04 | True |
| SVP75059-32 | SVP77079-15 | -15.81 | 0.15 | False |
| SVP75059-32 | Saiga | -18.10 | 0.04 | True |
| SVP75059-46 | SVP77076-1 | -14.75 | 0.25 | False |
| SVP75059-46 | SVP77076-38 | -15.01 | 0.22 | False |
| SVP75059-46 | SVP77076-4 | -15.29 | 0.19 | False |
| SVP75059-46 | SVP77078-30 | -12.18 | 0.59 | False |
| SVP75059-46 | SVP77079-15 | -9.96 | 0.87 | False |
| SVP75059-46 | Saiga | -12.25 | 0.58 | False |
| SVP77076-1 | SVP77076-38 | -0.26 | 1.00 | False |
| SVP77076-1 | SVP77076-4 | -0.54 | 1.00 | False |
| SVP77076-1 | SVP77078-30 | 2.58 | 1.00 | False |
| SVP77076-1 | SVP77079-15 | 4.79 | 1.00 | False |
| SVP77076-1 | Saiga | 2.50 | 1.00 | False |
| SVP77076-38 | SVP77076-4 | -0.28 | 1.00 | False |
| SVP77076-38 | SVP77078-30 | 2.83 | 1.00 | False |
| SVP77076-38 | SVP77079-15 | 5.05 | 1.00 | False |
| SVP77076-38 | Saiga | 2.76 | 1.00 | False |
| SVP77076-4 | SVP77078-30 | 3.12 | 1.00 | False |
| SVP77076-4 | SVP77079-15 | 5.33 | 1.00 | False |
| SVP77076-4 | Saiga | 3.04 | 1.00 | False |
| SVP77078-30 | SVP77079-15 | 2.22 | 1.00 | False |
| SVP77078-30 | Saiga | -0.07 | 1.00 | False |
| SVP77079-15 | Saiga | -2.29 | 1.00 | False |

## 6. Interpretación biológica

- **Agresividad relativa de las cepas:** la cepa F39 mostró la mayor severidad media
  (24.30), y la cepa F348 la menor
  (5.72). La diferencia entre cepas fue
  estadísticamente significativa (F = 43.22),
  lo que indica diferencias reales de agresividad entre aislados.
- **Susceptibilidad relativa de los genotipos:** el genotipo SVP72005-20-3-1 presentó la mayor
  severidad media y el genotipo SVP72017-17-5-10 la menor. La diferencia entre genotipos fue
  significativa (F = 8.75),
  reflejando un gradiente de susceptibilidad relativa en el material evaluado.
- **Interacción genotipo × cepa:** la interacción genotipo × cepa no resultó estadísticamente significativa al 5%, lo que sugiere que los efectos de cepa y genotipo actúan de forma aproximadamente aditiva. Cuando una interacción es significativa, el
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
