# Resumen del diseño RCBD (BDCA)

- **Factor:** Tratamiento (cuatro niveles: R, T0, T1, T2).- **Factor:** Bloque (nueve niveles: B1, B2, B3, B4, B5, B6, B7, B8, B9).- **Unidades experimentales:** 36 observaciones.- **Distribución por tratamiento:** {'T2': 9, 'R': 9, 'T0': 9, 'T1': 9}.- **Distribución por bloque:** {'B1': 4, 'B2': 4, 'B3': 4, 'B4': 4, 'B5': 4, 'B6': 4, 'B7': 4, 'B8': 4, 'B9': 4}.- **Balance:** Un cultivo por celda trt × bloque (estándar RCBD).- **Variable respuesta:** rendimiento (yield) del cultivo en mg.


## Auditoría de carga
---
Total de filas: 36
Filas únicas (sin duplicados): 36
Filas con valores NA: 0

Balance por celda:
- Total de celdas: 36
- Celdas balanceadas: 36
- Celdas desbalanceadas: 0

Conteo por celda (trt x block):
  - R x B1: 1
  - R x B2: 1
  - R x B3: 1
  - R x B4: 1
  - R x B5: 1
  - R x B6: 1
  - R x B7: 1
  - R x B8: 1
  - R x B9: 1
  - T0 x B1: 1
  - T0 x B2: 1
  - T0 x B3: 1
  - T0 x B4: 1
  - T0 x B5: 1
  - T0 x B6: 1
  - T0 x B7: 1
  - T0 x B8: 1
  - T0 x B9: 1
  - T1 x B1: 1
  - T1 x B2: 1
  - T1 x B3: 1
  - T1 x B4: 1
  - T1 x B5: 1
  - T1 x B6: 1
  - T1 x B7: 1
  - T1 x B8: 1
  - T1 x B9: 1
  - T2 x B1: 1
  - T2 x B2: 1
  - T2 x B3: 1
  - T2 x B4: 1
  - T2 x B5: 1
  - T2 x B6: 1
  - T2 x B7: 1
  - T2 x B8: 1
  - T2 x B9: 1

Conteo por tratamiento: {'T2': 9, 'R': 9, 'T0': 9, 'T1': 9}
Conteo por bloque: {'B1': 4, 'B2': 4, 'B3': 4, 'B4': 4, 'B5': 4, 'B6': 4, 'B7': 4, 'B8': 4, 'B9': 4}

Extremos de rendimiento:
- Mínimo: 4.38
- Máximo: 6.54
- Media: 5.801944444444444
- Std: 0.509469300600113

Número de cultivos con rendimiento negativo: 0

¿Auditoría aprobada? **Sí**

## Exploración descriptiva (EDA)
---

### Resumen descriptivo por tratamiento
| trt   |   n |   media |   desviacion_estandar |   error_estandar |   ic95_inferior |   ic95_superior |   minimo |   maximo |
|:------|----:|--------:|----------------------:|-----------------:|----------------:|----------------:|---------:|---------:|
| R     |   9 |   5.942 |                 0.465 |            0.155 |           5.585 |           6.3   |     5.06 |     6.54 |
| T0    |   9 |   5.31  |                 0.452 |            0.151 |           4.963 |           5.657 |     4.38 |     5.82 |
| T1    |   9 |   5.868 |                 0.46  |            0.153 |           5.514 |           6.222 |     5.04 |     6.45 |
| T2    |   9 |   6.088 |                 0.335 |            0.112 |           5.83  |           6.345 |     5.63 |     6.48 |

### Figuras exploratorias generadas
- `eda_boxplot_por_trt`: "eda_boxplot_por_trt.png"
- `eda_histogramas_por_trt`: "eda_histogramas_por_trt.png"
- `eda_qqplot_por_trt`: "eda_qqplot_por_trt.png"

## Verificación de supuestos
---

### Tabla de supuestos
| supuesto                             |   estadistico |    p_valor | cumple   |
|:-------------------------------------|--------------:|-----------:|:---------|
| Normalidad (Shapiro-Wilk)            |      0.981668 |   0.800165 | True     |
| Homocedasticidad (Levene)            |      0.102213 |   0.958159 | True     |
| Independencia serial (Durbin-Watson) |      2.07443  | nan        | True     |

### Decisión: PARAMETRICA
**Justificación:** Los residuos son normales (Shapiro-Wilk p=0.8002), homocedásticos (Levene p=0.9582) y muestran baja autocorrelación serial (Durbin-Watson = 2.0744). El ANOVA RCBD como tabla F es adecuado.

## ANOVA clásico de bloques RCBD
---

### Tabla ANOVA (tipo II)
| fuente      |   sum_sq |   df |        F |        PR(>F) |   eta2_parcial |
|:------------|---------:|-----:|---------:|--------------:|---------------:|
| tratamiento | 3.1295   |    3 |  28.7728 |   4.04874e-08 |         0.7824 |
| bloque      | 5.08494  |    8 |  17.5317 |   2.79121e-08 |         0.8539 |
| Residual    | 0.870128 |   24 | nan      | nan           |       nan      |

### Tamaños de efecto (eta2 parcial)
- **tratamiento**: 0.7824
- **bloque**: 0.8539

## Modelo mixto lineal (LMM) de bloques RCBD
---

### Efectos fijos (tratamiento)
| efecto       |   coeficiente |   error_estandar |       t |   p_valor |   ic95_inferior |   ic95_superior |
|:-------------|--------------:|-----------------:|--------:|----------:|----------------:|----------------:|
| Intercept    |        5.9422 |           0.1438 | 41.3241 |    0      |          5.6604 |          6.2241 |
| C(trt)[T.T0] |       -0.6322 |           0.0898 | -7.0435 |    0      |         -0.8082 |         -0.4563 |
| C(trt)[T.T1] |       -0.0744 |           0.0898 | -0.8294 |    0.4069 |         -0.2504 |          0.1015 |
| C(trt)[T.T2] |        0.1456 |           0.0898 |  1.6216 |    0.1049 |         -0.0304 |          0.3215 |
| Group Var    |        4.1329 |           2.5304 |  1.6333 |    0.1024 |         -0.8268 |          9.0925 |

### Varianzas e ICC
- Varianza de bloque: 0.1498
- Varianza residual: 0.0363
- ICC: 0.8052

### Limitación de aditividad
Dado un solo cultivo por celda trt × bloque, el término de interacción no puede ser estimado. Por lo tanto, la aditividad (efecto aditivo puro) es una suposición no testable; la inferencia se basa en el modelo de bloques RCBD sin interacción.

## Comparaciones múltiples post-hoc Tukey HSD
---

### Resultados Tukey HSD (referencia: R)
| par      |   diferencia_medias |   p_valor_ajustado |   ic95_inferior |   ic95_superior | significativo   | vs_referencia_R   |
|:---------|--------------------:|-------------------:|----------------:|----------------:|:----------------|:------------------|
| R vs T0  |          -0.632222  |         0.0194982  |     -1.18319    |      -0.0812508 | True            | True              |
| R vs T1  |          -0.0744444 |         0.982927   |     -0.625416   |       0.476527  | False           | True              |
| R vs T2  |           0.145556  |         0.890102   |     -0.405416   |       0.696527  | False           | True              |
| T0 vs T1 |           0.557778  |         0.0463564  |      0.00680639 |       1.10875   | True            | False             |
| T0 vs T2 |           0.777778  |         0.00305469 |      0.226806   |       1.32875   | True            | False             |
| T1 vs T2 |           0.22      |         0.702924   |     -0.330971   |       0.770971  | False           | False             |

### Comparaciones vs Referencia R
| par     |   diferencia_medias |   p_valor_ajustado |   ic95_inferior |   ic95_superior | significativo   | vs_referencia_R   |
|:--------|--------------------:|-------------------:|----------------:|----------------:|:----------------|:------------------|
| R vs T0 |          -0.632222  |          0.0194982 |       -1.18319  |      -0.0812508 | True            | True              |
| R vs T1 |          -0.0744444 |          0.982927  |       -0.625416 |       0.476527  | False           | True              |
| R vs T2 |           0.145556  |          0.890102  |       -0.405416 |       0.696527  | False           | True              |

## Documentación de variables derivadas
---

### EDA -> medias, IC95%
**Fuente:** `pipeline/bdca/eda.py -> resumen_descriptivo`
**Fórmula/Justificación:** Cada tratamiento: media de yield, IC95 via t de Student (df=n-1)

### ANOVA -> eta2 parcial
**Fuente:** `pipeline/bdca/modelos.py -> _calcular_eta2_parcial`
**Fórmula/Justificación:** Eta2 parcial = SS_efecto / (SS_efecto + SS_residual)

### LMM -> ICC
**Fuente:** `pipeline/bdca/modelos.py -> lmm_bloques`
**Fórmula/Justificación:** ICC = var_bloque / (var_bloque + var_residual) para estimar correlación intra-clase de bloque

### Post-hoc -> Tukey
**Fuente:** `pipeline/bdca/comparaciones.py -> posthoc_tukey`
**Fórmula/Justificación:** Tukey HSD por pares de tratamientos (statsmodels.pairwise_tukeyhsd) con columna 'vs_referencia_R'

### Supuestos -> tabla
**Fuente:** `pipeline/bdca/supuestos.py -> analisis_supuestos`
**Fórmula/Justificación:** Tabla de supuestos de normalidad (Shapiro-Wilk), homocedasticidad (Levene), independencia serial (Durbin-Watson)