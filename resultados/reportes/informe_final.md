## Resumen ejecutivo

Este informe resume el análisis estadístico reproducible de la actividad antifúngica de extractos de tomillo (*Thymus vulgaris*) obtenidos por tres técnicas de extracción (maceración, Soxhlet y ultrasonido) frente a 31 aislados de *Fusarium* spp., todos ensayados a 5 mg/mL. Se evaluó el rendimiento de extracción, la inhibición del crecimiento micelial y la inhibición de la producción de conidias, se compararon las técnicas, se agruparon los aislados por su perfil de susceptibilidad y se generó un ranking de técnicas. La técnica mejor puntuada fue **Maceración** (score compuesto 0.667). El análisis multivariado identificó **2** grupos de aislados.

Todos los resultados numéricos se guardaron como tablas en `resultados/tablas/` y las figuras en `resultados/figuras/`.


## 1. Calidad de datos (auditoría)

Se auditó el dataset consolidado tidy (279 filas, 7 columnas) antes de cualquier transformación: completitud, duplicados, columnas constantes, tipos, valores inconsistentes y atípicos por rango intercuartílico. Ningún atípico se eliminó automáticamente (regla de integridad del proyecto); los valores se conservaron y su efecto se evaluó en el análisis.

| variable | tipo_dato | n_no_nulos | pct_faltantes | columna_constante | n_filas_duplicadas | n_valores_inconsistentes | n_atipicos_iqr |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Metodo de extraccion | str | 279 | 0.0 | No | 0 | 0 | 0 |
| Aislamiento | str | 279 | 0.0 | No | 0 | 0 | 0 |
| Replica | str | 279 | 0.0 | No | 0 | 0 | 0 |
| Crecimiento micelial (mm) | int64 | 279 | 0.0 | No | 0 | 0 | 0 |
| %INH micelial | float64 | 279 | 0.0 | No | 0 | 0 | 1 |
| Conidias (log10/ml) | float64 | 279 | 0.0 | No | 0 | 0 | 8 |
| %INH conidias | float64 | 279 | 0.0 | No | 0 | 0 | 9 |

Filas duplicadas exactas: **0**. Atípicos por IQR (1.5xIQR) flagueados, sin eliminar: {'%INH micelial': 1, 'Conidias (log10/ml)': 8, '%INH conidias': 9}.


## 2. Diseño experimental (inferido y caveats)

El diseño se infiere como un DCA factorial con dos factores fijos (técnica de extracción × aislado de Fusarium) y tres réplicas biológicas por combinación; la unidad experimental es la caja Petri. La concentración (5 mg/mL) es constante y NO es factor experimental. Caveat: cada porcentaje de inhibición se calculó contra un único control C4 compartido por las tres réplicas del aislado; por lo tanto, las réplicas de %INH no son totalmente independientes (pseudorreplicación del control). El análisis del crecimiento crudo (mm) no presenta este problema.

| atributo | valor |
| --- | --- |
| Tipo de diseño | DCA factorial método × aislado |
| Factores | método de extracción (fijo); aislado (fijo en ANOVA, aleatorio en LMM) |
| Niveles de método | 3 |
| Niveles de aislado | 31 |
| Número de tratamientos (método × aislado) | 93 |
| Réplicas biológicas por celda | 3 |
| Unidad experimental | Caja Petri |
| Concentración ensayada (mg/mL) | 5.0 |
| Variables de respuesta | Crecimiento micelial (mm); Inhibición micelial (%); Conidias (log10/mL); Inhibición de conidias (%) |
| Diseño balanceado | Sí |


### 2.1 Validación de la fórmula del %INH (integridad de datos)

Se reconstruyó el %INH a partir de los controles C4 y las mediciones del bioensayo para confirmar la consistencia de los datos de entrada. La validación es informativa y no reemplaza las respuestas reportadas usadas en la inferencia.

| variable | n_verificadas | max_diff_abs | n_discrepancias | estado | nota |
| --- | --- | --- | --- | --- | --- |
| porcentaje_inhibicion_micelial | 279 | 0.0 | 0 | ok | La fórmula (1 - C1/C4) x 100 coincide con el %INH reportado (tolerancia 1e-6). |
| porcentaje_inhibicion_conidias | 279 | 0.0 | 0 | ok | La fórmula (1 - C1/C4) x 100 coincide con el %INH reportado (tolerancia 1e-6). |

La fórmula reconstruida `%INH = (1 - C1/C4) × 100` coincide con el %INH reportado por el laboratorio (diferencias máximas del orden de 1e-10 y ninguna discrepancia por encima de la tolerancia 1e-6). Esto confirma la consistencia interna de los datos: los controles C4 se incorporaron al dataset maestro (columnas `control_crecimiento_mm` y `control_conidias_log10`) y el %INH usado en la inferencia es el reportado por el investigador. Para las conidias, la fórmula se validó sobre la escala log10 directamente (el laboratorio reportó la reducción en log10), tal como se verificó en los datos.



## 3. Supuestos de los modelos

Se verificaron normalidad de residuos (Shapiro-Wilk), homocedasticidad (Levene y Bartlett, por método) e independencia (Durbin-Watson) sobre el modelo OLS factorial método × aislado.

- **Crecimiento micelial (mm)**: Shapiro-Wilk p=0.0328, Levene p=0.0628, Bartlett p=0.0819, Durbin-Watson=2.71.
- **Inhibición micelial (%)**: Shapiro-Wilk p=0.0788, Levene p=0.0982, Bartlett p=0.0760, Durbin-Watson=2.70.
- **Conidias (log10/mL)**: Shapiro-Wilk p=0.0000, Levene p=0.0000, Bartlett p=0.0000, Durbin-Watson=2.46.
- **Inhibición de conidias (%)**: Shapiro-Wilk p=0.0000, Levene p=0.0001, Bartlett p=0.0000, Durbin-Watson=2.41.


## 4. Análisis seleccionados y por qué

### 4.1 Rendimiento de extracción (ANOVA de una vía)

**Modelo seleccionado:** anova_ols. eta2=0.9443; omega2=0.9172; Kruskal-Wallis p=0.0273.

| fuente | sum_sq | df | F | PR(>F) | eta2_parcial | omega2_parcial |
| --- | --- | --- | --- | --- | --- | --- |
| metodo_extraccion | 1698.3726 | 2.0 | 50.8412 | 0.0002 | 0.9443 | 0.9172 |
| Residual | 100.2163 | 6.0 |  |  |  |  |

Justificación: El ANOVA de una via es adecuado si los residuos son normales y homocedasticos. Ambos supuestos se cumplen; se usa la tabla F como inferencia principal. Kruskal-Wallis confirma diferencias entre metodos.

### 4.2 Crecimiento micelial (mm) (unidad: mm)

**Modelo seleccionado:** factorial_no_parametrico. Kruskal-Wallis (método): H=132.866, p=0.0000. Scheirer-Ray-Hare: método H=132.866 (p=0.0000); interacción H=30.881 (p=0.9993).

| fuente | sum_sq | df | F | PR(>F) | eta2_parcial | omega2_parcial |
| --- | --- | --- | --- | --- | --- | --- |
| metodo_extraccion | 16593.4265 | 2.0 | 244.0724 | 0.0 | 0.7241 | 0.7201 |
| aislamiento | 5256.7885 | 30.0 | 5.1548 | 0.0 | 0.454 | 0.3648 |
| metodo_extraccion:aislamiento | 4233.0179 | 60.0 | 2.0754 | 0.0001 | 0.401 | 0.2071 |
| Residual | 6322.6667 | 186.0 |  |  |  |  |

Justificación: Los supuestos del ANOVA factorial no se cumplen (Shapiro-Wilk p=0.0328; Levene p=0.0628). Se utiliza la via no parametrica (Kruskal-Wallis por metodo y Scheirer-Ray-Hare para la interaccion) como inferencia principal, conservando la tabla ANOVA y los tamanos de efecto como referencia descriptiva. Para %INH micelial, el efecto techo (muchos valores = 100) explica la violacion de normalidad.

### 4.3 Inhibición micelial (%) (unidad: %)

**Modelo seleccionado:** factorial_ols.

| fuente | sum_sq | df | F | PR(>F) | eta2_parcial | omega2_parcial |
| --- | --- | --- | --- | --- | --- | --- |
| metodo_extraccion | 53840.3359 | 2.0 | 254.6034 | 0.0 | 0.7325 | 0.7285 |
| aislamiento | 21504.0667 | 30.0 | 6.7793 | 0.0 | 0.5223 | 0.4441 |
| metodo_extraccion:aislamiento | 13887.0985 | 60.0 | 2.189 | 0.0 | 0.4139 | 0.2241 |
| Residual | 19666.4716 | 186.0 |  |  |  |  |

Justificación: Los supuestos de normalidad y homocedasticidad se cumplen (Shapiro-Wilk y Levene p>0.05); el ANOVA factorial es la via adecuada. Nota: las replicas de %INH comparten el control C4 (pseudorreplicacion), lo que puede inflar levemente la precision.

### 4.4 Conidias (log10/mL) (unidad: log10(conidias/mL))

**Modelo seleccionado:** factorial_no_parametrico. Kruskal-Wallis (método): H=151.225, p=0.0000. Scheirer-Ray-Hare: método H=151.225 (p=0.0000); interacción H=44.272 (p=0.9360).

| fuente | sum_sq | df | F | PR(>F) | eta2_parcial | omega2_parcial |
| --- | --- | --- | --- | --- | --- | --- |
| metodo_extraccion | 238.945 | 2.0 | 263.5756 | 0.0 | 0.7392 | 0.7354 |
| aislamiento | 98.9316 | 30.0 | 7.2753 | 0.0 | 0.5399 | 0.4645 |
| metodo_extraccion:aislamiento | 160.537 | 60.0 | 5.9028 | 0.0 | 0.6557 | 0.5436 |
| Residual | 84.3093 | 186.0 |  |  |  |  |

Justificación: Los supuestos del ANOVA factorial no se cumplen (Shapiro-Wilk p=0.0000; Levene p=0.0000). Se utiliza la via no parametrica (Kruskal-Wallis por metodo y Scheirer-Ray-Hare para la interaccion) como inferencia principal, conservando la tabla ANOVA y los tamanos de efecto como referencia descriptiva. Para %INH micelial, el efecto techo (muchos valores = 100) explica la violacion de normalidad.

### 4.5 Inhibición de conidias (%) (unidad: %)

**Modelo seleccionado:** factorial_no_parametrico. Kruskal-Wallis (método): H=152.936, p=0.0000. Scheirer-Ray-Hare: método H=152.936 (p=0.0000); interacción H=45.794 (p=0.9121).

| fuente | sum_sq | df | F | PR(>F) | eta2_parcial | omega2_parcial |
| --- | --- | --- | --- | --- | --- | --- |
| metodo_extraccion | 46545.0725 | 2.0 | 251.2568 | 0.0 | 0.7299 | 0.7259 |
| aislamiento | 12279.5895 | 30.0 | 4.4191 | 0.0 | 0.4161 | 0.321 |
| metodo_extraccion:aislamiento | 37082.6752 | 60.0 | 6.6726 | 0.0 | 0.6828 | 0.5795 |
| Residual | 17228.1588 | 186.0 |  |  |  |  |

Justificación: Los supuestos del ANOVA factorial no se cumplen (Shapiro-Wilk p=0.0000; Levene p=0.0001). Se utiliza la via no parametrica (Kruskal-Wallis por metodo y Scheirer-Ray-Hare para la interaccion) como inferencia principal, conservando la tabla ANOVA y los tamanos de efecto como referencia descriptiva. Para %INH micelial, el efecto techo (muchos valores = 100) explica la violacion de normalidad.

### 4.6 Sensibilidad LMM - Crecimiento micelial (mm)

Modelo mixto `crecimiento_micelial_mm ~ método + (1|aislamiento)`. ICC (aislado)=0.2552; p del método (LRT)=0.0000.

| efecto | coeficiente | error_estandar | t | p_valor | ic95_inferior | ic95_superior |
| --- | --- | --- | --- | --- | --- | --- |
| Intercept | 7.4839 | 0.9673 | 7.7369 | 0.0 | 5.588 | 9.3798 |
| C(metodo_extraccion)[T.soxhlet] | 14.7204 | 0.9606 | 15.324 | 0.0 | 12.8376 | 16.6032 |
| C(metodo_extraccion)[T.ultrasonido] | 17.6129 | 0.9606 | 18.3351 | 0.0 | 15.7301 | 19.4957 |
| Group Var | 0.3426 | 0.1241 | 2.761 | 0.0058 | 0.0994 | 0.5859 |

### 4.7 Sensibilidad LMM - Inhibición micelial (%)

Modelo mixto `porcentaje_inhibicion_micelial ~ método + (1|aislamiento)`. ICC (aislado)=0.3210; p del método (LRT)=0.0000.

| efecto | coeficiente | error_estandar | t | p_valor | ic95_inferior | ic95_superior |
| --- | --- | --- | --- | --- | --- | --- |
| Intercept | 86.349 | 1.8833 | 45.8487 | 0.0 | 82.6577 | 90.0404 |
| C(metodo_extraccion)[T.soxhlet] | -26.7564 | 1.7127 | -15.6226 | 0.0 | -30.1132 | -23.3995 |
| C(metodo_extraccion)[T.ultrasonido] | -31.5844 | 1.7127 | -18.4415 | 0.0 | -34.9412 | -28.2275 |
| Group Var | 0.4728 | 0.1597 | 2.9606 | 0.0031 | 0.1598 | 0.7858 |

### 4.8 Sensibilidad LMM - Conidias (log10/mL)

Modelo mixto `conidias_log10_ml ~ método + (1|aislamiento)`. ICC (aislado)=0.2045; p del método (LRT)=0.0000.

| efecto | coeficiente | error_estandar | t | p_valor | ic95_inferior | ic95_superior |
| --- | --- | --- | --- | --- | --- | --- |
| Intercept | 5.2683 | 0.1377 | 38.2658 | 0.0 | 4.9984 | 5.5381 |
| C(metodo_extraccion)[T.soxhlet] | 1.7724 | 0.1463 | 12.1144 | 0.0 | 1.4856 | 2.0591 |
| C(metodo_extraccion)[T.ultrasonido] | 2.1101 | 0.1463 | 14.4229 | 0.0 | 1.8234 | 2.3969 |
| Group Var | 0.257 | 0.1007 | 2.5529 | 0.0107 | 0.0597 | 0.4544 |

### 4.9 Sensibilidad LMM - Inhibición de conidias (%)

Modelo mixto `porcentaje_inhibicion_conidias ~ método + (1|aislamiento)`. ICC (aislado)=0.0867; p del método (LRT)=0.0000.

| efecto | coeficiente | error_estandar | t | p_valor | ic95_inferior | ic95_superior |
| --- | --- | --- | --- | --- | --- | --- |
| Intercept | 29.2394 | 1.7463 | 16.7432 | 0.0 | 25.8166 | 32.6622 |
| C(metodo_extraccion)[T.soxhlet] | -24.8071 | 2.179 | -11.3848 | 0.0 | -29.0779 | -20.5363 |
| C(metodo_extraccion)[T.ultrasonido] | -29.4089 | 2.179 | -13.4968 | 0.0 | -33.6797 | -25.1382 |
| Group Var | 0.0949 | 0.0563 | 1.6843 | 0.0921 | -0.0155 | 0.2053 |



## 5. Comparaciones múltiples

### 5.1 Crecimiento micelial (mm) (Dunn (FDR))

Letras CLD (métodos que comparten letra no difieren, p>=0.05): **Maceración=c; Soxhlet=b; Ultrasonido=a**.

| par | estadistico_z | p_valor | p_valor_ajustado | letras |
| --- | --- | --- | --- | --- |
| maceracion vs soxhlet | -8.8304 | 0.0 | 0.0 | sin letra comun |
| maceracion vs ultrasonido | -10.8314 | 0.0 | 0.0 | sin letra comun |
| soxhlet vs ultrasonido | -2.001 | 0.0454 | 0.0454 | sin letra comun |

### 5.2 Inhibición micelial (%) (Tukey HSD)

Letras CLD (métodos que comparten letra no difieren, p>=0.05): **Maceración=a; Soxhlet=b; Ultrasonido=b**.

| par | diferencia_medias | p_valor_ajustado | ic95_inferior | ic95_superior | significativo | letras |
| --- | --- | --- | --- | --- | --- | --- |
| maceracion vs soxhlet | -26.7564 | 0.0 | -31.6371 | -21.8757 | Sí | sin letra comun |
| maceracion vs ultrasonido | -31.5844 | 0.0 | -36.4651 | -26.7037 | Sí | sin letra comun |
| soxhlet vs ultrasonido | -4.828 | 0.0533 | -9.7087 | 0.0527 | No | b |

### 5.3 Conidias (log10/mL) (Dunn (FDR))

Letras CLD (métodos que comparten letra no difieren, p>=0.05): **Maceración=c; Soxhlet=b; Ultrasonido=a**.

| par | estadistico_z | p_valor | p_valor_ajustado | letras |
| --- | --- | --- | --- | --- |
| maceracion vs soxhlet | -8.7954 | 0.0 | 0.0 | sin letra comun |
| maceracion vs ultrasonido | -11.8408 | 0.0 | 0.0 | sin letra comun |
| soxhlet vs ultrasonido | -3.0454 | 0.0023 | 0.0023 | sin letra comun |

### 5.4 Inhibición de conidias (%) (Dunn (FDR))

Letras CLD (métodos que comparten letra no difieren, p>=0.05): **Maceración=a; Soxhlet=b; Ultrasonido=c**.

| par | estadistico_z | p_valor | p_valor_ajustado | letras |
| --- | --- | --- | --- | --- |
| maceracion vs soxhlet | 9.0219 | 0.0 | 0.0 | sin letra comun |
| maceracion vs ultrasonido | 11.8359 | 0.0 | 0.0 | sin letra comun |
| soxhlet vs ultrasonido | 2.814 | 0.0049 | 0.0049 | sin letra comun |



## 6. Análisis multivariado (susceptibilidad)


- Varianza explicada por PC1 y PC2: 42.2% y 29.2%.
- Coeficiente cofenético (Ward): 0.7356.
- Número óptimo de clusters KMeans (silhouette): 2.
- Categorías biológicas: terciles del score compuesto de susceptibilidad (promedio z de inhibición micelial y de conidias) etiquetadas como Alta / Moderada / Baja susceptibilidad relativa. No se utiliza el término 'resistente' por no existir un criterio validado.

### 6.1 Categorías y clusters por aislado

| aislamiento | score_susceptibilidad | categoria_susceptibilidad | cluster_kmeans | inhib_micelial_maceracion | inhib_micelial_soxhlet | inhib_micelial_ultrasonido | inhib_conidias_maceracion | inhib_conidias_soxhlet | inhib_conidias_ultrasonido |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H5N | 0.9236 | Alta susceptibilidad relativa | 1 | 100.0 | 62.2222 | 62.2222 | 28.4427 | 21.1068 | 11.0253 |
| HC3 | 0.8174 | Alta susceptibilidad relativa | 1 | 100.0 | 65.6863 | 60.2941 | 79.0043 | -3.4199 | 9.2208 |
| H6B | 0.7961 | Alta susceptibilidad relativa | 1 | 92.8205 | 76.4103 | 63.0769 | 24.6851 | 12.5105 | 7.6406 |
| H8N | 0.6564 | Alta susceptibilidad relativa | 1 | 92.7928 | 74.7748 | 61.2613 | 26.771 | 7.9965 | 7.4316 |
| HC26 | 0.5603 | Alta susceptibilidad relativa | 1 | 100.0 | 62.5731 | 62.5731 | 34.4385 | 11.7595 | 0.3095 |
| H9N | 0.5221 | Alta susceptibilidad relativa | 1 | 100.0 | 67.619 | 59.5238 | 25.7596 | 7.926 | 4.3593 |
| H4B | 0.5054 | Alta susceptibilidad relativa | 1 | 93.1373 | 74.5098 | 60.7843 | 24.7595 | 8.532 | 1.1292 |
| FU2 (UCMU21) | 0.4296 | Alta susceptibilidad relativa | 1 | 100.0 | 69.7917 | 43.2292 | 78.5132 | 2.2062 | -2.494 |
| HC28 | 0.4119 | Alta susceptibilidad relativa | 1 | 100.0 | 61.3333 | 57.3333 | 36.2069 | 16.0035 | -5.924 |
| HC9 | 0.3591 | Alta susceptibilidad relativa | 1 | 100.0 | 80.7292 | 60.4167 | 10.7884 | 0.3688 | -0.6455 |
| FUSARIUM JULIAN H20 | 0.358 | Moderada susceptibilidad relativa | 1 | 78.4314 | 70.5882 | 69.1176 | 16.715 | 2.4155 | 7.343 |
| FUSARIUM MARCE 1.2 | 0.3421 | Moderada susceptibilidad relativa | 1 | 91.1111 | 63.3333 | 50.0 | 20.5227 | 12.7249 | 9.6401 |
| H4G | 0.2603 | Moderada susceptibilidad relativa | 1 | 91.358 | 66.6667 | 56.1728 | 20.1044 | 9.8782 | 0.9574 |
| H4N | 0.2302 | Moderada susceptibilidad relativa | 1 | 82.7586 | 65.5172 | 63.2184 | 24.1004 | 6.0606 | 1.4205 |
| H6N | 0.176 | Moderada susceptibilidad relativa | 1 | 80.0 | 52.0 | 60.0 | 28.4288 | 8.2899 | 9.3316 |
| H3N | 0.1713 | Moderada susceptibilidad relativa | 1 | 100.0 | 58.6667 | 57.3333 | 27.154 | 3.0461 | 0.6963 |
| FU1 | 0.1114 | Moderada susceptibilidad relativa | 1 | 100.0 | 56.25 | 50.0 | 23.9879 | 6.3738 | 3.919 |
| H2N | 0.1053 | Moderada susceptibilidad relativa | 1 | 69.1358 | 60.4938 | 67.9012 | 27.0237 | 7.7626 | 0.3321 |
| HC23 | -0.209 | Moderada susceptibilidad relativa | 1 | 81.6092 | 58.6207 | 55.1724 | 11.1111 | 3.3656 | 1.1526 |
| HC15 | -0.2302 | Moderada susceptibilidad relativa | 1 | 100.0 | 53.3333 | 53.3333 | 20.2709 | -0.0 | -5.024 |
| HC10 | -0.2333 | Baja susceptibilidad relativa | 0 | 75.5556 | 51.1111 | 55.5556 | 25.2841 | 13.7311 | -6.6288 |
| HC17 | -0.2495 | Baja susceptibilidad relativa | 1 | 75.1111 | 59.1111 | 56.0 | 17.9253 | -0.9115 | 3.342 |
| H11G | -0.3142 | Baja susceptibilidad relativa | 0 | 100.0 | 55.1282 | 37.1795 | 78.425 | -9.8166 | -9.493 |
| H1N | -0.4654 | Baja susceptibilidad relativa | 0 | 62.6667 | 52.0 | 52.0 | 27.1224 | 3.0384 | 1.6979 |
| HC16 | -0.5788 | Baja susceptibilidad relativa | 0 | 80.4598 | 58.6207 | 57.4713 | 16.7883 | 4.9635 | -20.8759 |
| HC5 | -0.6451 | Baja susceptibilidad relativa | 0 | 65.2174 | 55.0725 | 63.7681 | 16.4646 | 2.8283 | -15.9091 |
| HC27 | -0.6724 | Baja susceptibilidad relativa | 0 | 70.6667 | 48.0 | 41.3333 | 19.5256 | 3.2613 | 2.2448 |
| H8G | -0.8173 | Baja susceptibilidad relativa | 0 | 100.0 | 47.7778 | 38.8889 | 78.5982 | -23.9165 | -15.4093 |
| HC20 | -0.8587 | Baja susceptibilidad relativa | 0 | 62.3188 | 44.9275 | 50.7246 | 2.6042 | 0.9549 | 1.3455 |
| HC6 | -0.955 | Baja susceptibilidad relativa | 0 | 67.9012 | 38.2716 | 37.037 | 20.5473 | 2.4675 | 1.8843 |
| HC19 | -1.5075 | Baja susceptibilidad relativa | 0 | 63.7681 | 36.2319 | 34.7826 | 14.3478 | -4.1063 | -9.2754 |

### 6.2 Cruce cluster x categoría

| Alta susceptibilidad relativa | Baja susceptibilidad relativa | Moderada susceptibilidad relativa |
| --- | --- | --- |
| 0 | 10 | 0 |
| 10 | 1 | 10 |


## 7. Ranking de técnicas

El score compuesto promedia las métricas normalizadas (min-max 0-1) de rendimiento, inhibición micelial e inhibición de conidias; mayores valores indican mejor desempeño global.

| ranking | metodo_extraccion | rendimiento_medio_pct | inhib_micelial_medio_pct | inhib_conidias_medio_pct | rendimiento_norm | inhib_micelial_norm | inhib_conidias_norm | score_compuesto |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | maceracion | 12.068 | 86.349 | 29.239 | 0.0 | 1.0 | 1.0 | 0.667 |
| 2 | soxhlet | 43.38 | 59.593 | 4.432 | 1.0 | 0.153 | 0.156 | 0.436 |
| 3 | ultrasonido | 17.053 | 54.765 | -0.17 | 0.159 | 0.0 | 0.0 | 0.053 |


## 8. Interpretación biológica


Los tres extractos (5 mg/mL) mostraron actividad antifúngica variable. La existencia de un efecto techo en la inhibición micelial (numerosos valores de 100 %) indica que, a esta concentración, muchos aislados fueron completamente inhibidos; por ello la inferencia principal sobre esta variable se apoyó en la vía seleccionada por los supuestos y en los modelos mixtos como análisis de sensibilidad.

La inhibición de conidias se expresa en escala log10; valores negativos indican que el extracto indujo una mayor esporulación relativa en esa celda. Esta característica debe considerarse al interpretar diferencias entre métodos.

La variabilidad entre aislados (ICC) fue evaluada con modelos mixtos; una ICC alta indica que el aislado explica una proporción importante de la variación total, lo que justifica la consideración de la susceptibilidad diferencial (análisis multivariado).

De acuerdo con el score compuesto, el ordenamiento de las técnicas fue: Maceración (0.667), Soxhlet (0.436), Ultrasonido (0.053).



## 9. Conclusiones


1. **Rendimiento**: el método de extracción afectó significativamente el rendimiento; Soxhlet mostró el mayor valor medio (ver sección 4.1).
2. **Inhibición micelial**: a 5 mg/mL se observó un efecto fuerte; las diferencias entre métodos y aislados fueron evaluadas según la vía seleccionada por los supuestos y con modelos mixtos de sensibilidad.
3. **Conidias**: la variable es continua en escala log10; se modeló con un modelo lineal sobre log10, sin aplicar Poisson/NB (los conteos crudos no están disponibles).
4. **Susceptibilidad**: los aislados se agruparon en categorías de susceptibilidad relativa (Alta / Moderada / Baja) a partir del score compuesto y de la clusterización.
5. **Ranking**:  La técnica Maceración encabezó el score compuesto, impulsada por su actividad antifúngica.


## 10. Limitaciones

- **Control compartido (pseudorreplicación)**: cada %INH se calculó contra un único control C4 compartido por las tres réplicas del aislado. Las réplicas de %INH no son totalmente independientes; el análisis del crecimiento crudo en mm no presenta este problema y los modelos mixtos con aislado aleatorio mitigan parcialmente la dependencia.
Los controles C4 ahora están explícitos en el dataset maestro (columnas `control_crecimiento_mm` y `control_conidias_log10`, una por aislado y compartidos por las 3 réplicas); su uso en futuras iteraciones permitiría recalcular %INH con otros controles o modelar directamente el crecimiento frente al control.
- **Escala log10 del %INH de conidias**: el laboratorio reportó la reducción de conidias en escala log10; no equivale a la reducción porcentual de conteos crudos.
- **Sin dosis-respuesta**: solo se ensayaron 5 mg/mL; no es posible estimar EC50/EC90 ni extrapolar a otras concentraciones.
- **Conidias como log10 continuas**: los conteos crudos no están disponibles; no aplica la rama Poisson/NB del pipeline (documentada para datasets futuros).
- **Tamaño muestral**: 3 réplicas biológicas por celda limitan la potencia; los clusters del análisis multivariado deben interpretarse con cautela por el tamaño de muestra.
