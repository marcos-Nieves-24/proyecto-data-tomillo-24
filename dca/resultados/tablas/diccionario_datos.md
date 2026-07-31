# Diccionario de datos - Master dataset (pipeline reproducible)

Fuente: dca/resultados/database/consolidado_tidy.xlsx (hoja 'Consolidado') y
dca/resultados/database/rendimiento_extraccion.csv.

## Hoja Bioensayo

| Variable | Unidad | Descripción |
|----------|--------|-------------|
| metodo_extraccion | - | Técnica de extracción (maceración, soxhlet, ultrasonido). Factor fijo. |
| aislamiento | - | Código del aislado de Fusarium spp. (31 aislados). |
| replica | - | Réplica biológica (1-3). |
| crecimiento_micelial_mm | mm | Crecimiento micelial a 5 mg/mL. |
| porcentaje_inhibicion_micelial | % | Inhibición calculada contra el control C4 del propio aislado. |
| conidias_log10_ml | log10(conidias/mL) | Concentración de conidias en escala log10 (transformada por el laboratorio). |
| porcentaje_inhibicion_conidias | % | Reducción de conidias en escala log10; puede ser negativa (mayor esporulación). |
| control_crecimiento_mm | mm | Control C4 del aislado: crecimiento micelial sin extracto (extraído del Excel crudo del laboratorio). Compartido por las 3 réplicas del aislado. |
| control_conidias_log10 | log10(conidias/mL) | Control C4 del aislado: conidias sin extracto en escala log10 (Excel crudo del laboratorio). Compartido por las 3 réplicas del aislado. |

## Hoja Rendimiento

| Variable | Unidad | Descripción |
|----------|--------|-------------|
| metodo_extraccion | - | Técnica de extracción. |
| peso_material_seco_g | g | Peso de material vegetal seco. |
| peso_extracto_obtenido_g | g | Peso de extracto obtenido. |
| rendimiento_pct | % | Rendimiento = (peso_extracto / peso_material_seco) × 100. |
| replica_biologica | - | Réplica biológica (1-3). |

## Variables derivadas y transformaciones

- Concentración única ensayada: 5 mg/mL (no es factor experimental).
- %INH micelial: 100 indica inhibición completa (efecto techo).
- %INH conidias: calculado por el laboratorio sobre la escala log10.
- Cada %INH se calculó contra un control C4 compartido por las 3 réplicas
  del aislado (pseudorreplicación del control; ver limitaciones).
- Las columnas control_* se extraen del Excel crudo del laboratorio
  (datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx, hojas por técnica) y se propagan
  a las 3 réplicas de cada aislado.
