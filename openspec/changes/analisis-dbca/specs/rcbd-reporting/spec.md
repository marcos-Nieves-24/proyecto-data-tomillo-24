# RCBD Reporting Specification

## Purpose

Generates the full DBCA deliverable at the DCA standard: an educational
notebook `bdca/analisis_bdca.ipynb`, complete results under
`bdca/resultados/`, and a real static page `pagina/bdca/index.html` with the
hub listing bdca as active.

## Requirements

### Requirement: Educational notebook for the bdca design

`generar_notebook_pipeline.py` MUST produce `bdca/analisis_bdca.ipynb` with a
single-response structure (yield): each major section combines markdown
explanation, statistical rationale, code, output and interpretation in
neutral professional Spanish, mirroring the DCA notebook standard without
reusing its four-block storytelling.

#### Scenario: Notebook generated and executable

- GIVEN `PIPELINE_DISENO=bdca`
- WHEN `generar_notebook_pipeline.py` runs
- THEN `bdca/analisis_bdca.ipynb` is created
- AND it executes end-to-end, with every section containing explanation and interpretation

### Requirement: Results artifacts under bdca/resultados/

The pipeline MUST write the DCA-standard results under `bdca/resultados/`:
tables (audit, assumptions, ANOVA, LMM variances/ICC, Tukey pairs, means and
CI) in `tablas/`, figures (boxplot, residual diagnostics, Tukey plot) in
`figuras/`, markdown and HTML report in `reportes/`, and a formatted Excel
summary in `excel/`.

#### Scenario: Full results produced

- GIVEN the bdca pipeline run completes
- WHEN output directories are inspected
- THEN each of `database`, `tablas`, `figuras`, `reportes` and `excel` contains its expected artifacts

#### Scenario: DCA results untouched

- GIVEN a bdca run completes
- WHEN `dca/resultados/` is compared before and after
- THEN no DCA artifact is created, modified or removed

### Requirement: Static page and hub activation

`generar_pagina.py` MUST render `pagina/bdca/index.html` with the bdca
results (replacing the "Próximamente" placeholder) using interactive Plotly
figures, and the hub `pagina/index.html` MUST mark the bdca card as active
while leaving dca active and factorial as future.

#### Scenario: bdca page renders with results

- GIVEN a completed bdca run
- WHEN `generar_pagina.py` runs for the bdca design
- THEN `pagina/bdca/index.html` contains results and figures instead of the placeholder

#### Scenario: Hub statuses updated

- GIVEN the hub is regenerated
- WHEN `pagina/index.html` is inspected
- THEN the dca and bdca cards are marked active ("Disponible")
- AND the factorial card remains "Próximamente"
