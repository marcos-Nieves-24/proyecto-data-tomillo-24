# Pipeline Design Registry Specification

## Purpose

The pipeline is hardcoded to the DCA design (`DISENO = "dca"` in
`pipeline/config.py`, with output directories created at import time). This
capability introduces a design registry (`DISENOS`) and resolves the active
experimental design from the `PIPELINE_DISENO` environment variable (default
`dca`) before any directory is created, enabling the `bdca` design without
altering DCA behavior when the variable is unset.

## Requirements

### Requirement: Design registry and active-design resolution

The pipeline MUST define a `DISENOS` registry in `pipeline/config.py`
enumerating the supported experimental designs (`dca`, `bdca`), each with its
own input and output paths. The active design MUST be resolved from the
`PIPELINE_DISENO` environment variable and MUST default to `dca` when the
variable is unset or holds an unrecognized value.

#### Scenario: Default design when variable is unset

- GIVEN the environment variable `PIPELINE_DISENO` is not set
- WHEN `pipeline` is imported
- THEN the active design resolves to `dca`
- AND all paths point to the existing `dca/` and `datos_crudos/dca/` directories

#### Scenario: Explicit selection of the bdca design

- GIVEN `PIPELINE_DISENO=bdca` is set in the environment
- WHEN `pipeline` is imported
- THEN the active design resolves to `bdca`
- AND all paths point to `bdca/` and `datos_crudos/bdca/`

#### Scenario: Unrecognized design value

- GIVEN `PIPELINE_DISENO` holds an unrecognized value (e.g. `foo`)
- WHEN `pipeline` is imported
- THEN the pipeline falls back to the `dca` design without raising an error
- AND a warning describing the fallback is emitted

### Requirement: Design resolution precedes directory creation

The pipeline MUST resolve the active design BEFORE creating output
directories, so that directories are created under the resolved design's
paths and no DCA directory is touched when another design is active.

#### Scenario: Directories created for the active design

- GIVEN the active design resolves to `bdca`
- WHEN `pipeline` is imported
- THEN `bdca/resultados/{database,tablas,figuras,reportes,excel}` are created

#### Scenario: DCA directories untouched by a bdca run

- GIVEN the active design is `bdca`
- WHEN the pipeline runs to completion
- THEN no file under `dca/` or `datos_crudos/dca/` is created or modified

### Requirement: DCA behavior unchanged by default

With `PIPELINE_DISENO` unset (or set to `dca`), the pipeline MUST produce the
same DCA formulas, narratives, notebooks and artifacts as before the registry
change.

#### Scenario: DCA regeneration regression guard

- GIVEN `PIPELINE_DISENO` is unset
- WHEN `generar_notebook_pipeline.py` and `generar_pagina.py` run
- THEN the regenerated `dca/analisis_dca.ipynb` and `pagina/dca/index.html` match the previous versions
