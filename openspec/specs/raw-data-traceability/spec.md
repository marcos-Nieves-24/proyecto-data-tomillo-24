# Raw Data Traceability Specification

## Purpose

Moves the Jenkyn RCBD CSV into the design's raw-data folder and documents its
provenance so every transformation stays reproducible, per AGENTS.md section 4
and the project's data-integrity rules (raw data is never modified).

## Requirements

### Requirement: Raw CSV relocated to datos_crudos/bdca/

The CSV `DBCA_Jenkyn_control_mildeo.csv` MUST be moved from the repository
root to `datos_crudos/bdca/` via `git mv`, preserving its content
byte-for-byte, and MUST NOT be modified.

#### Scenario: File moved with intact content

- GIVEN the CSV at the repository root
- WHEN the file is relocated with `git mv`
- THEN the file exists at `datos_crudos/bdca/DBCA_Jenkyn_control_mildeo.csv`
- AND its byte content is identical to the original

#### Scenario: Rollback preserves the source

- GIVEN a traceability problem with the moved file
- WHEN the change is rolled back
- THEN the CSV is restored to its original location and content via `git mv`

### Requirement: Documentation of source and provenance

AGENTS.md section 4 and the docstrings of `config.py`, `limpiar.py` and
`informe.py` MUST document the bdca raw source (path and schema `plot, trt,
block, yield`) and the provenance rules applying to it.

#### Scenario: Docs updated for bdca

- GIVEN the CSV has been relocated
- WHEN AGENTS.md section 4 and the module docstrings are inspected
- THEN they reference `datos_crudos/bdca/` and its schema
- AND they state that raw values are never modified

### Requirement: Transformation provenance

Every derived variable used in the bdca analysis MUST document its source
variables, formula and reason, and the audit tables MUST record them.

#### Scenario: Derived values documented

- GIVEN the analysis produces derived values (e.g. means, ICC, adjusted p-values)
- WHEN the report and audit tables are inspected
- THEN each transformation states its source variables, formula and reason
