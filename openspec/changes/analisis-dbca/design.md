# Design: DBCA (RCBD) Yield Analysis — Jenkyn mildew control trial

## Technical Approach

Add `bdca` as a second experimental design without touching DCA logic. A `DISENOS` registry in `pipeline/config.py` enumerates designs and resolves the active one from `PIPELINE_DISENO` (default `dca`, unknown → warning + fallback `dca`) **before** the import-time directory creation. The bdca analysis lives in a new `pipeline/bdca/` subpackage organized by phase (mirroring `pipeline/`), reusing `config.py` export helpers (`guardar_tabla`, `save_figure_pub`, `exportar_excel`, `fijar_semilla`). Generators select per-design templates; the raw CSV moves to `datos_crudos/bdca/` with documented provenance. DCA paths are untouched when `PIPELINE_DISENO` is unset. Maps to proposal §Proposed Solution and all four delta specs.

## Architecture Decisions

### Decision: Design registry and resolution timing

| Option | Tradeoff | Decision |
|---|---|---|
| if/elif on env var | works, but unmaintainable for a 3rd design (factorial) | **`DISENOS` dict registry** |
| Resolve lazily on first use | dirs are created at import (config.py L53-54) — too late | **Resolve at import, before `mkdir` loop** |
| Fail hard on unknown value | breaks batch runs on typos | **`warnings.warn` + fallback to `dca`** |

### Decision: Where bdca analysis code lives

| Option | Tradeoff | Decision |
|---|---|---|
| Conditionals inside existing DCA modules (proposal's letter) | interleaves DCA paths → regression risk, bloats review budget | **New `pipeline/bdca/` subpackage** |
| Single flat `pipeline/bdca.py` | fewer files, but violates per-phase organization rule (config.yaml) | **One module per phase** |

Existing DCA modules receive **docstring-only** updates (`config.py`, `limpiar.py`, `informe.py` — provenance per spec).

### Decision: Statistical model for yield

| Option | Tradeoff | Decision |
|---|---|---|
| Classical ANOVA `yield ~ C(trt)+C(block)` only | ignores block variance; educational | **LMM primary + classical ANOVA complement** |
| Dunnett vs R | fewer comparisons, new dependency | **Tukey HSD (6 pairs) + explicit T0/T1/T2-vs-R interpretation** |
| CLD letters (DCA parity) | not required by spec; DCA CLD helper is hardcoded to `metodo_extraccion` | **Omit CLD; means+95% CI figure with significance annotation** |

### Decision: Generator routing

| Option | Tradeoff | Decision |
|---|---|---|
| Keep single generator, branch templates | one script, additive | **Both generators read `PIPELINE_DISENO` and pick template sets** |
| Hub placeholder overwrites bdca page | placeholder loop writes `disenos[1:]` blindly | **Placeholder only for `estado: futuro` designs; bdca becomes `activo`** |

## Data Flow

```
DBCA_Jenkyn_control_mildeo.csv (datos_crudos/bdca/)
  → pipeline/bdca/cargar.py   (audit: 36 rows, 9/trt, 4/block, no NA/dups)
  → pipeline/bdca/eda.py      (descriptivos + IC95% por trt)
  → pipeline/bdca/supuestos.py (Shapiro-Wilk, Levene, DW sobre residuos del ANOVA de bloques)
  → pipeline/bdca/modelos.py  (ANOVA clásico + LMM REML con ICC)
  → pipeline/bdca/comparaciones.py (Tukey 6 pares + interpretación vs R)
  → pipeline/bdca/informe.py  → bdca/resultados/{tablas,figuras,reportes,excel}
generar_notebook_pipeline.py → bdca/analisis_bdca.ipynb
generar_pagina.py             → pagina/bdca/index.html  + hub activo
```

## File Changes

| File | Action | Description |
|---|---|---|
| `pipeline/config.py` | Modify | Add `DISENOS` + `_resolver_diseno_activo()`; derive paths from `DISENOS[DISENO]`; docstring provenance |
| `pipeline/bdca/{__init__,cargar,eda,supuestos,modelos,comparaciones,informe}.py` | Create | Phase modules for single-response RCBD analysis |
| `pipeline/limpiar.py`, `pipeline/informe.py` | Modify | Docstring provenance of `datos_crudos/bdca/` (no behavior change) |
| `generar_notebook_pipeline.py` | Modify | `construir_notebook_dca()` (renamed) + `construir_notebook_bdca()`; env-driven `SALIDA` |
| `generar_pagina.py` | Modify | `--diseno` choices `dca,bdca`; `main_bdca()`; hub placeholder only for `futuro` |
| `DBCA_Jenkyn_control_mildeo.csv` | Move | `mv` + `git add` → `datos_crudos/bdca/` (byte-for-byte; file is untracked, so `git mv` is not usable) |
| `AGENTS.md` | Modify | §4 documents bdca source + schema `plot, trt, block, yield` |
| `bdca/analisis_bdca.ipynb`, `bdca/resultados/*`, `pagina/bdca/index.html` | Create | Generated deliverables |

## Interfaces / Contracts

```python
# pipeline/config.py
DISENOS = {
    "dca":  {"nombre": ..., "dir_diseno": RAIZ/"dca",  "dir_crudos": RAIZ/"datos_crudos"/"dca"},
    "bdca": {"nombre": ..., "dir_diseno": RAIZ/"bdca", "dir_crudos": RAIZ/"datos_crudos"/"bdca"},
}
DISENO = _resolver_diseno_activo()   # env PIPELINE_DISENO, default "dca", warn+fallback
DIR_DISENO, DIR_CRUDOS = DISENOS[DISENO]["dir_diseno"], DISENOS[DISENO]["dir_crudos"]
# existing derived constants (EXCEL_CRUDO, DIR_TABLAS, ...) unchanged derivation

# pipeline/bdca/modelos.py — signatures
def anova_bloques(df) -> dict        # smf.ols("yield ~ C(trt) + C(block)"), anova_lm typ=2
def lmm_bloques(df) -> dict          # smf.mixedlm("yield ~ C(trt)", groups=df["block"], reml=True)
                                     # ICC = var_block/(var_block+var_residual)  (cov_re[0,0]/scale)
def posthoc_tukey(df) -> pd.DataFrame  # pairwise_tukeyhsd; columna "vs_referencia_R" por par
```

LMM/ICC and Tukey follow the proven DCA patterns (`analisis_sensibilidad_lmm`, `_tukey_hsd`); no new dependencies. Balanced design (1 obs/cell) → no interaction term; additivity documented as untestable.

## Testing Strategy

Testing is `none`; verification is pipeline-execution smoke tests (config.yaml `verify`):

| Layer | What to Test | Approach |
|---|---|---|
| Smoke | `import pipeline.config` with env unset / `bdca` / `foo` | `python -c "import pipeline.config"` (must import the config module directly: `pipeline/__init__.py` does not import `config.py`); assert warning+fallback on `foo` |
| Regression guard | DCA regeneration unchanged | Run generators with env unset; `git diff` clean on `dca/analisis_dca.ipynb`, `pagina/dca/index.html` |
| E2E | bdca notebook executes end-to-end | `PIPELINE_DISENO=bdca` run; audit tables: 36 rows, 9/trt, 4/block, no NA/dups; all 5 result dirs populated |
| Isolation | bdca run touches nothing under `dca/` | `git status` before/after + timestamps |

## Threat Matrix

`N/A` — no routing, shell, subprocess, VCS/PR automation, executable-file classification, or process-integration boundary. (Moving the raw CSV is a one-time implementation step, not a runtime boundary; no RED tests required.)

## Migration / Rollout

Phased: (1) registry in `config.py` (additive, default `dca`) → (2) `mv` CSV to `datos_crudos/bdca/` + `git add` + traceability docs → (3) `pipeline/bdca/` modules → (4) generator templates + hub → (5) run bdca pipeline → (6) verify DCA regression. Rollback: unset env → DCA-only; revert the CSV move (restore original path); delete `bdca/` + `pagina/bdca/index.html`; regenerate hub. No data migration.

## Open Questions

- [ ] Notebook execution for verification: `nbconvert --execute` or manual run?
- [ ] Keep raw treatment labels (`R/T0/T1/T2`) or add Spanish display labels in figures only?
- [ ] Confirmed: no CLD letters for bdca post-hoc (spec omits them).
