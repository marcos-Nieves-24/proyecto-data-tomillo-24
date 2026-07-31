# Proposal: DBCA (RCBD) Yield Analysis — Jenkyn mildew control trial

**id**: `analisis-dbca`
**titulo**: Add a complete randomized block design (DBCA/RCBD) analysis to the reproducible pipeline, matching the DCA standard, without breaking the existing DCA pipeline.
**resumen**: Introduce `bdca/` as a second experimental design: move the raw Jenkyn CSV into `datos_crudos/bdca/`, parameterize the pipeline via a design registry (`DISENOS` dict + `PIPELINE_DISENO` env var, default `dca`), and build an educational notebook + results + Vercel page for a single-response RCBD analysis (yield).

## Problem / Context

- Raw `DBCA_Jenkyn_control_mildeo.csv` untracked at root: 36×4 (`plot, trt, block, yield`), balanced — trt R/T0/T1/T2 (9 each), block B1–B9 (4 each); no NA/duplicates; yield 4.38–6.54 (mean 5.80).
- Pipeline hardcoded to DCA: `DISENO="dca"` in `config.py`, DCA-only formulas (`supuestos.py`), narratives (`diseno.py`, `informe.py`), generators. No design registry.
- `pagina/bdca/index.html` is a "Próximamente" placeholder; hub lists bdca as future.

## Proposed Solution

- **Registry**: `DISENOS` dict in `config.py`; `DISENO_ACTIVO` from env `PIPELINE_DISENO` (default `dca`), resolved before directory creation. No DCA logic changes.
- **Raw data**: `git mv` CSV → `datos_crudos/bdca/`; update AGENTS.md §4 + docstrings (`config.py`, `limpiar.py`, `informe.py`).
- **Model** (yield, single response): EDA → assumptions → **classical block ANOVA** `yield ~ C(trt)+C(block)` (educational complement) + **LMM** `yield ~ C(trt)` with `groups=block` (REML), report **ICC**; **Tukey HSD** with explicit pair-vs-**R** (reference) interpretation; no Dunnett (no new deps). Report raw means + 95% CI. Document limitation: additivity not testable (1 obs/cell, no interaction).
- **Generators**: DBCA templates in `generar_notebook_pipeline.py` (→ `bdca/analisis_bdca.ipynb`) and `generar_pagina.py` (→ `pagina/bdca/index.html`), single-response structure (not DCA 4-block storytelling). Hub: mark bdca active.
- **N/A**: isolate multivariate, technique ranking, conidia Poisson/NB, C4 controls.

## Scope

### In Scope
- `config.py` design registry + env var
- DBCA modules/narratives (eda/supuestos/modelos/comparaciones/informe)
- `bdca/analisis_bdca.ipynb` + results in `bdca/resultados/` (tablas/figuras/reportes/excel)
- `pagina/bdca/index.html` + hub activation
- Raw CSV move + traceability

### Out of Scope
- Dunnett (no new deps); factorial (`FACTORIAL_Snijders...csv`); susceptibility, dose-response, PCA/clustering, ranking; any DCA behavior/artifact change

## Impact on Existing DCA (must not break)

| Area | Impact | Mitigation |
|------|--------|------------|
| `pipeline/config.py` | Modified | Additive registry; default `dca`; dirs created after design resolution |
| `supuestos.py`, `modelos.py`, `comparaciones.py`, `informe.py`, `diseno.py`, `limpiar.py` | Modified | Design-conditional branches; DCA paths unchanged when env unset |
| `generar_notebook_pipeline.py`, `generar_pagina.py` | Modified | New per-design templates; DCA output unchanged |
| `pagina/` | Modified | bdca page replaces placeholder; hub active |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| DCA regression | Med | Default env = `dca`; smoke-run DCA notebook; 400-line review budget |
| Block model interpretation | Low | Random block (LMM, ICC) + classical ANOVA complement documented |
| CSV move breaks traceability | Low | `git mv` + AGENTS.md §4 + docstrings |
| Additivity untestable | Certain | Documented limitation |

## Decisions

1. LMM random block (REML) + ICC primary; classical block ANOVA as complement.
2. R = reference; Tukey HSD + explicit pairs-vs-R; no Dunnett.
3. Response: raw yield; means + 95% CI.
4. Full deliverable (notebook + results + page), DCA standard.
5. Registry via `DISENOS` dict + env var; DCA untouched.
6. Additivity limitation documented.

## Rollback Plan

- `PIPELINE_DISENO=dca` (or unset) → DCA-only behavior; no DCA code paths removed.
- Revert `git mv` if traceability issues.
- Delete `bdca/` + artifacts; hub placeholder regenerated via `generar_pagina.py --hub`.
- Source CSV preserved byte-for-byte.

## Dependencies

- statsmodels (LMM REML, ANOVA), scipy/pingouin (Tukey) — already in stack. No new packages.

## Work Phases (estimated)

1. Design registry in `config.py` + env resolution (small)
2. CSV move + traceability docs (small)
3. DBCA analysis modules (large)
4. Notebook + page templates + hub (medium)
5. Run pipeline for `bdca` (small)
6. Verify DCA regression + DBCA outputs (small)

## Success Criteria

- [ ] `python -c "import pipeline"` passes; DCA notebook regenerates unchanged
- [ ] `PIPELINE_DISENO=bdca` produces notebook + full results
- [ ] `pagina/bdca/index.html` renders; hub shows bdca active
- [ ] CSV in `datos_crudos/bdca/` with documented traceability
- [ ] Report includes ICC, Tukey pairs vs R, means + 95% CI, additivity limitation
