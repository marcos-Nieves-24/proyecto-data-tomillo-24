# Informe de Archive — `analisis-dbca`

## Cambio

- **ID**: `analisis-dbca`
- **Título**: Análisis de rendimiento DBCA (RCBD) — ensayo de control de mildiú de Jenkyn
- **Rama**: `feat/analisis-dbca`
- **Commits verificados**: 138ec7c, 748b574 (PR1 registro+CSV); 160854a (PR2 pipeline/bdca); 440f2cb, 7766f1d (PR3 notebook+página); dd33231, 0855d64 (PR4 auditoría+ejecución)

## Estado final

| Aspecto | Estado |
|---------|--------|
| Especificaciones delta | 4 specs sincronizadas a baseline |
| Requerimientos | 14/14 (PASS) |
| Escenarios | 23/23 (PASS) |
| Tareas | 18/18 marcadas `[x]`, 0 pendientes |
| Verificación | PASS (sin issues CRITICAL) |
| Archivo | CERRADO |

---

## Sincronización de specs delta a baseline

Los cuatro specs delta son specs completas (no deltas parciales) y no existían
specs baseline previas en `openspec/specs/`. Por convención del skill
`sdd-archive` (si la main spec no existe, la delta spec se copia directamente),
cada spec delta pasó a ser la spec baseline del dominio:

| Dominio | Acción | Spec baseline creada |
|---------|--------|----------------------|
| `pipeline-design-registry` | Creada (copia íntegra) | `openspec/specs/pipeline-design-registry/spec.md` |
| `raw-data-traceability` | Creada (copia íntegra) | `openspec/specs/raw-data-traceability/spec.md` |
| `rcbd-yield-analysis` | Creada (copia íntegra) | `openspec/specs/rcbd-yield-analysis/spec.md` |
| `rcbd-reporting` | Creada (copia íntegra) | `openspec/specs/rcbd-reporting/spec.md` |

Los requisitos baseline resultantes: REG-1..REG-3 (6 escenarios), RAW-1..RAW-3
(4 escenarios), YLD-1..YLD-5 (9 escenarios), REP-1..REP-3 (4 escenarios).
Total: 14 requerimientos / 23 escenarios.

Verificación de integridad: `diff -r` entre `openspec/changes/analisis-dbca/specs`
y `openspec/specs` reporta sincronización byte-idéntica (SYNC OK).

## Conservación del registro histórico

Por instrucción explícita del orquestador, el directorio del cambio
`openspec/changes/analisis-dbca/` **se conserva en su lugar** como registro
histórico completo y NO se movió a `openspec/changes/archive/`. No se eliminó
ni modificó ningún artefacto del cambio (proposal, specs, design, tasks,
verify-report).

> Nota de desviación frente al procedimiento estándar del skill `sdd-archive`
> (que mueve el directorio a `openspec/changes/archive/YYYY-MM-DD-{change-name}/`):
> la instrucción del orquestador de preservar el cambio in-situ es la que prevalece.
> El archive no fue destructivo y la fuente de verdad (specs baseline) quedó
> actualizada en `openspec/specs/`.

## Cierre del cambio

- Todas las tareas de `tasks.md` están marcadas `[x]` (18/18); no hay tareas
  pendientes ni checkboxes obsoletos que reconciliar.
- El `verify-report.md` registra PASS en 14/14 requerimientos y 23/23
  escenarios, sin issues CRITICAL. Las advertencias (artefactos no
  deterministas, FutureWarning de seaborn) no bloquean el archive.
- El ciclo SDD para `analisis-dbca` queda completo: planificado, especificado,
  diseñado, implementado, verificado y archivado.

## Fuente de verdad actualizada

Los siguientes specs baseline reflejan ahora el comportamiento implementado:

- `openspec/specs/pipeline-design-registry/spec.md`
- `openspec/specs/raw-data-traceability/spec.md`
- `openspec/specs/rcbd-yield-analysis/spec.md`
- `openspec/specs/rcbd-reporting/spec.md`

## Alcance

Solo se tocó `openspec/`. No se modificaron `pipeline/`, `generar_*.py`,
`bdca/`, `dca/`, `datos_crudos/` ni `pagina/`. No se ejecutó `git add`/`commit`.

## Trazabilidad

- **Engram**: `mem_save` con topic_key `sdd/analisis-dbca/archive-report`, type
  `architecture`, project `proyecto-tomillo`, capture_prompt false.
- **Archivo**: `openspec/changes/analisis-dbca/archive-report.md` (este informe).

## Riesgos

- Ninguno funcional. Única observación: la no utilización del directorio
  `openspec/changes/archive/` (decisión del orquestador para conservar el
  registro in-situ); puede convenir a futuro decidir si ese directorio se
  emplea para otros cambios.
- `tasks.md.bak` permanece en el directorio del cambio (registrado como
  SUGGESTION en el verify-report); se conserva como parte del registro.

## Siguiente paso recomendado

Ninguno — ciclo SDD completo. Listo para el próximo cambio.
