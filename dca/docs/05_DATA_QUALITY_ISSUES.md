# Data Quality Issues — Proyecto Tomillo × Fusarium spp.

> **Audit date:** 2026-07-28
> **Source file:** `datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx`
> **Status:** Open — each issue requires investigator resolution before analysis

---

## Issues Summary

| ID | Severity | Category | Issue |
|----|----------|----------|-------|
| DQ01 | **CRITICAL** | Structural | Consensus sheet (DATOS_CONSESUS) is 96.8% empty |
| DQ02 | **CRITICAL** | Structural | Soxhlet & ultrasonido data only exist in per-method sheets, not in consensus |
| DQ03 | **CRITICAL** | Completeness | No conidia data anywhere in DATOS_CONSESUS |
| DQ04 | **CRITICAL** | Completeness | No derived inhibition variables in DATOS_CONSESUS |
| DQ05 | **CRITICAL** | Design | Control (C4) measured only once per isolate, shared across 3 reps |
| DQ06 | **CRITICAL** | Completeness | Intermediate concentrations (1, 0.2 mg/mL) only tested in maceración |
| DQ07 | **HIGH** | Consistency | ULTRASONIDO sheet contains "SOXHLET (5 días)" label — copy-paste error |
| DQ08 | **HIGH** | Consistency | Variable column layout between sheets (columns shifted) |
| DQ09 | **HIGH** | Ambiguity | Extra unnamed columns in every per-method sheet |
| DQ10 | **HIGH** | Ambiguity | id and bloque columns entirely empty in DATOS_CONSESUS |
| DQ11 | **HIGH** | Completeness | C2 and C3 concentration data missing for 6 isolates in MACERACIÓN |
| DQ12 | **MODERATE** | Consistency | Grupo classification differs between MAC/SOX and ULTRASONIDO |
| DQ13 | **MODERATE** | Ambiguity | "RADIO en mm" — is this radius or diameter? |
| DQ14 | **MODERATE** | Naming | Isolate naming inconsistency: FU2 vs FU2 (UCMU21) |
| DQ15 | **LOW** | Documentation | "#" column in ULTRASONIDO has non-contiguous numbering |
| DQ16 | **LOW** | Data | Control growth values identical across replicates within each isolate |
| DQ17 | **LOW** | Completeness | Extra columns only have data for reps 1–2, never rep 3 |

---

## DQ01 — Consensus sheet is critically incomplete

**Observation:** `DATOS_CONSESUS` has 558 rows × 11 columns but only 2 columns
(`metodo_extraccion`, `aislamiento`, `concentracion`, `rep_biologica`, `control`)
are fully populated. `crecimiento_micelio_mm` has data in only 186/558 rows (maceración only).
`ptj_inhibicion`, `conidia_ml`, and `ptj_reduccion_conidias` are **entirely empty**.

**Risk:** If analysis uses DATOS_CONSESUS as the primary data source, it will find
no conidia data, no inhibition percentages, and no data for 2 of 3 extraction methods.

**Resolution:** Investigator must confirm whether the per-method sheets are the
authoritative data source, and whether DATOS_CONSESUS was simply never populated.

---

## DQ02 — Soxhlet and ultrasonido data only in per-method sheets

**Observation:** `crecimiento_micelio_mm` is non-null only for `maceracion` (186 rows).
For `soxhlet` and `ultrasonido`, all 372 rows have `NaN` for all measurement columns.

**Risk:** The "consensus" sheet is misleading — it suggests data exist for all three
methods when they do not.

**Resolution:** The per-method sheets (MACERACIÓN, SOXHLET, ULTRASONIDO) contain the
actual raw data. The pipeline must extract and reshape these sheets, not rely on
DATOS_CONSESUS.

---

## DQ03 — No conidia data in DATOS_CONSESUS

**Observation:** `conidia_ml` and `ptj_reduccion_conidias` have 0 non-null values
across all 558 rows.

**Impact:** Objective 3 (conidia production) cannot be addressed from the consensus sheet.
Data exists in per-method sheets but must be extracted.

---

## DQ04 — No derived inhibition in DATOS_CONSESUS

**Observation:** `ptj_inhibicion` has 0 non-null values.

**Impact:** No pre-computed inhibition percentages exist in the consensus sheet.
Inhibition must be calculated from raw growth data using control values.

---

## DQ05 — Control measured only once per isolate (pseudoreplication risk)

**Observation:** In all three per-method sheets, C4 (control) has only 31/93 non-null
values. This means the control mycelial radius or conidia count was measured **once
per isolate** and then shared across all 3 replicates.

**Risk:** 
- If the same control value is used to calculate %INH for all 3 replicates,
  the inhibition percentages within an isolate are **not independent** — they share
  a common denominator.
- This creates a form of **technical pseudoreplication**: the 3 biological replicates
  at the treatment level are compared against a single control measurement.

**Statistics implication:** Standard ANOVA treats each replicate's inhibition as
independent, but the shared control introduces correlation. A mixed model
(isolate as random intercept) may partially address this, but the effective
sample size for control comparisons is 31, not 93.

**Resolution needed:** 
1. Confirm that each isolate had one control plate shared across 3 treatment plates.
2. Clarify whether the 3 treatment replicates are truly independent biological
   replicates (separate inoculations) or technical repeats of the same plate.
3. If all 3 treatment replicates used the same control, this must be modeled
   explicitly (e.g., using a paired structure or modeling raw growth, not %INH).

---

## DQ06 — Intermediate concentrations only in maceración

**Observation:** 
- MACERACIÓN tests 4 concentration levels: 0 (C4), 0.2 (C3), 1 (C2), 5 (C1) mg/mL
- SOXHLET tests only 2 levels: 0 (C4) and 5 (C1) mg/mL
- ULTRASONIDO tests only 2 levels: 0 (C4) and 5 (C1) mg/mL
- DATOS_CONSESUS only has 0 and 5 mg/mL

**Impact on dose-response (Objective 4):**
- Dose-response curves (EC50, EC90) require ≥3 concentration levels.
- MACERACIÓN has 4 levels → dose-response possible.
- SOXHLET and ULTRASONIDO have only 2 levels → **cannot fit 4-parameter logistic**
  or estimate EC50/EC90 for these methods.
- If the 0.2 and 1 mg/mL data exist for SOXHLET and ULTRASONIDO but were not
  transcribed, the investigator must provide them.

---

## DQ07 — Copy-paste error in ULTRASONIDO sheet headers

**Observation:** Row 1 of the ULTRASONIDO sheet contains `"SOXHLET (5 días)"`
in the RADIO and CONIDIA section headers. This is clearly a copy-paste error
from the SOXHLET sheet.

**Impact:** Low for data analysis (the data values themselves are not affected),
but reflects on documentation quality and reproducibility.

---

## DQ08 — Inconsistent column layout across method sheets

**Observation:** The three per-method sheets have different column structures:

| Sheet | Isolate col | Grupo col | Replicate col | C1 | C2 | C3 | C4 | %INH C1 | Extra |
|-------|------------|-----------|--------------|----|----|----|----|---------|-------|
| MAC | Col 0 | Col 1 | Col 2 | Col 3 | Col 4 | Col 5 | Col 6 | Col 7 | Col 8 |
| SOX | Col 0 | Col 1 | Col 2 | Col 3 | Col 4 (C4!) | — | — | Col 5 | Col 6 |
| ULT | Col 1 | Col 2 | Col 3 | Col 4 | — | — | Col 5 | Col 6 | Col 7 |

Additionally, ULTRASONIDO has an extra `#` column (Col 0) that shifts all other
columns by 1, and the conidia block starts at Col 8 in SOXHLET, Col 9 in MACERACIÓN,
Col 8 in ULTRASONIDO.

**Impact:** Any code that parses these sheets must handle the structural differences.
A single unified parser for all three will produce errors if column indices
are hardcoded.

---

## DQ09 — Extra unnamed columns (ambiguous)

**Observation:** Every per-method sheet has unnamed numeric columns whose purpose
is unclear:

| Sheet | Radio extra | Conidia extra | Data pattern |
|-------|------------|--------------|--------------|
| MACERACIÓN | Col 8 | Col 15 | Present in reps 1–2, absent in rep 3 |
| SOXHLET | Col 6 | Col 11 | Same pattern |
| ULTRASONIDO | Col 7 | Col 12 | Same pattern |

**Values (radio extra, MACERACIÓN):** Range 0–100%, mean ~47%, 62/93 non-null.
Correlated with but not equal to `%INH C1`.

**Hypotheses (unconfirmed):**
1. %INH for C2 or C3 (but values don't match the formula with available data)
2. %INH calculated with a different control reference (e.g., a control from a
   different replicate or an external reference)
3. A different normalization or correction factor

**Resolution:** Must consult the investigator to determine what these columns
represent before including them in analysis.

---

## DQ10 — `id` and `bloque` entirely empty

**Observation:** Both columns in DATOS_CONSESUS have 0/558 non-null values.

**Impact:**
- Cannot verify if the dataset is blocked or randomized.
- Cannot confirm whether `Grupo para actividad` from the per-method sheets
  corresponds to `bloque`.
- No unique row identifier — difficult to merge or trace individual observations.

---

## DQ11 — C2 and C3 missing for 6 isolates in MACERACIÓN

**Observation:** C2 (1 mg/mL) and C3 (0.2 mg/mL) in MACERACIÓN have 75/93
non-null values, meaning ~6 isolates are missing these concentrations.

**Affected isolates:** Need to check which specific isolates have missing C2/C3
data. This reduces the dose-response resolution for those isolates.

**Resolution:** Confirm whether these were not tested or if the data is missing
from the sheet.

---

## DQ12 — Grupo classification inconsistency in ULTRASONIDO

**Observation:** In MACERACIÓN and SOXHLET, column 1 contains `Grupo X - HXX`
labels. In ULTRASONIDO, column 2 contains isolate names instead of group labels.

**Impact:** The group classification cannot be consistently extracted from
ULTRASONIDO using the same parsing logic as MACERACIÓN and SOXHLET.

---

## DQ13 — "RADIO en mm": radius or diameter?

**Observation:** The header uses `RADIO` which is Spanish for "radius". However,
biological growth measurements in poison plate assays are typically reported as
**colony diameter**, not radius.

**Values observed:** 0–75 mm. A typical 90 mm Petri dish has a usable diameter of
~80–85 mm. If these are diameters, they range from 0 (complete inhibition) to
75 mm (near full plate coverage). This is biologically plausible for diameter.
If these are radii, the colony diameters would be 0–150 mm, which exceeds a
standard Petri dish — **implausible**.

**Conclusion:** Despite the header saying "RADIO", the values are likely colony
**diameters** in mm.

**Resolution:** Confirm with the investigator. This affects %INH calculation only
if a ratio-based formula is used (ratio cancels out, so %INH is unaffected).
But absolute growth values in mm need correct interpretation.

---

## DQ14 — Isolate naming inconsistency

**Observation:** Isolate `FU2` appears as both `FU2` and `FU2 (UCMU21)` in
different rows/locations in the dataset.

**Impact:** Depending on string matching, these could be treated as two different
isolates when they are the same. Must be standardized.

**Other naming of note:** `FUSARIUM JULIAN H20` and `FUSARIUM MARCE 1.2` are
named more descriptively than the HC/H codes. This could complicate grouping.

---

## DQ15 — ULTRASONIDO "#" column non-contiguous numbering

**Observation:** The `#` column has values: 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37.

**Impact:** Minimal. The numbers seem to be sequential but skip some values
(3, 14, 15, 19, 35 missing). This could be a remnant from a larger master list
or data entry order.

---

## DQ16 — Control growth values identical across replicates

**Observation:** In DATOS_CONSESUS, for each isolate at 0 mg/mL (control),
all 3 replicates have the **identical** growth value (std = 0 for all isolates).

**Example:** HC3 control = 68 mm for reps 1, 2, and 3.

**Interpretation:**
- The same control measurement was entered for all 3 replicates.
- This is consistent with DQ05 (single control measurement shared across reps).
- When used to calculate % inhibition, the 3 %INH values per isolate share a
  common control value, making them statistically dependent.

---

## DQ17 — Extra columns only in reps 1–2

**Observation:** The extra unnamed columns consistently have 62/93 non-null
entries, meaning they are populated for replicates 1 and 2 (2 × 31 = 62) but
never for replicate 3.

**Pattern:** Always present in RÉPLICA 1 and RÉPLICA 2, absent in RÉPLICA 3.

**Possible explanations:**
1. The extra column is a secondary calculation not needed for rep 3
2. The extra column was added after rep 3 was already measured
3. Data entry was incomplete

---

## DQ18 — Conidia values are log-transformed, not raw counts

**Observation:** The header clearly states "Log Conidias/ml". Values range from
0 to 8.65. A value of 7.0 corresponds to 10⁷ = 10,000,000 conidia/mL.

**Impact:** Since values are already log-transformed:
- Cannot model as Poisson or negative binomial (not raw counts).
- The data is already log-normal or approximately normal on the log scale.
- A linear model on the log scale is the natural approach.
- The "reduction" percentages computed from log values are not standard
  percent reductions of raw counts. They represent (1 − log₁₀(T) / log₁₀(C)) × 100,
  which is NOT a conventional inhibition percentage.

**Important:** `%INH C1` for conidia likely uses the formula:
  %INH = (1 − log₁₀(treatment) / log₁₀(control)) × 100

This is NOT the same as percent reduction in raw conidia count. The investigator
should confirm whether this is intentional or whether the inhibition should be
computed from back-transformed counts.

---

## Notes for Analysis Pipeline

1. **Do NOT use DATOS_CONSESUS as primary source.** Parse per-method sheets instead.
2. **Reconstruct C4 control values** for each isolate (once per isolate, shared across reps).
3. **Calculate %INH** from raw growth, not from pre-computed values.
4. **Model raw growth** (not %INH) to avoid denominator-induced correlation, or
   use a mixed model that accounts for shared controls.
5. **Conidia** are log₍₁₀₎-transformed — handle accordingly.
6. **Dose-response (EC50)** is only possible for MACERACIÓN data (4 concentrations).
7. **Validate** the mapping between "RADIO en mm" values and actual colony
   measurement protocol.
