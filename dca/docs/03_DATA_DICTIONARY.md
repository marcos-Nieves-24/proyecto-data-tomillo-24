# Data Dictionary — Proyecto Tomillo × Fusarium spp.

> **Document version:** Audit, 2026-07-28
> **Source file:** `datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx`
> **Status:** Preliminary — based on structural inspection, pending investigator confirmation

---

## 1. Sheet Inventory

| Sheet | Rows | Cols | Content |
|-------|------|------|---------|
| `RENDIMIENTOS` | 9 | 4 | Extraction yield data |
| `DATOS_CONSESUS` | 558 | 11 | Consensus dataset (sparse — see notes) |
| `MACERACIÓN` | 95 | 16 | Maceration method — full experimental data |
| `SOXHLET` | 95 | 12 | Soxhlet method — full experimental data |
| `ULTRASONIDO` | 95 | 13 | Ultrasound method — full experimental data |

---

## 2. Sheet: `RENDIMIENTOS`

**Purpose:** Extraction yield comparison across 3 techniques.

| Column | Type | Description | Values |
|--------|------|-------------|--------|
| `TOMILLO - METANOL` | str | Extraction technique | `MACERACIÓN`, `SOXHLET`, `ULTRASONIDO` |
| `PESO MATERIAL SECO (g)` | int | Dry plant material weight used | 30 (SOXHLET), 40 (MACERACIÓN, ULTRASONIDO) |
| `PESO EXTRACTO OBTENIDO (g)` | float | Mass of extract recovered | 4.5–15.2 |
| `RENDIMIENTO` | float | Yield = (extract / dry material) × 100 | 11.3–50.6 |

**Structure:** 3 methods × 3 replicates = 9 rows. Complete.

---

## 3. Sheet: `DATOS_CONSESUS`

**Purpose:** Consensus/compiled dataset intended to contain all experimental observations.

**⚠️ Critical: this sheet is incomplete.** Only mycelial growth data for the maceration method is populated. All other columns (`ptj_inhibicion`, `conidia_ml`, `ptj_reduccion_conidias`) are entirely empty.

| Column | Type | Description | Completeness |
|--------|------|-------------|-------------|
| `id` | float | (Empty) — intended row identifier | 0/558 non-null |
| `bloque` | float | (Empty) — intended experimental block | 0/558 non-null |
| `metodo_extraccion` | str | Extraction technique | 558/558: `maceracion`, `soxhlet`, `ultrasonido` |
| `aislamiento` | str | Fusarium isolate code | 558/558: 31 unique isolates |
| `concentracion` | str | Extract concentration | 558/558: `0 mg/ml`, `5 mg/ml` |
| `rep_biologica` | int | Biological replicate number | 558/558: 1, 2, 3 |
| `control` | str | Whether this is a control | 558/558: `Si`, `No` |
| `crecimiento_micelio_mm` | float | Mycelial growth diameter (mm) | 186/558 — only maceración method |
| `ptj_inhibicion` | float | Inhibition percentage (derived) | **0/558** — entirely empty |
| `conidia_ml` | float | Conidia per mL | **0/558** — entirely empty |
| `ptj_reduccion_conidias` | float | Conidia reduction percentage (derived) | **0/558** — entirely empty |

**Structure:**
- 31 isolates × 3 methods × 2 concentrations × 3 reps = 558 rows
- Only maceración has growth data (31 × 2 × 3 = 186)
- Soxhlet and ultrasonido growth data are entirely absent

### 3.1 Concentration Codes

| Label | Actual concentration | Notes |
|-------|---------------------|-------|
| `0 mg/ml` | 0 (control) | Corresponds to `C4 (CTROL)` in per-method sheets |
| `5 mg/ml` | 5 mg/mL | Corresponds to `C1 (5 mg/ml)` in per-method sheets |

**⚠️ Only 2 of the 3–4 planned concentrations are present in DATOS_CONSESUS.**

---

## 4. Sheet: `MACERACIÓN`

**Purpose:** Raw experimental data for the maceration extraction technique.

**Structure:** Multi-header format with 3 header rows and data in row 3+.

### 4.1 RADIO block (mycelial growth, mm)

| Column | Header label | Content | Completeness |
|--------|-------------|---------|-------------|
| Col 0 | `Aislado` | Isolate code (filled once per 3-rep group) | 31/93 |
| Col 1 | `Grupo para actividad` | Experimental group classification | 31/93 |
| Col 2 | `ACTIVIDAD (…) (RADIO en mm)` | Replicate label: `RÉPLICA 1, 2, 3` | 93/93 |
| Col 3 | `C1 \n(5 mg/ml)` | Mycelial radius at 5 mg/mL (mm) | 93/93 |
| Col 4 | `C2 \n(1 mg/ml)` | Mycelial radius at 1 mg/mL (mm) | 75/93 |
| Col 5 | `C3 \n(0.2 mg/ml)` | Mycelial radius at 0.2 mg/mL (mm) | 75/93 |
| Col 6 | `C4 \n(CTROL)` | Control (0 mg/mL) mycelial radius (mm) | 31/93 |
| Col 7 | `%INH\nC1` | Inhibition % at C1 relative to C4 | 93/93 |
| Col 8 | *(unnamed)* | **Ambiguous** — extra derived measure | 62/93 |

### 4.2 CONIDIA block (log conidia/mL)

| Column | Header label | Content | Completeness |
|--------|-------------|---------|-------------|
| Col 9 | `ACTIVIDAD (…) (PRODUCCIÓN CONIDIAS) (Log Conidias/ml)` | Replicate label: `RÉPLICA 1, 2, 3` | 93/93 |
| Col 10 | `C1 \n(5 mg/ml)` | Log conidia/mL at 5 mg/mL | 93/93 |
| Col 11 | `C2 \n(1 mg/ml)` | Log conidia/mL at 1 mg/mL | 75/93 |
| Col 12 | `C3 \n(0.2 mg/ml)` | Log conidia/mL at 0.2 mg/mL | 75/93 |
| Col 13 | `C4 \n(CTROL)` | Log conidia/mL control | 31/93 |
| Col 14 | `%INH\nC1` | Conidia reduction % at C1 | 93/93 |
| Col 15 | *(unnamed)* | **Ambiguous** — extra derived measure | 62/93 |

### 4.3 Concentration mapping (MACERACIÓN)

| Code | Concentration | RADIO column | CONIDIA column |
|------|-------------|-------------|----------------|
| C1 | 5 mg/mL | Col 3 | Col 10 |
| C2 | 1 mg/mL | Col 4 | Col 11 |
| C3 | 0.2 mg/mL | Col 5 | Col 12 |
| C4 | 0 mg/mL (control) | Col 6 | Col 13 |

---

## 5. Sheet: `SOXHLET`

**Purpose:** Raw experimental data for the Soxhlet extraction technique.

### 5.1 Structure differences vs MACERACIÓN

| Feature | MACERACIÓN | SOXHLET |
|---------|-----------|---------|
| Concentrations in RADIO | C1, C2, C3, C4 | C1, C4 only |
| Concentrations in CONIDIA | C1, C2, C3, C4 | C1, C4 only |
| RADIO columns | Cols 3–8 | Cols 3–6 |
| CONIDIA columns | Cols 10–15 | Cols 8–11 |
| Extra columns exist | Yes | Yes |

### 5.2 Column mapping

| Position | Label | Content | Completeness |
|----------|-------|---------|-------------|
| Col 0 | `Aislado` | Isolate code | 31/93 |
| Col 1 | `Grupo para actividad` | Experimental group | 31/93 |
| Col 2 | `ACTIVIDAD (…) (RADIO en mm)` | Replicate label | 93/93 |
| Col 3 | `C1 \n(5 mg/ml)` | Mycelial radius at 5 mg/mL | 93/93 |
| Col 4 | `C4 \n(CTROL)` | Control radius | 31/93 |
| Col 5 | `%INH\nC1` | Inhibition % | 93/93 |
| Col 6 | *(unnamed)* | **Ambiguous extra** | 62/93 |
| Col 7 | `… (CONIDIAS)` | Replicate label (conidia block) | 93/93 |
| Col 8 | `C1 \n(5 mg/ml)` | Log conidia/mL at 5 mg/mL | 93/93 |
| Col 9 | `C4 \n(CTROL)` | Log conidia/mL control | 31/93 |
| Col 10 | `%INH\nC1` | Conidia reduction % | 93/93 |
| Col 11 | *(unnamed)* | **Ambiguous extra** | 62/93 |

---

## 6. Sheet: `ULTRASONIDO`

**Purpose:** Raw experimental data for the ultrasound extraction technique.

### 6.1 Structure differences vs other sheets

| Feature | MACERACIÓN / SOXHLET | ULTRASONIDO |
|---------|---------------------|-------------|
| Col 0 | Aislado | `#` (sequential number) |
| Col 1 | Grupo | Aislado (isolate code) |
| Col 2 | Replicate | Grupo para actividad |
| Col 3 | C1 (RADIO) | Replicate label |
| Col 4 | C2 (RADIO) | C1 (5 mg/mL) RADIO |
| Col 5 | C3 (RADIO) | C4 (CTRL) RADIO |
| Col 6 | C4 (CTRL) RADIO | %INH C1 RADIO |
| Col 7 | %INH C1 RADIO | **Ambiguous extra** (RADIO) |
| ⚠️ | **Row 1 header** says "SOXHLET (5 días)" | **Copy-paste error — should be "ULTRASONIDO (5 días)"** |

### 6.2 Concentrations

Same as SOXHLET: only C1 (5 mg/mL) and C4 (control). No intermediate C2 (1 mg/mL) or C3 (0.2 mg/mL).

### 6.3 Column mapping

| Col | Label | Content | Completeness |
|-----|-------|---------|-------------|
| 0 | `#` | Sequential isolate number (2–37, non-contiguous) | 31/93 |
| 1 | `Aislado` | Isolate code | 31/93 |
| 2 | `Grupo para actividad` | Experimental group | 31/93 |
| 3 | `… (RADIO en mm)` | Replicate label | 93/93 |
| 4 | `C1 \n(5 mg/ml)` | Mycelial radius at 5 mg/mL | 93/93 |
| 5 | `C4 \n(CTROL)` | Control radius | 31/93 |
| 6 | `%INH\nC1` | Inhibition % | 93/93 |
| 7 | *(unnamed)* | **Ambiguous extra** (RADIO) | 62/93 |
| 8 | `… (CONIDIAS)` | Replicate label (conidia) | 93/93 |
| 9 | `C1 \n(5 mg/ml)` | Log conidia/mL at 5 mg/mL | 93/93 |
| 10 | `C4 \n(CTROL)` | Log conidia/mL control | 31/93 |
| 11 | `%INH\nC1` | Conidia reduction % | 93/93 |
| 12 | *(unnamed)* | **Ambiguous extra** (conidia) | 62/93 |

---

## 7. Variable Classification

### 7.1 Experimental (measured) variables

| Variable | Sheet(s) | Nature |
|----------|---------|--------|
| Mycelial growth radius (mm) | MAC/SOX/ULT — Cols C1–C4 | Direct measurement |
| Log conidia/mL | MAC/SOX/ULT — Cols C1–C4 | Measured (log-transformed by lab) |
| Extract weight (g) | RENDIMIENTOS | Direct measurement |
| Dry material weight (g) | RENDIMIENTOS | Direct measurement |

### 7.2 Derived variables

| Variable | Formula (assumed) | Source |
|----------|------------------|--------|
| `RENDIMIENTO` | (extract weight / dry weight) × 100 | RENDIMIENTOS |
| `%INH C1` (RADIO) | (1 − C1 / C4) × 100 | MAC/SOX/ULT |
| `%INH C1` (CONIDIA) | (1 − C1_conidia / C4_conidia) × 100 | MAC/SOX/ULT |
| `ptj_inhibicion` | (derived, formula unknown) | DATOS_CONSESUS (all empty) |
| `ptj_reduccion_conidias` | (derived, formula unknown) | DATOS_CONSESUS (all empty) |

### 7.3 Classification/design variables

| Variable | Role | Levels |
|----------|------|--------|
| `metodo_extraccion` | Fixed factor | 3: maceracion, soxhlet, ultrasonido |
| `aislamiento` | Fixed factor (or random?) | 31 isolates |
| `concentracion` | Fixed factor | 2–4 levels depending on sheet |
| `rep_biologica` | Biological replicate | 3 per combination |
| `control` | Indicator | Si (control), No (treatment) |
| `Grupo para actividad` | Experimental group (1–4) | 4 groups, see Section 8 |
| `#` | Sequential number | ULTRASONIDO only |

---

## 8. Experimental Groups (`Grupo para actividad`)

The `Grupo para actividad` column assigns each isolate to one of four groups,
which may correspond to experimental batches, plates, or processing days.

**Group composition (31 isolates):**

| Group | Isolates | Count |
|-------|---------|-------|
| Grupo 1 | HC5, HC10, HC15, HC16, HC19, HC20, HC23, HC27, HC28 | 9 |
| Grupo 2 | H1N, H2N, H3N, H4N, H5N, H6N, H9N, FU1, FU2(UCMU21) | 9 |
| Grupo 3 | HC3, HC17, HC26, H8N, HC9, FU2(var), FUSARIUM M.1.2, H8N | ~8 |
| Grupo 4 | HC6, H4B, H6B, H4G, H8G, H11G, HC9, FUS. JULIAN H20 | ~8 |

**⚠️ Note:** Groups are NOT consistent across method sheets. In ULTRASONIDO,
col 2 contains isolate names instead of group labels.

---

## 9. Ambiguous / Unknown Columns

### 9.1 Extra column in RADIO block

In all three method sheets, the RADIO block has an unnamed column after `%INH C1`:
- MACERACIÓN: Col 8 (62/93 non-null)
- SOXHLET: Col 6 (62/100 non-null)
- ULTRASONIDO: Col 7 (62/100 non-null)

**Hypotheses:** Unknown. Could be:
- % inhibition for a different concentration not shown in the column header
- % inhibition calculated using a different control value
- A different measurement or derived index entirely

**Resolution needed:** Consult investigator.

### 9.2 Extra column in CONIDIA block

- MACERACIÓN: Col 15 (62/93 non-null)
- SOXHLET: Col 11 (62/93 non-null)
- ULTRASONIDO: Col 12 (62/93 non-null)

Same uncertainty as above.

### 9.3 `ptj_inhibicion` and `ptj_reduccion_conidias` in DATOS_CONSESUS

These columns are entirely empty in `DATOS_CONSESUS`. They exist in the consensus
sheet schema but contain no data. In the per-method sheets, `%INH C1` is present
and populated. Whether `ptj_inhibicion` should equal `%INH C1` is unclear.

### 9.4 `id` and `bloque` in DATOS_CONSESUS

Both entirely empty. `id` was likely intended as a unique row identifier. `bloque`
may have been intended to encode experimental blocks (possibly corresponding to
`Grupo para actividad` from the per-method sheets).

---

## 10. Units

| Variable | Unit | Notes |
|----------|------|-------|
| Mycelial growth | mm (radius or diameter?) | Header says "RADIO" (radius) |
| Conidia | Log₁₀(conidia/mL) | Lab-transformed, not raw counts |
| Concentration | mg/mL | Weight/volume |
| Yield | % (w/w) | (extract / dry material) × 100 |
| Inhibition | % | (1 − treatment / control) × 100 |
| Conidia reduction | % | Same formula |
