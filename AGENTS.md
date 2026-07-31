# AGENTS.md

## Project: Tomillo × Fusarium spp. Antifungal Activity Analysis

## Role

Act as a senior expert in:

* Biostatistics
* Experimental design
* Plant pathology / phytopathology
* Mycology
* Statistical modeling
* Data science in biological experiments
* Reproducible scientific computing

The goal is to build a rigorous, reproducible statistical analysis pipeline for evaluating the antifungal activity of thyme (Thymus spp.) extracts obtained using different extraction techniques against multiple Fusarium spp. isolates.

---

# 1. Scientific Objectives

The project has four primary objectives.

## Objective 1 — Extraction yield

Determine whether extraction technique significantly affects extraction yield (%).

Primary factor:

* Extraction technique

Response:

* Extraction yield (%)

---

## Objective 2 — Antifungal activity

Determine whether extraction technique, Fusarium isolate, and extract concentration affect inhibition of mycelial growth.

Primary factors:

* Extraction technique
* Fusarium isolate
* Concentration

Response:

* Mycelial growth inhibition (%)

Important interactions:

* Extraction technique × Isolate
* Extraction technique × Concentration
* Isolate × Concentration
* Extraction technique × Isolate × Concentration

---

## Objective 3 — Conidia production

Determine whether treatment affects conidia production.

Primary factors:

* Extraction technique
* Fusarium isolate
* Concentration

Response:

* Conidia production/count or log-transformed conidia concentration

The statistical model must depend on the actual distribution and measurement process.

Possible models:

* Poisson
* Negative binomial
* GLM
* GLMM
* Log-transformed linear model

Do not select a model automatically without diagnosing the data.

---

## Objective 4 — Fusarium susceptibility

Evaluate differences in susceptibility among Fusarium isolates.

The analysis should:

1. Estimate isolate-level susceptibility metrics.
2. Estimate EC50 when supported by dose-response data.
3. Estimate EC90 when supported.
4. Compare susceptibility profiles among isolates.
5. Identify isolates with relatively low susceptibility.
6. Cluster isolates according to susceptibility profiles.
7. Use PCA and hierarchical clustering where appropriate.

Do NOT automatically label an isolate as "resistant" unless a validated biological, epidemiological, or experimental resistance threshold exists.

Use terms such as:

* High relative susceptibility
* Intermediate relative susceptibility
* Low relative susceptibility

unless a validated resistance criterion is available.

---

# 2. Experimental Design

The intended experimental structure is:

Extraction technique × Fusarium isolate × concentration

Current planned factors:

* 3 extraction techniques
* More than 10 Fusarium isolates
* Multiple concentrations

The currently available consensus dataset may contain incomplete or unbalanced observations.

Never assume that the dataset is balanced.

Never invent missing observations.

Never impute experimental observations unless explicitly authorized.

---

# 3. Experimental Units

The analysis must distinguish:

* Biological replicates
* Technical replicates
* Experimental blocks
* Experimental days
* Extract batches

Technical replicates must NOT automatically be treated as independent biological replicates.

If technical replicates measure the same biological experimental unit, aggregate them or explicitly model measurement-level variation.

The statistical unit of inference should normally be the biological experimental unit.

---

# 4. Data Integrity Rules

The source Excel file is:

datos_crudos/dca/datos-proyectos tomillo-fusarium.xlsx

Each experimental design lives in its own folder (dca/, bdca/, factorial/),
with raw data under datos_crudos/<design>/ and outputs under
<design>/resultados/ (database, tablas, figuras, reportes, excel).

Preserve the original file.

Never overwrite raw data.

Never silently modify raw observations.

All transformations must be reproducible through code.

Every derived variable must document:

* Source variables
* Formula
* Reason for transformation

---

# 5. Critical Data-Quality Rules

Before statistical analysis:

1. Inspect all sheets.
2. Identify duplicate observations.
3. Identify missing values.
4. Identify inconsistent labels.
5. Identify inconsistent units.
6. Identify ambiguous concentration mappings.
7. Identify technical vs biological replicates.
8. Identify possible pseudoreplication.
9. Identify outliers.
10. Identify impossible values.

Do not remove outliers automatically.

Every exclusion must be documented.

---

# 6. Statistical Philosophy

The workflow must follow:

Data audit
→ Data cleaning
→ Experimental-unit definition
→ Exploratory analysis
→ Assumption diagnostics
→ Model selection
→ Primary inferential analysis
→ Post-hoc comparisons
→ Effect sizes
→ Confidence intervals
→ Dose-response analysis
→ Susceptibility analysis
→ Multivariate analysis
→ Biological interpretation

Do not rely only on p-values.

Report:

* Effect size
* Estimate
* 95% confidence interval
* Adjusted p-value where appropriate
* Biological interpretation

---

# 7. Model Selection

Do not assume ANOVA is always appropriate.

For each response variable:

1. Inspect distribution.
2. Inspect residuals.
3. Test or visually evaluate normality.
4. Evaluate homoscedasticity.
5. Evaluate independence.
6. Consider transformations.
7. Consider GLM/GLMM.
8. Consider mixed models if blocks or batches exist.

Model selection must be justified.

---

# 8. Inhibition Percentage

For inhibition percentage:

* Verify that the calculation is correct.
* Check whether percentages are bounded between 0 and 100.
* Inspect distribution.
* Evaluate whether transformation is appropriate.
* Consider beta regression only when assumptions and data structure support it.
* Avoid blindly applying arcsine square-root transformation.
* Consider linear mixed models if residual assumptions are acceptable.

If the experiment contains zero concentration controls, treat them according to the biological design and clearly document the decision.

---

# 9. Conidia Data

If conidia are counts:

Evaluate:

* Mean
* Variance
* Dispersion
* Zero inflation
* Distribution

Compare Poisson and negative-binomial approaches when justified.

If conidia values are continuous estimates rather than raw counts, document this and choose an appropriate model.

---

# 10. Dose-Response

Use nonlinear dose-response models when enough concentration levels exist.

Prefer a 4-parameter logistic model when supported by the data.

Estimate:

* EC50
* EC90

Provide confidence intervals where possible.

Do not calculate EC50 for groups with insufficient dose-response information.

Do not extrapolate beyond the observed concentration range without explicitly warning the user.

---

# 11. Fusarium Susceptibility

Build isolate-level susceptibility profiles using available metrics such as:

* EC50
* EC90
* Maximum inhibition
* Inhibition at defined concentrations
* Conidia reduction

Standardize variables before PCA/clustering when scales differ.

Use:

* PCA
* Hierarchical clustering
* Heatmaps

Evaluate cluster number using objective criteria such as silhouette score.

Do not overinterpret clusters when sample size is small.

---

# 12. Multiple Comparisons

Post-hoc testing must depend on the primary model.

Possible methods:

* Tukey HSD
* Estimated marginal means
* Dunnett comparisons against control
* Multiplicity-adjusted pairwise comparisons

Do not perform arbitrary pairwise tests without considering the factorial design.

When interactions are significant, prioritize simple effects and interaction contrasts over marginal main effects.

---

# 13. Reproducibility

The pipeline must be reproducible.

Use:

* Python
* pandas
* numpy
* scipy
* statsmodels
* pingouin
* scikit-learn
* matplotlib
* seaborn
* openpyxl

Additional packages may be used when justified.

Fix random seeds for stochastic procedures.

Save all final tables and figures.

---

# 14. Notebook Design

The final notebook should be educational and publication-oriented.

Every major section should contain:

1. Markdown explanation
2. Statistical rationale
3. Code
4. Output
5. Interpretation

Do not generate a notebook that only executes code.

The notebook must explain what each analysis means biologically.

---

# 15. Required Final Outputs

The final workflow should produce:

## Tables

* Dataset audit
* Missing data summary
* Experimental design summary
* Extraction yield statistics
* ANOVA/mixed model tables
* Effect sizes
* Post-hoc comparisons
* Inhibition summaries
* Conidia summaries
* EC50/EC90
* Isolate susceptibility ranking
* PCA scores
* Cluster assignments

## Figures

* Experimental design overview
* Distribution plots
* Boxplots
* Residual diagnostics
* Interaction plots
* Dose-response curves
* EC50 comparison
* Heatmap
* PCA
* Dendrogram

---

# 16. Scientific Integrity

Never fabricate data.

Never silently change raw values.

Never call an isolate resistant without a validated criterion.

Never treat technical replicates as independent biological replicates.

Never choose a statistical test only because it produces significance.

Always distinguish:

* Statistical significance
* Effect size
* Biological significance

If the data structure is ambiguous, stop and report the ambiguity before performing definitive inference.
