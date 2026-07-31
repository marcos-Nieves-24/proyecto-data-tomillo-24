# RCBD Yield Analysis Specification

## Purpose

Implements the statistical analysis for the DBCA (RCBD) Jenkyn mildew control
trial with a single response (yield): data audit, assumption checking,
classical block ANOVA plus a linear mixed model with random block (ICC),
Tukey HSD post-hoc with explicit pair-vs-reference (R) interpretation, and
raw means with 95% CI. Additivity is documented as untestable (one
observation per cell).

## Requirements

### Requirement: RCBD data audit

The pipeline MUST load the bdca raw CSV (36 rows: `plot, trt, block, yield`)
and MUST validate balance: one observation per `trt × block` cell, 9 per
treatment, 4 per block, no missing values or duplicates. Detected anomalies
MUST be reported, never imputed.

#### Scenario: Balanced dataset loads cleanly

- GIVEN the bdca CSV with 36 balanced observations
- WHEN the audit runs
- THEN a table reports 36 rows, no NAs, no duplicates, 9 per trt (R/T0/T1/T2) and 4 per block (B1–B9)

#### Scenario: Anomaly detection without imputation

- GIVEN an audit finding of a missing or duplicate value
- WHEN the audit completes
- THEN the anomaly is reported in the audit table
- AND the pipeline does not impute or drop observations

### Requirement: Assumption checking with documented justification

The pipeline MUST check normality of residuals (Shapiro-Wilk),
homoscedasticity and independence for the block model, and MUST document the
results and their consequences before selecting the inferential route, per
AGENTS.md model-selection rules.

#### Scenario: Assumptions satisfied

- GIVEN residuals pass normality and homoscedasticity checks
- WHEN the assumptions phase completes
- THEN the classical ANOVA F-test is reported as valid and justified

#### Scenario: Assumption violation

- GIVEN one or more assumptions fail
- WHEN the assumptions phase completes
- THEN the violation is reported together with the chosen compensating route
- AND the justification explains why that route is appropriate

### Requirement: Block model — classical ANOVA and LMM with ICC

The pipeline MUST fit the classical block ANOVA `yield ~ C(trt) + C(block)`
as an educational complement and MUST fit the LMM `yield ~ C(trt)` with
`groups=block` (REML) reporting the ICC as the primary block-aware model.

#### Scenario: Both models fit on balanced data

- GIVEN the audited balanced bdca dataset
- WHEN the models phase runs
- THEN the ANOVA table includes trt and block terms
- AND the LMM reports fixed effects, block variance, residual variance and ICC

#### Scenario: Additivity limitation documented

- GIVEN one observation per `trt × block` cell
- WHEN the report is generated
- THEN a limitation note states that additivity is not testable (no interaction term estimable)

### Requirement: Post-hoc Tukey HSD with reference interpretation

The pipeline MUST run Tukey HSD for all treatment pairs and MUST explicitly
interpret every pair against the reference treatment R. Dunnett MUST NOT be
used (no new dependencies).

#### Scenario: All pairs interpreted against the reference

- GIVEN a significant treatment effect in the model
- WHEN the post-hoc phase runs
- THEN pairwise comparisons cover all 6 treatment pairs
- AND each pair T0/T1/T2 vs R is interpreted explicitly with adjusted p-values

### Requirement: Raw means and 95% confidence intervals

The pipeline MUST report raw yield means and 95% CI per treatment, computed
from observed data only.

#### Scenario: Means and intervals reported

- GIVEN the audited dataset
- WHEN the summary tables are produced
- THEN each treatment row includes mean, standard error and 95% CI of raw yield
