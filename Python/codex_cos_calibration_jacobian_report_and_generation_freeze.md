# Codex Task: Final COS Calibration, Implied-State Jacobian, Report Documentation, and Sample-Generation Freeze

## Repository and starting point

Repository: `andlalu/P-ONE`  
Branch: `master`  
Current reviewed refactor commit: `2834ac02d9840d2ecb31bbc0e863a6a38ae39aa7`

Primary code scope:

```text
Python/Models/Heston
Python/OptionPricing
Python/OptionData
Python/Estimation/ISCGMM
Python/DGPSimulation
Python/Scripts
Python/Scripts/configs
```

Report:

```text
Reports/18062026/draft_proposalJun26.tex
```

Read this specification in full before changing code. This task is the final numerical and configuration gate before the production Monte Carlo samples are generated. The task must leave the repository with a calibrated, documented and reproducible production-generation configuration.

Do not generate the full Monte Carlo sample set in this task. Generate exactly one end-to-end sample only after the calibration and code changes have passed their tests.

---

# 1. Objectives

The task has six linked objectives:

1. correct the remaining issues identified after the IS-CGMM refactor;
2. calibrate the maturity-specific fixed effective COS widths and COS term counts;
3. implement and validate the implied-volatility Jacobian with respect to the implied variance state;
4. benchmark the available scalar state-implying methods and draw a documented conclusion;
5. update the report concisely with the fixed-basis method and Jacobian derivation;
6. freeze and validate the production sample-generation configuration.

The final output must include an explicit recommendation and the selected numerical values. Do not stop at producing benchmark tables without drawing a conclusion.

---

# 2. Correct remaining refactor issues

## 2.1 Distinguish basis compatibility from identical `n_cos`

The mathematical compatibility requirement is:

```text
same basis convention
same maturity grid
same maturity-specific effective widths
compatible Heston COS implementation family/version
```

Generation and estimation are allowed to use different COS term counts.

This is intentional:

- panel generation may use a more accurate, larger `n_cos`;
- estimation may use a smaller `n_cos` if the calibration shows that it remains sufficiently accurate;
- this avoids requiring the estimator to reproduce exactly the same numerical approximation used to generate the panel.

Refactor the current metadata check so that it does not always require equal `n_cos`.

Use explicit metadata names such as:

```text
generation_n_cos
estimation_n_cos
basis_convention
pricing_implementation
maturities
effective_widths
maturity_tolerance
```

A method such as the following is acceptable:

```python
assert_compatible_metadata(
    metadata,
    *,
    require_same_n_cos: bool = False,
)
```

Exact `n_cos` equality may remain available for regression tests or exact-reproduction runs.

Panel metadata must record the generation value. Estimation output must record both:

```text
panel generation n_cos
estimation n_cos
```

Do not overwrite the generation value with the estimation value.

## 2.2 Validate an explicit `dt`

The criterion currently verifies equal panel spacing but can accept a supplied `dt` that differs from the observed interval.

Change the logic so that:

1. the panel interval is inferred and checked for constancy;
2. when `config.dt` is omitted, the inferred interval is used;
3. when `config.dt` is supplied, it must equal the inferred interval within `spacing_tolerance`;
4. a descriptive error is raised on a mismatch.

Do not silently evaluate a weekly panel using another transition horizon.

## 2.3 Strengthen canonical panel validation

Extend `OptionPanelDate` and `OptionPanel` validation to cover:

- finite observed IVs;
- positive observed IVs where required by the estimation convention;
- finite observed prices when present;
- finite rates and dividend yields;
- option types restricted to `call` and `put`;
- `spot` consistent with `exp(log_spot)` within a documented tolerance;
- finite true variance and non-negative true variance when present;
- strictly increasing panel times;
- unique and consistently ordered date indices;
- no empty date slices;
- aligned optional vectors.

Keep the panel object general enough for the future particle-filter estimator.

## 2.4 Repository hygiene

Remove tracked:

```text
__pycache__/
*.pyc
```

files from version control. The existing generic `.gitignore` entries should remain. Remove redundant file-specific cache entries where this simplifies the file without changing behaviour.

---

# 3. Establish one shared production COS-basis configuration

## Problem

The generation configuration currently embeds COS settings while the estimator uses a separate basis JSON. This permits the two configurations to drift.

## Required change

Create one authoritative production basis file, for example:

```text
Python/Scripts/configs/heston_cos_basis_production.json
```

It should contain:

```json
{
  "basis_convention": "fixed_effective_width",
  "pricing_implementation": "heston_cos_fixed_basis_v1",
  "maturities": [...],
  "effective_widths": [...],
  "generation_n_cos": ...,
  "estimation_n_cos": ...,
  "maturity_tolerance": ...
}
```

The exact final numerical values are to be filled from the calibration exercise in this task.

Both:

```text
sample generation
IS-CGMM estimation
```

must load this same file or the same validated shared object.

Do not duplicate the calibrated widths in several JSON files.

The POC-only file:

```text
Python/Scripts/configs/is_cgmm_poc_test_basis.json
```

may remain as an explicitly named test fixture, but it must not be mistaken for the production configuration.

---

# 4. Fixed effective COS-width calibration exercise

## 4.1 Purpose

Select:

- one fixed effective COS half-width for each production maturity;
- a production `generation_n_cos`;
- a production `estimation_n_cos`.

The selected basis must be sufficiently accurate over the state and parameter region the estimator is allowed to visit. It is not enough to calibrate only at the true DGP parameter vector.

The result of this exercise is a methodological conclusion, not merely a diagnostic output.

## 4.2 Production contracts

Use the exact production contract grid:

```text
maturities: 1/12, 1/4, 1/2 years
log-forward moneyness: -0.15, -0.075, 0, 0.075, 0.15
OTM option convention used by the panel generator
r = 0.02
q = 0.00
```

## 4.3 State and parameter domain

Cover at least:

```text
initial variance over the intended state-search interval
the true DGP parameter vector
the POC starting vector
moderate perturbations around the true vector
risk-neutral parameter combinations induced by the intended optimiser bounds
```

Do not use a naive full Cartesian product of extreme parameter bounds if it creates economically meaningless combinations. Instead:

1. document the admissible calibration domain;
2. use a reproducible space-filling or structured design over that domain;
3. enforce the model restrictions;
4. report excluded combinations and the reason for exclusion;
5. assess whether the current production optimiser bounds should be tightened.

The calibration must explicitly check the implied risk-neutral parameters:

```text
kappa_q
vbar_q
sigma_v
rho
```

If the current optimiser bounds permit regions in which no economical fixed basis can achieve the target accuracy at reasonable cost, recommend and implement tighter production bounds with a concise numerical justification. Do not widen the bounds.

## 4.4 Independent high-accuracy reference

Construct a trusted high-accuracy reference.

Preferred approach:

1. implement or use an independent numerical Heston Fourier-inversion price for a representative validation subset;
2. use it to validate a very high-resolution COS reference;
3. use the validated high-resolution COS setup for the full calibration grid.

At minimum, the reference must not simply be the same candidate configuration being tested.

The existing variance-scaled reference path may be used as an additional cross-check, but it must not be treated as unquestioned ground truth.

Document:

```text
reference method
reference numerical tolerances
reference integration or COS settings
cross-check results
```

## 4.5 Candidate search

Evaluate a documented grid of:

```text
maturity-specific effective widths
n_cos values
```

A reasonable starting search set is:

```text
effective widths: 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00
n_cos: 64, 128, 256, 512, 1024
```

Refine locally around the best candidates when required.

Because the optimal width may differ by maturity, select widths separately by maturity while keeping the implementation and evaluation criteria common.

## 4.6 Error measures

Report at least:

```text
absolute price error
relative price error where meaningful
absolute Black implied-volatility error
vega-scaled price error
failure/invalid-IV count
runtime
```

Report distributional summaries and the maximum over the production calibration domain.

Use the observation-error scale as context. The baseline quote noise is `0.0015` in implied-volatility units. Numerical pricing/inversion error should be materially smaller.

Use the following initial acceptance targets unless the results demonstrate that a slightly different threshold is more defensible:

```text
generation:
    maximum absolute IV error <= 2e-6

estimation:
    maximum absolute IV error <= 1e-5

both:
    zero invalid prices or IV inversions over the accepted production domain
```

These correspond to numerical errors far below the 0.0015 baseline observation-error scale.

Do not hide isolated failures behind average-error measures.

## 4.7 Required conclusion

Select the lowest-cost specification satisfying the tolerances.

The conclusion must state:

```text
effective width for each maturity
generation_n_cos
estimation_n_cos
accepted state/parameter domain
maximum observed error
99th percentile error
runtime comparison
whether optimiser bounds were changed
```

Write machine-readable outputs under a stable diagnostics directory, for example:

```text
outputs/calibration/fixed_cos/
```

Include:

```text
calibration_summary.csv
calibration_selected.json
reference_crosscheck.csv
```

Optional plots should be concise and directly useful for reviewing the choice.

Update `heston_cos_basis_production.json` and the production generation/estimation profiles with the selected values.

---

# 5. Implied-state Jacobian

## 5.1 Object to differentiate

For date `i`, define the model-IV vector:

```text
sigma_model_i(v; theta_Q)
```

and its Jacobian with respect to the scalar initial variance state:

```text
J_i(v; theta_Q)
    = partial sigma_model_i(v; theta_Q) / partial v
```

This is a `q_i x 1` vector.

The date-level residual and loss are:

```text
r_i(v) = sigma_obs_i - sigma_model_i(v)
L_i(v) = r_i(v)' r_i(v)
```

so that:

```text
dL_i/dv = -2 J_i(v)' r_i(v)
```

## 5.2 Numerical Jacobian

Implement an independent finite-difference Jacobian for validation and fallback.

Requirements:

- central difference in the interior;
- one-sided difference near the variance bounds;
- configurable relative and absolute step sizes;
- step-size clipping to the admissible variance interval;
- clear handling of failed IV inversions;
- no silent replacement of non-finite derivatives with zero.

Use an explicit name such as:

```python
finite_difference_iv_variance_jacobian(...)
```

This function is a numerical reference and diagnostic, not automatically the preferred production implementation.

## 5.3 Semi-analytical Jacobian under the fixed basis

Implement the semi-analytical initial-variance sensitivity.

Under a fixed COS basis, the affine term is:

```text
exp(A_k + B_k v)
```

and therefore:

```text
partial/partial v exp(A_k + B_k v)
    = B_k exp(A_k + B_k v)
```

Use this to compute the Heston price derivative with respect to the initial variance while reusing the prepared basis and affine coefficients.

For an interior Black implied-volatility solution, use the implicit relation:

```text
partial sigma_model_j / partial v
    =
    (partial HestonPrice_j / partial v)
    / BlackVega_j
```

Requirements:

- return price, IV, price derivative and IV derivative together where this avoids duplicated work;
- use the already prepared fixed basis;
- define a minimum acceptable Black vega;
- return a clear diagnostic or controlled failure when the vega is too small;
- account explicitly for clipping or bound behaviour;
- do not apply this formula to the variance-scaled reference path.

Use a descriptive API such as:

```python
price_iv_and_initial_variance_jacobian_fixed_basis(...)
```

## 5.4 Validation

Validate the semi-analytical Jacobian against the finite-difference Jacobian over:

```text
all three maturities
the production moneyness grid
the accepted state range
the accepted risk-neutral parameter calibration domain
```

Report:

```text
absolute Jacobian error
relative Jacobian error where meaningful
low-vega cases
boundary cases
runtime
```

Use several finite-difference step sizes to demonstrate that the numerical comparison is in its stable region.

Add regression tests with explicit tolerances.

---

# 6. State-implying solver benchmark and decision

## 6.1 Methods

Benchmark at least:

1. the current bounded golden-section search;
2. SciPy bounded Brent minimisation;
3. a bounded derivative-aware method using the residual vector and Jacobian, preferably `scipy.optimize.least_squares`;
4. the derivative-aware method with the finite-difference Jacobian;
5. the derivative-aware method with the semi-analytical Jacobian.

The old MATLAB implementation used a bounded nonlinear least-squares formulation. This is relevant context, but do not assume it is automatically best in the present implementation.

## 6.2 Benchmark panels

Use:

```text
clean panels
low-iid panels
spatially correlated panels
persistent-factor panels
```

from one or more deliberately small test samples. Do not generate the full production sample set.

Include dates that are:

```text
typical
high variance
low variance
close to a state-search boundary
affected by price capping or low vega where available
```

## 6.3 Measures

Compare:

```text
implied variance estimate
date-level loss
difference from a high-accuracy scalar reference
function evaluations
Jacobian evaluations
pricing calls
runtime
failure rate
boundary-hit behaviour
sensitivity to warm starts
```

## 6.4 Required conclusion

Codex must draw and document a conclusion.

Select the production state-implying method based on:

```text
accuracy
robustness
runtime
clarity
```

Do not select a derivative method merely because it is faster on clean panels.

If the semi-analytical Jacobian method is selected:

- make it the explicit production default;
- retain golden section as a reference/fallback only when justified;
- record fallback use in diagnostics.

If golden section or bounded Brent remains the production default:

- retain the Jacobian implementation and tests;
- explain numerically why the derivative method was not selected.

Update the production configuration with the selected solver and solver-specific tolerances.

---

# 7. Configuration freeze and simplification

## 7.1 Rename the production generation configuration

The current file:

```text
clean_generation_run_001.json
```

contains both clean and noisy panel specifications. Rename it to something that reflects its actual scope, for example:

```text
generation_run_001.json
```

Update all call sites and documentation.

## 7.2 Do not silently create noise defaults

At present an omitted noise block can be interpreted through default-noise construction.

For production generation:

```text
noise omitted or null => no noisy panels
explicit noise block => exactly the configured noisy panels
```

Do not silently introduce default noisy panels.

Tests may use explicit convenience constructors, but dissertation production runs must be fully explicit.

## 7.3 One source of truth for COS settings

The production generation config should reference the shared production COS-basis file rather than repeat widths and term counts.

The estimator profile should do the same.

## 7.4 Explicit output format

Do not select CSV or Parquet merely according to installed optional dependencies for production runs.

Add an explicit setting:

```text
panel_format: parquet
```

Production generation must fail clearly when the requested format is unavailable.

CSV fallback may remain for unit tests and small debugging exercises, but it must be explicit.

## 7.5 Do not override `return_daily` silently

The loader currently forces `return_daily=True`.

Either:

- respect the value in the configuration; or
- remove it from user configuration and set it once in a clearly documented production profile.

Do not accept one value in JSON and silently replace it.

Document whether daily arrays include burn-in observations.

## 7.6 Validate or remove apparent choices

Review:

```text
option_type_rule
pricing_method
iv_method
```

If only one value is supported:

- validate that exact value; or
- replace the apparent free choice with a named implementation constant.

Do not retain configuration fields that suggest unsupported alternatives.

## 7.7 Resolve paths reproducibly

Define how relative paths are interpreted.

Prefer one of:

```text
relative to repository root
relative to the configuration file
```

Do not depend on the process working directory.

Record the resolved absolute output root in the run metadata while keeping portable relative paths where practical.

## 7.8 Record the actual execution environment

The resolved run metadata should contain:

```text
Git commit SHA
configuration hash
Python version
NumPy/SciPy/pandas/pyarrow versions
requested worker count
resolved worker count
thread environment
panel format
basis configuration identifier/hash
```

`workers: null` must not remain the only record when a concrete worker count was used.

## 7.9 Simplify noisy-result manifest generation

Review the current flow in which noisy panels are generated within `generate_one_sample`, while a later helper may rescan or call generation again with `skip_existing=True` to construct a manifest.

Prefer returning the noisy-panel results as part of the per-sample result or writing the noisy manifest directly from the original generation results.

Do not rerun generation logic merely to discover outputs that were already created.

## 7.10 Atomic and resumable output

Before large-scale generation:

- write output files atomically where practical;
- write one per-sample status/result record;
- do not mark a sample complete until the path, clean panel, metadata and requested noisy panels are all complete;
- make `skip_existing` verify completeness rather than only file existence;
- preserve errors for failed samples.

---

# 8. One explicit generation CLI

Provide one clear command-line entry point, for example:

```text
Python/Scripts/generate_heston_samples.py
```

It should accept:

```text
--config
--sample-start
--sample-end
--workers
--skip-existing
--paths-only
--panels-only
--log-level
```

The semantics of `sample-end` must be documented as exclusive.

The script should:

1. load and validate the complete configuration;
2. resolve the shared COS-basis file;
3. resolve the output root;
4. run the requested samples;
5. write the resolved run configuration;
6. write clean and noisy manifests;
7. return a non-zero exit status when any requested sample fails.

Do not require users to import `Scripts.generation` manually to launch a production run.

---

# 9. One-sample final dry run

After all calibration, tests and configuration changes have completed, run exactly:

```text
sample 000
workers = 1
```

using the frozen production generation configuration.

Generate:

```text
one Heston state path
one clean option panel
one low-iid observed panel
one spatially correlated observed panel
one persistent-factor observed panel
one persistent-factor path file
all metadata sidecars
all manifests
```

Validate:

```text
array lengths
finite and non-negative variances
panel row counts
contract grid
option-type convention
price and IV finiteness
no-arbitrage price bounds
COS-basis metadata
generation/estimation basis compatibility
noise seeds
factor-file presence
capping counts
resolved configuration hash
```

Run one short clean implied-state inversion against this sample using the production estimation basis. This is not a full six-parameter C-GMM estimation run.

Do not generate samples 001--099 in this task.

---

# 10. Report updates

Report file:

```text
Reports/18062026/draft_proposalJun26.tex
```

## 10.1 Critical writing-style requirement

Any added text must match:

```text
the language style
sentence construction style
terminology
English level
degree of formality
paragraph density
mathematical notation
```

of the existing report.

This is a critical acceptance condition.

Before writing, inspect the surrounding text in:

```text
Section: Implied-State C-GMM Estimation
Appendix: Implementation of the IS C-GMM Estimator
```

Do not rewrite surrounding paragraphs merely to make the new text sound more polished. Do not introduce a noticeably different academic voice. Do not use generic AI-style headings, inflated claims, excessive signposting or long explanatory lists.

The additions must be concise and clearly reviewable.

## 10.2 Main-text addition

Add at most one short paragraph to the main implied-state section explaining that:

- option prices are evaluated using maturity-specific fixed effective COS truncation intervals;
- the same width convention is used in panel generation and state inversion;
- generation and estimation may use different numbers of COS terms after numerical validation;
- implementation and calibration details are deferred to the appendix.

Do not insert detailed numerical derivations in the main text.

## 10.3 Appendix fixed-basis subsection

Add a concise subsection under Appendix `\ref{app:is_cgmm_implementation}` covering:

1. the fixed interval for maturity `\mathcal T`;
2. the corresponding COS frequency grid;
3. the fact that the basis does not depend on candidate variance;
4. why this allows basis and affine-coefficient reuse during state inversion;
5. the selected widths and generation/estimation term counts;
6. a compact table of calibration errors and runtime;
7. a brief explanation of the numerical selection rule.

Use the report's existing notation. Do not reintroduce lower-case option maturity `\tau`; option maturity remains `\mathcal T`.

## 10.4 Jacobian derivation

Add a concise derivation of the implied-state Jacobian.

At minimum define:

```text
J_i(v; theta_Q)
    = partial sigma_i_model(v; theta_Q) / partial v
```

Show that under the fixed affine COS basis:

```text
partial exp(A_k + B_k v) / partial v
    = B_k exp(A_k + B_k v)
```

Then derive:

```text
partial sigma_model_j / partial v
    =
    (partial P_model_j / partial v)
    / Black-Scholes vega_j
```

for interior implied-volatility solutions.

For the date-level least-squares loss show:

```text
partial L_i(v) / partial v
    =
    -2 J_i(v)' [
        sigma_i_obs - sigma_i_model(v)
    ]
```

Also state briefly that a finite-difference Jacobian is retained for numerical validation.

Discuss low-vega and boundary cases in one concise paragraph.

Do not derive the gradient of the outer C-GMM criterion in this task.

## 10.5 State-solver description

Update the existing numerical-minimisation discussion to reflect the benchmark conclusion and selected production state solver.

State:

```text
solver selected
bounds and warm-start principle
fallback, if any
diagnostics recorded
```

Do not overstate efficiency based on a single sample.

## 10.6 Report reviewability

Provide:

```text
a unified diff for the report
a short list of inserted paragraphs/subsections
a short list of equations and labels added
```

Avoid unrelated formatting or wording changes.

Compile the LaTeX report and resolve all introduced errors, labels and references.

Do not add new citations unless they are already present in the bibliography and directly necessary.

---

# 11. Tests

Retain all existing behavioural tests unless the requested change explicitly makes one obsolete.

Add or update tests for:

1. basis compatibility with different generation and estimation `n_cos`;
2. strict compatibility of convention, maturities and widths;
3. explicit-`dt` agreement and mismatch rejection;
4. strengthened panel validation;
5. fixed-basis production metadata round-trip;
6. numerical Jacobian interior accuracy;
7. numerical Jacobian one-sided boundary behaviour;
8. semi-analytical Jacobian against finite differences;
9. low-vega handling;
10. price/IV derivative reuse of the prepared basis;
11. state-solver benchmark harness determinism;
12. selected production state solver;
13. explicit noise-null semantics;
14. explicit panel-format semantics;
15. shared production basis config loading;
16. resolved worker-count metadata;
17. atomic/completeness-aware `skip_existing`;
18. one-sample end-to-end generation;
19. clean and all three noisy-panel outputs;
20. report compilation where the environment supports it.

Run the complete Python test suite.

Do not weaken existing numerical tolerances solely to obtain a passing result.

---

# 12. Required deliverables

Provide:

```text
code changes
updated production configuration files
shared calibrated COS-basis JSON
calibration scripts
calibration CSV/JSON outputs
Jacobian validation outputs
state-solver benchmark outputs
one-sample dry-run outputs
updated report
report unified diff
complete test results
```

In the completion summary state explicitly:

```text
selected widths
generation_n_cos
estimation_n_cos
selected state solver
Jacobian implementation used in production
accepted calibration domain
maximum pricing/IV error
one-sample validation result
remaining risks before launching samples 000--099
```

---

# 13. Non-goals

Do not implement:

- the gradient/Jacobian of the outer C-GMM criterion;
- second-step C-GMM;
- the particle filter;
- AWS orchestration;
- the full 100-sample generation run;
- the full clean or noisy IS-CGMM estimation batch;
- new observation-error designs;
- broad report rewriting.

---

# 14. Acceptance criteria

The task is complete only when:

1. generation and estimation share one fixed-width basis specification;
2. generation and estimation may use separately calibrated `n_cos` values;
3. the selected fixed widths and term counts are supported by a reproducible numerical study;
4. Codex has drawn and documented a clear calibration conclusion;
5. the numerical and semi-analytical state Jacobians agree within documented tolerances;
6. the state-solver benchmark has produced and implemented a clear conclusion;
7. production configuration contains no placeholder COS values;
8. no hidden noise defaults or output-format fallbacks remain in the production path;
9. one explicit generation CLI exists;
10. run metadata records the actual numerical and execution configuration;
11. one end-to-end sample and all requested noisy variants pass validation;
12. the report concisely documents the fixed basis, Jacobian derivation and selected state solver;
13. all report additions match the existing writing style and are clearly reviewable;
14. the full Python test suite passes;
15. no samples beyond sample 000 have been generated.
