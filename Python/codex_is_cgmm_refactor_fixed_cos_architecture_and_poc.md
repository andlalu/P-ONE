# Codex Task: Refactor the Heston IS-CGMM Pipeline and Introduce Fixed COS-Basis Architecture

## Repository and reference point

Repository: `andlalu/P-ONE`  
Reference commit: `449483940d02511cb7996b3f0117770a1feecae4`

Primary scope:

```text
Python/Estimation/ISCGMM
Python/DGPSimulation
Python/OptionPricing
```

This task should make the IS-CGMM implementation clean, minimal, reusable and ready for a local clean-panel proof of concept. It must also avoid creating method-specific duplicates that would later need to be repeated for the particle-filter estimator.

The current first-step C-GMM criterion is structurally good. Preserve its overall flow and mathematical behaviour unless a change is explicitly requested below.

---

## Guiding principles

1. Prefer one canonical representation of a model parameter set and one canonical representation of an option panel.
2. Do not create PF-specific or IS-specific copies of shared market-data structures.
3. Keep model mathematics in shared model modules, not inside pricing or estimation packages when the same mathematics is used in both.
4. Use descriptive, applied names. Avoid generic names such as `types`, `_ScalarMinResult`, `_ContractGroup` and `_DatePricingCache` when the object has a clear domain meaning.
5. Remove unused compatibility code after confirming it has no call sites.
6. Preserve existing tests except where an explicitly requested functional change makes a test obsolete.
7. Do not weaken numerical tolerances or delete tests merely to make the suite pass.
8. Run the complete Python test suite after the refactor, not only the IS-CGMM tests.
9. Introduce the fixed effective COS-basis architecture before any production panel samples are generated.
10. Preserve every requirement in this specification; do not silently narrow the scope when moving files or resolving dependencies.

---

# 1. Shared Heston parameter model

## Problem

Heston parameters are currently represented separately in:

```text
Python/DGPSimulation/types.py
Python/OptionPricing/types.py
Python/Estimation/ISCGMM/types.py
```

The estimation package defines another joint parameter class and then converts it into the physical and risk-neutral parameter classes. This duplicates model concepts and makes ownership unclear.

## Required change

Create a shared package:

```text
Python/Models/
    __init__.py
    Heston/
        __init__.py
        parameters.py
        affine_transform.py
```

In `parameters.py`, define clear shared parameter containers:

```python
HestonParameters
HestonPhysicalParameters
HestonRiskNeutralParameters
```

Recommended responsibilities:

### `HestonParameters`

Canonical joint parameter vector used by estimation:

```text
eta
kappa
vbar
sigma_v
rho
eta_v
r
q
```

It should provide:

```python
validate()
kappa_q
vbar_q
to_physical()
to_risk_neutral()
```

### `HestonPhysicalParameters`

Parameters required by the physical-measure simulator and transition law.

### `HestonRiskNeutralParameters`

Parameters required by the option-pricing code.

Do not use inheritance merely to avoid repeating six field declarations. Explicit parameter projections are preferable because the physical and risk-neutral objects have different meanings.

## Migration

Update simulation, option pricing and estimation imports to use the shared classes. Remove method-specific parameter dataclasses once all call sites have migrated.

Do not retain duplicate aliases indefinitely. A very short compatibility layer is acceptable only if required to migrate the repository in the same change, but the final result should have one authoritative definition for each parameter representation.

---

# 2. Split `ISCGMM/types.py`

## Problem

`types.py` currently contains:

- parameter classes;
- parameter transformations;
- panel data structures;
- configuration objects;
- result objects.

This is too broad and obscures ownership.

## Required target structure

```text
Python/Estimation/ISCGMM/
    config.py
    parameter_transform.py
    results.py
    implied_state.py
    heston_transition_cf.py
    cgmm_criterion.py
    estimate.py
    tests/
```

### `config.py`

Move:

```text
ImpliedStateConfig
CcfQuadratureConfig
CgmmConfig
OptimizerConfig
LoggingConfig
```

### `parameter_transform.py`

Move or implement:

```python
to_free(...)
from_free(...)
free_parameter_bounds(...)
```

The transformation should continue to enforce:

```text
kappa > 0
vbar > 0
sigma_v > 0
kappa_q > 0
rho in (-1, 1)
```

### `results.py`

Move or implement:

```text
ImpliedStateResult
CriterionDiagnostics
FirstStepEstimate
OptimizerStageResult
```

Result classes should be serialisable to plain dictionaries without embedding live estimator objects.

Remove `types.py` after all imports have migrated.

---

# 3. One canonical option-panel representation

## Problem

The repository currently has a dense `OptionPanel` under option pricing and a separate `EstimationPanel`/`PanelDate` representation under IS-CGMM. The latter is useful, but the ownership is wrong and the duplication will become worse when the particle filter is added.

## Required change

Create:

```text
Python/OptionData/
    __init__.py
    panel.py
    io.py
```

Define one canonical date-sliced representation:

```python
OptionPanelDate
OptionPanel
```

The canonical panel should be sufficiently general for:

- clean Monte Carlo panels;
- contaminated panels;
- IS-CGMM;
- the future particle-filter estimator.

Suggested fields:

```text
date_index
time
spot
log_spot
strikes
maturities
option_types
observed_iv
observed_price
clean_iv
rates
dividend_yields
true_variance
metadata
```

Optional fields should remain optional rather than generating separate subclasses for clean, noisy, IS or PF use cases.

## Dense pricing object

If the existing dense option-pricing `OptionPanel` remains useful for vectorised pricing, rename it to something explicit such as:

```python
OptionPriceCube
```

It should be treated as an internal pricing output rather than the canonical market/panel data structure.

## Loading

Move CSV/Parquet loading from `Estimation/ISCGMM/panel.py` into `OptionData/io.py`.

The estimator should receive an `OptionData.OptionPanel` directly.

## Remove unused adapter

`panel_from_option_panel` appears to be a compatibility adapter for an older dense NPZ schema.

Before deleting it:

1. perform a repository-wide call-site search;
2. confirm it is unused by scripts and tests;
3. remove it if unused;
4. remove any test that exists solely for this unused adapter.

Do not leave unused compatibility functions in the repository.

---

# 4. Refactor `implied_state.py` for applied readability

## Preserve the mathematics

The implied state remains the scalar variance value that minimises the cross-sectional sum of squared IV errors at each date.

The path should continue to:

1. start the first date from the long-run variance;
2. warm-start each later date from the previously implied variance;
3. search locally first;
4. fall back to the full admissible variance interval if the local solution hits an artificial local boundary;
5. return path-level diagnostics.

## Rename internal objects and functions

Use descriptive applied names. The following mapping is preferred:

| Current | Preferred |
|---|---|
| `_ScalarMinResult` | `_ImpliedVarianceSolveResult` |
| `_CacheStats` | `_PricingCacheDiagnostics` |
| `_ContractGroup` | `_MaturityPricingBatch` |
| `_DatePricingCache` | `_DateImpliedStatePricingContext` |
| `_bounded_minimize_scalar` | `_bounded_golden_section_search` |
| `_make_date_cache` | `_prepare_date_pricing_context` |
| `_coefficients_for_group` | `_get_cos_coefficients` |
| `_model_iv_for_date` | `_compute_model_implied_volatilities` |
| `_date_objective` | `_implied_variance_fit_loss` |
| `_warm_bounds` | `_local_variance_search_bounds` |
| `_minimize_date_variance` | `_imply_variance_for_date` |

Equivalent names are acceptable if they are equally clear.

## Explain the maturity grouping

Retain grouping by maturity because it allows:

- one affine coefficient solve per maturity;
- vectorised pricing across strikes;
- validation of rates and dividend yields within a maturity;
- reuse of COS coefficients.

The code and docstrings should make this explicit. Avoid generic “contract group” terminology.

## Improve file flow

Put the applied entry point and workflow near the top of the file, followed by the lower-level pricing and search helpers.

The file should read conceptually as:

```text
imply_heston_variance_path
    prepare date pricing contexts
    loop over dates
        determine local variance bounds
        imply variance for date
            compute model IVs
            evaluate IV fit loss
    return path and diagnostics
```

Keep only intermediary structures that provide a genuine computational or readability benefit.

---

# 5. State-implying solver: retain golden section, prepare for derivative work later

## Current decision

Keep bounded golden-section search for this refactor. Do not replace it with Brent or a gradient-based solver in the production path yet.

## Semi-analytical derivative feasibility

A useful derivative is feasible under a fixed COS basis.

For a date-level loss

```text
L(v) = sum_j [IV_obs_j - IV_model_j(v)]^2
```

the derivative is

```text
dL/dv = -2 sum_j residual_j * d IV_model_j / dv
```

For an interior Black implied-volatility solution,

```text
d IV_model_j / dv
    = (d HestonPrice_j / dv) / BlackVega_j
```

Under a fixed COS grid, Heston prices depend on the initial variance through the affine term

```text
exp(A + B v)
```

so the price derivative can be computed semi-analytically by multiplying by `B`.

However, the current pricing path may allow the effective COS truncation width to change with `v`. In that case the grid, coefficients and payoff expansion also change, so the above derivative is incomplete. Clipped prices, clipped IVs and very low Black vegas also introduce non-smooth or unstable points.

The refactor must therefore support the fixed effective COS-basis architecture described in the next section. This makes the later semi-analytical variance derivative mathematically well defined without requiring derivatives of the truncation rule itself.

## Requirement for this task

Do not add a new derivative-based production solver now. Instead:

1. keep the code structured so that a future method such as
   `compute_model_iv_and_variance_sensitivity(...)`
   can be added cleanly;
2. document in the relevant docstring that the semi-analytical derivative assumes a fixed effective COS basis;
3. factor the fixed-basis pricing kernel so that the affine coefficient `B` and the already-computed basis can be reused by a later sensitivity implementation;
4. do not make hidden assumptions about truncation-width derivatives;
5. leave a concise technical note in code comments or documentation, not a large TODO framework.

A later, separate numerical task should compare golden section, bounded Brent and a safeguarded derivative-aided method on clean and contaminated dates.

Do not implement the gradient of the outer C-GMM criterion in this task. The code should remain modular enough to add implied-state sensitivities, transition-CCF sensitivities and kernel sensitivities later, but that work is explicitly outside the present refactor.

---


# 6. Fixed effective COS-basis architecture before sample generation

## Decision and scope

No production clean or contaminated option-panel samples should be generated before the COS truncation convention is fixed.

The current variance-dependent effective truncation width complicates the implied-state derivative because a change in candidate variance can also change:

- the truncation interval;
- the Fourier frequencies;
- the payoff coefficients;
- the cached affine coefficients.

Introduce a fixed effective COS-basis architecture now so that clean-panel generation and implied-state estimation use the same pricing convention.

This refactor should implement the architecture and validation. It must **not** choose the final calibrated numerical widths. The actual width and `n_cos` calibration will be a separate numerical task after this refactor.

## Shared fixed-basis configuration

Create one clear shared configuration object for the fixed COS basis, for example:

```python
FixedCosBasisConfig
```

It should represent at least:

```text
maturity grid
effective truncation width for each maturity
number of COS terms
numerical matching tolerance for maturities
```

Prefer aligned arrays or explicit records over an unvalidated floating-point dictionary.

The configuration must validate:

- strictly positive effective widths;
- strictly positive maturities;
- unique maturity entries;
- a one-to-one match between maturities and widths;
- `n_cos > 1`;
- requested panel maturities can be matched unambiguously.

Do not use one generic field whose meaning changes between fixed and variance-dependent truncation.

## Same convention in generation and estimation

Both the option-panel generation code and `implied_state.py` must consume the same fixed-basis specification.

The required convention is:

- the effective truncation width is fixed for a given maturity;
- it does not change with the simulated or candidate variance;
- the COS grid and payoff expansion are therefore fixed for that maturity;
- affine coefficients can be precomputed and reused for all candidate variance values at a fixed parameter vector.

Generation and estimation may later use different `n_cos` values for a controlled accuracy study, but they must use the same fixed-width convention and explicitly recorded maturity-specific widths. Do not introduce an undocumented difference between the generation and estimation pricing approaches.

## Pricing API

Refactor the COS-pricing API so that fixed effective widths are supplied explicitly.

The fixed-basis path must not call an automatic function such as:

```python
effective_truncation_width(variance, maturity, ...)
```

during candidate-variance evaluation.

For a fixed parameter vector and maturity, prepare and reuse:

```text
COS frequency grid
payoff expansion terms
Heston affine A/B coefficients
```

Only the affine evaluation in the initial variance should change during the inner state search.

If the existing variance-dependent truncation path is retained for a high-accuracy reference or a future calibration comparison, it must:

- have an explicit descriptive name;
- not remain the silent default;
- have verified call sites and tests;
- not be used by the production panel-generation or implied-state paths.

Otherwise remove it after confirming it is unused.

## Sample and result metadata

Store the complete COS-basis specification in:

- generated clean-panel metadata;
- contaminated-panel metadata inherited from the clean panel;
- estimation result metadata;
- diagnostic output required to reproduce a job.

At minimum record:

```text
basis convention = fixed effective width
maturity-specific effective widths
n_cos
pricing implementation/version identifier
```

When an estimator loads a generated panel, validate the basis metadata against the estimation configuration. A clear error is preferable to silently estimating with a different truncation convention.

## Derivative readiness

With a fixed basis, the dependence of a COS price on the initial variance is through the affine term:

```text
exp(A + B v)
```

A later sensitivity implementation can therefore use:

```text
d/dv exp(A + B v) = B exp(A + B v)
```

The refactor should keep the fixed-basis pricing calculation factored so that:

- the basis is computed once;
- the affine `B` coefficients are available;
- price and variance sensitivity can later be evaluated together without rebuilding the grid.

Do not yet add a derivative-based state solver or criterion gradient.

## No arbitrary production widths in this task

Do not hard-code final production widths merely to complete the refactor.

Tests may use small explicit fixture widths.

The subsequent fixed-COS calibration task will:

1. define a high-accuracy pricing reference;
2. test candidate widths and `n_cos` values over the intended moneyness, maturity, variance and parameter ranges;
3. select production generation and estimation settings;
4. validate price and IV errors against tolerances below the observation-noise scale;
5. implement and finite-difference-check the initial-variance price/IV sensitivity;
6. benchmark golden section, bounded Brent and a safeguarded derivative-aided state solve;
7. only then generate the production clean and contaminated panels.

---

# 7. Shared affine/Riccati machinery

## Problem

The transition CCF correctly differs from the option-pricing CCF, but shared affine mathematics currently lives in the option-pricing package.

## Required change

Move the generic constant-Riccati solver into:

```text
Python/Models/Heston/affine_transform.py
```

Both:

```text
OptionPricing/heston_ccf_solver.py
Estimation/ISCGMM/heston_transition_cf.py
```

should import the shared Riccati implementation.

The option-pricing solver should remain a thin risk-neutral wrapper.  
The transition solver should remain a thin physical-transition wrapper with the non-zero terminal variance loading required by the joint return/variance transform.

Do not blindly replace both analytic implementations with one formula if this changes complex-log branch behaviour. Consolidate only the genuinely shared machinery, and preserve analytic-versus-RK and pricing regression tests.

---

# 8. Keep one CCF quadrature method

Remove Gauss-Legendre support.

Retain Gauss-Hermite only.

Replace the generic quadrature configuration with something explicit, for example:

```python
@dataclass(frozen=True)
class CcfQuadratureConfig:
    dimension: int = 2
    order: int = 3
    scale: float = 1.0
```

Remove:

```text
scheme
gauss_legendre branch
Gauss-Legendre-specific tests
```

Keep the current tensor-product construction and normalisation for Gaussian integration.

---

# 9. Require equally spaced estimation observations

The Monte Carlo estimation panel is weekly and equally spaced. The current median-interval fallback can silently accept irregular observations and then use one interval for every transition.

Replace `_panel_dt` with a helper such as:

```python
_infer_constant_transition_interval(...)
```

It should:

1. compute all consecutive time differences;
2. require them to be positive;
3. verify they are equal within a clear numerical tolerance;
4. return the common interval;
5. raise a descriptive error for irregular panels.

A user-supplied `dt` may remain as an explicit override, but the validation logic should remain clear.

Do not add transition-specific irregular-time CCF handling in this task.

---

# 10. Production first-step estimator

## Replace the coordinate-search fallback

Replace the current simple coordinate search in `estimate.py` with bounded derivative-free Powell optimisation using SciPy.

The criterion object should remain independent of the optimiser.

## Starting-value strategy

Do not perform a broad grid search or a global search before Powell.

Add a small deterministic candidate-start screen:

- one configured base start;
- at most two modest perturbations of the base start;
- perturbations should be approximately 10–15% in economically meaningful coordinates;
- evaluate the criterion once at each candidate;
- choose the best candidate and start Powell there.

The default number of screened starts must not exceed three.

The base start and bounds must be supplied through configuration. Do not hard-code the simulation truth inside the estimator implementation.

For the clean Monte Carlo POC, provide a configuration profile with a base start reasonably close to the DGP and moderate bounds. Suggested POC bounds in natural parameter space are:

```text
eta:       [0.0, 10.0]
kappa:     [2.0, 12.0]
vbar:      [0.005, 0.080]
sigma_v:   [0.15, 0.80]
rho:       [-0.90, -0.10]
kappa_q:   [0.50, 6.00]
```

These are POC defaults, not universal model-validity constraints. Keep them configurable.

Convert these bounds consistently into the free parameterisation used by Powell.

## Staged Powell

Implement two stages:

### Stage 1

- bounded Powell;
- moderate function-evaluation limit;
- loose tolerances;
- purpose: move reliably into the relevant basin.

### Stage 2

- restart from the Stage-1 solution;
- tighter tolerances;
- larger but still controlled evaluation limit.

Then perform one final criterion evaluation with full diagnostics.

## Failure handling

Expected numerical failures from:

- invalid transformed parameters;
- failed implied-state inversion;
- non-finite pricing;
- non-finite CCF evaluation;

should return a large finite penalty and increment a counter.

Unexpected programming errors must propagate. Do not catch `Exception` broadly.

## Estimation result

`FirstStepEstimate` should contain at least:

```text
candidate starts and their criteria
selected start
initial criterion
estimated parameters
free parameters
final criterion
success
status/message
iterations
function evaluations
penalty evaluations
stage results
total runtime
final criterion diagnostics
```

Provide a clean conversion to a JSON-serialisable dictionary.

---

# 11. Minimal but useful progress logging

Diagnostics stored in result objects are not sufficient for long cloud jobs. Add lightweight execution logging using Python's standard `logging` module.

Do not add an external logging dependency.

## Required INFO-level events

Log:

1. estimation job started;
2. panel loaded, including sample/scenario identifiers and panel dimensions;
3. candidate-start screening completed, including the three criterion values;
4. Powell Stage 1 started;
5. Powell Stage 1 completed;
6. Powell Stage 2 started;
7. Powell Stage 2 completed;
8. final diagnostic evaluation completed;
9. result written;
10. job completed or failed.

## Periodic progress

During Powell, emit one concise progress message:

- every configurable `N` criterion evaluations, default `N=10`;
- include evaluation count, current best criterion, elapsed time and penalty count.

Do not log full parameter vectors on every evaluation. It is acceptable to include the current best parameter vector only at stage boundaries.

## DEBUG level

DEBUG may include:

- selected transformed parameter vectors;
- implied-state summary statistics;
- boundary-hit and failure rates;
- component timing summaries.

Do not log:

- every date-level variance solve;
- every option price;
- every CCF node;
- full implied-state arrays;
- full objective history for all successful jobs.

Detailed arrays for a deliberately selected diagnostic job or a failed job may be written separately to NPZ/Parquet by the calling script. Do not make this the default for every successful batch job.

Logging must work naturally with stdout/CloudWatch when the estimator is later run on AWS.

---

# 12. Preserve the C-GMM criterion structure

The existing `cgmm_criterion.py` is well structured. Preserve the high-level sequence:

```text
imply states
construct transition states
evaluate CCF residuals
construct Gaussian instrument kernel
evaluate the quadratic form
return diagnostics
```

Only make the requested changes:

- new imports;
- shared parameter and panel types;
- Gauss-Hermite-only quadrature;
- constant-spacing validation;
- clearer naming where needed.

Do not redesign the criterion or add second-step C-GMM in this task.

---

# 13. Tests

## Preserve existing tests

All behavioural tests should remain unless a requested functional change makes one obsolete.

Specifically:

- remove tests solely for Gauss-Legendre;
- remove tests solely for the unused old-panel adapter if the adapter is deleted;
- update imports and names after file moves;
- do not remove or weaken other tests.

## Add tests for

1. shared Heston parameter projections;
2. physical and risk-neutral validation;
3. canonical panel loading from CSV/Parquet;
4. constant-spacing acceptance;
5. irregular-spacing rejection;
6. Gauss-Hermite quadrature node/weight shape and determinism;
7. Powell reproducibility;
8. Powell bound compliance;
9. no more than three starting candidates;
10. penalty handling for expected numerical failures;
11. periodic logging without excessive per-evaluation output;
12. result serialisation;
13. analytic transition CCF against RK4;
14. option-pricing regression after moving Riccati code;
15. clean implied-state recovery tests already present;
16. fixed COS-basis validation and maturity matching;
17. fixed-basis grid independence from candidate variance;
18. reuse of the same fixed-width convention in generation and estimation;
19. COS-basis metadata round-trip;
20. rejection of a panel/configuration COS-basis mismatch;
21. regression showing that the fixed-basis pricing path does not silently invoke variance-dependent truncation.

Run the complete Python test suite.

---

# 14. CLI / proof-of-concept readiness

Provide or update a command-line script that can:

```text
select an input panel
select clean/noisy IV column
truncate to a maximum number of dates
load the fixed COS-basis specification
load optimiser/config profile
run first-step IS-CGMM
write one JSON result
configure log level
```

The command should be suitable for:

```text
sample 000, 30 dates
sample 000, 100 dates
sample 000, full panel
```

Do not add AWS orchestration in this task.

---

# 15. Expected target structure

A reasonable final structure is:

```text
Python/
    Models/
        __init__.py
        Heston/
            __init__.py
            parameters.py
            affine_transform.py

    OptionData/
        __init__.py
        panel.py
        io.py

    DGPSimulation/
        __init__.py
        base.py
        heston_simulator.py
        variance_drawers.py
        types.py
        io.py
        tests/

    OptionPricing/
        ...

    Estimation/
        ISCGMM/
            __init__.py
            config.py
            parameter_transform.py
            results.py
            implied_state.py
            heston_transition_cf.py
            cgmm_criterion.py
            estimate.py
            tests/
```

`DGPSimulation/types.py` may continue to hold simulation-specific configuration and path-output classes, but it must no longer define a separate Heston parameter class after migration.

---

# 16. Non-goals

Do not implement in this task:

- second-step C-GMM;
- AWS/EC2/AWS Batch orchestration;
- particle-filter estimation;
- new DGP simulation logic;
- new observation-error designs;
- report plots;
- final numerical calibration of the maturity-specific COS widths;
- production clean or contaminated sample generation;
- production use of a derivative-based implied-state solver;
- implementation of the outer C-GMM criterion gradient/Jacobian;
- a global optimiser or a large initial grid search.

---

# 17. Acceptance criteria

The task is complete when:

1. one canonical Heston parameter model is used across simulation, pricing and estimation;
2. one canonical date-sliced option-panel model is available to both IS-CGMM and future PF work;
3. unused old-panel conversion code is removed;
4. `ISCGMM/types.py` is removed and responsibilities are split clearly;
5. `implied_state.py` reads in an applied, domain-oriented manner;
6. a shared fixed effective COS-basis configuration exists with maturity-specific widths;
7. clean-panel generation and implied-state estimation can consume the same fixed-width convention;
8. the production fixed-basis path does not change its grid with candidate variance;
9. COS-basis metadata is stored and panel/configuration mismatches are rejected;
10. the pricing structure is ready for a later semi-analytical initial-variance sensitivity without implementing the production derivative solver now;
11. shared Riccati machinery no longer lives only inside option pricing;
12. only Gauss-Hermite quadrature remains;
13. irregular observation spacing is rejected explicitly;
14. staged bounded Powell is the production first-step optimiser;
15. initial screening uses no more than three nearby deterministic starts;
16. moderate configurable parameter bounds are enforced;
17. concise progress logging is available during long runs;
18. one local clean-panel POC can be launched from the command line using an explicit fixed-basis test configuration;
19. all relevant existing tests remain and the complete test suite passes;
20. no unused duplicate classes, adapters or compatibility wrappers remain;
21. no arbitrary final production COS widths or production samples are introduced in this refactor.
