# Codex implementation specification: structural clean-up and naming consistency

## Repository and baseline

Repository:

```text
andlalu/P-ONE
```

Baseline reviewed commit:

```text
8d35cacfaa09f2b665ce2b6a7f03c43feca0188e
```

Work from the current branch HEAD. Before editing:

1. record the current Git SHA;
2. compare the affected files with the baseline above;
3. identify material changes since the baseline;
4. run a repository-wide call-site search before deleting or moving any function.

This task is a **structural clean-up**. Preserve the current numerical methods and research design.

---

# 1. Objective

Make the Python codebase internally consistent, easier to read and proportionate to a one-off doctoral research project.

The intended design is:

- concrete functions and classes rather than unused abstract interfaces;
- one clear production path for simulation, pricing, clean-panel generation and noise addition;
- files divided by logical concern;
- I/O placed consistently in I/O modules;
- configuration classes named and documented consistently;
- test-only helpers kept in tests;
- legacy production methods removed;
- concise mathematical comments where the implementation is not self-explanatory;
- no terminology that confuses two Powell passes with two-stage C-GMM estimation.

Do not turn the repository into a general-purpose framework.

---

# 2. Preserve numerical behaviour

Do not change:

- the Heston physical-measure simulation;
- Andersen QE formulas;
- the physical-to-risk-neutral parameter mapping;
- fixed-width COS pricing formulas;
- calibrated COS maturities, widths or term counts;
- Heston affine/Riccati formulas or complex-log branch convention;
- Black pricing, vega or implied-volatility inversion;
- clean-panel contract grid or option-type rules;
- any of the three noise formulas, parameters or seeds;
- price rounding and arbitrage-bound mechanics;
- analytical or finite-difference Jacobian formulas;
- implied-state objective or solver algorithms;
- the current production solver selection and fallback;
- the first-step C-GMM criterion;
- quadrature, instruments or transition characteristic function;
- parameter transforms, candidate starts or parameter bounds;
- the two consecutive Powell optimisation passes;
- output schemas unless a terminology-only column rename is explicitly required below.

The current production state solver remains:

```text
bounded_brent
```

with:

```text
golden_section
```

as fallback.

Do not investigate the analytical-Jacobian fallback cases in this task.

---

# 3. User directives to implement

The following directives should be followed closely.

## 3.1 Remove lightweight interface layers

> Clean up the interface functions for now and keep a simpler structure per project.

> Ditch `base.py` in `OptionPricing` — the interfaces in general are not used robustly and for this small-scale one-off research project the interfaces are unnecessary.

> Also in `DGPSimulation` ditch the base interface.

Delete:

```text
Python/OptionPricing/base.py
Python/DGPSimulation/base.py
```

Remove inheritance from:

- `CcfSolver`;
- `OptionPricer`;
- `OptionPriceCubeGenerator`;
- `OptionPriceCubeStore`;
- `ModelParams`;
- `SimulationConfig`;
- `SimulationPath`;
- `VarianceDrawer`;
- `PathSimulator`.

Do not replace them with:

- protocols;
- abstract base classes;
- registries;
- factories;
- dependency-injection containers.

Use concrete classes, explicit constructor types and ordinary functions.

Preserve the ability to inject a concrete variance drawer into `HestonPathSimulator` if currently used by tests or research code. Type it directly or with a small callable type alias; do not recreate an interface hierarchy.

---

## 3.2 Move physical-to-risk-neutral conversion to `HestonParameters`

Current function to remove from clean-panel code:

```python
def heston_q_from_p(
    params_p: HestonPhysicalParameters,
    eta_v: float,
) -> HestonRiskNeutralParameters:
    return HestonParameters(
        eta=params_p.eta,
        kappa=params_p.kappa,
        vbar=params_p.vbar,
        sigma_v=params_p.sigma_v,
        rho=params_p.rho,
        eta_v=eta_v,
        r=params_p.r,
        q=params_p.q,
    ).to_risk_neutral()
```

Add this class method to `HestonParameters`:

```python
@classmethod
def from_physical(
    cls,
    params: HestonPhysicalParameters,
    *,
    eta_v: float,
) -> HestonParameters:
    ...
```

Production use should then be explicit:

```python
params_q = HestonParameters.from_physical(
    params_p,
    eta_v=eta_v,
).to_risk_neutral()
```

Keep `to_physical`, `to_risk_neutral` and `with_rates`.

Add a short class-level explanation of:

```text
kappa_q = kappa - eta_v
vbar_q  = kappa * vbar / kappa_q
```

Do not duplicate this conversion logic elsewhere.

---

## 3.3 Add concise mathematical orientation to the Heston CCF

The CCF implementation should include a brief preamble explaining that Heston is exponential-affine in the current variance:

```math
\phi(u,\tau\mid V_t)
=
\exp\left(A(u,\tau)+B(u,\tau)V_t\right).
```

Explain in one or two sentences:

- `A` is the variance-independent affine coefficient;
- `B` loads the current variance;
- the branch-stable representation is retained to avoid complex-log branch jumps on wide transform grids.

Rename local variables:

```text
c       -> affine_a
d_term  -> affine_b
```

Return:

```python
CoefficientTensor(
    u_grid=u,
    maturities=t,
    cf_a=affine_a,
    cf_b=affine_b,
)
```

Keep the deterministic-variance `NotImplementedError`, but add a short comment explaining that the current closed form divides by `sigma_v**2` and that the deterministic limit requires a separate formula.

Do not change the Riccati mathematics.

---

## 3.4 Explain fixed COS basis reuse

Expand the documentation for `PreparedFixedCosBasis`.

It should briefly state that for one maturity the object stores:

- the COS frequency grid;
- payoff expansion terms;
- affine Heston `A` and `B` coefficients;
- the fixed effective width;
- the number of COS terms.

Explain that these objects are prepared once and reused across:

- observation dates;
- candidate latent variance values;
- state-solver evaluations.

Explain that fixed-basis reuse makes the variance derivative inexpensive:

```math
\frac{\partial}{\partial v}\exp(A+Bv)
=
B\exp(A+Bv).
```

Keep the English brief and simple.

---

# 4. Target package structure

Use the following structure unless a call-site audit identifies a strong reason for a small variation.

```text
Python/
├── DGPSimulation/
│   ├── config.py
│   ├── path.py
│   ├── heston_simulator.py
│   ├── variance_drawers.py
│   └── io.py
│
├── Models/
│   └── Heston/
│       ├── parameters.py
│       └── affine_transform.py
│
├── OptionPricing/
│   ├── config.py
│   ├── cos_basis.py
│   ├── cos_pricer.py
│   ├── heston_ccf.py
│   └── fft_pricer.py
│
├── OptionData/
│   ├── panel.py
│   ├── clean_panel.py
│   ├── noise_common.py
│   ├── noise_low_iid.py
│   ├── noise_spatial.py
│   ├── noise_persistent_factor.py
│   ├── add_noise.py
│   └── io.py
│
├── Estimation/
│   └── ISCGMM/
│       ├── config.py
│       ├── results.py
│       ├── implied_state.py
│       ├── implied_state_jacobian.py
│       ├── cgmm_criterion.py
│       ├── estimate.py
│       └── parameter_transform.py
│
└── Scripts/
    ├── experiment_config.py
    ├── generation.py
    └── ...
```

Do not create empty wrapper files merely to match this diagram.

---

# 5. DGP simulation clean-up

## 5.1 Delete historical benchmark

Delete:

```text
Python/DGPSimulation/benchmark.py
```

It is a historical microbenchmark and has no role in the current experiment.

## 5.2 Replace vague `types.py`

Move:

```text
HestonSimConfig
```

to:

```text
Python/DGPSimulation/config.py
```

Move:

```text
HestonPath
```

to:

```text
Python/DGPSimulation/path.py
```

Delete:

```text
Python/DGPSimulation/types.py
```

Update all imports.

This establishes a repository-wide convention:

```text
config.py   numerical/execution settings
path.py     simulated path containers
results.py  estimator outputs
io.py       persistence
```

## 5.3 Remove base-class inheritance

`HestonSimConfig` and `HestonPath` remain plain frozen dataclasses.

`HestonPathSimulator` remains a concrete dataclass or ordinary class.

`AndersenQeVarianceDrawer` and `EulerVarianceDrawer` remain concrete implementations with:

```python
draw_next_variance(...)
```

No abstract parent class.

---

# 6. Option-pricing clean-up

## 6.1 Delete `OptionPricing/base.py`

Remove all references and inheritance.

Concrete classes should expose their actual methods directly:

- `CosOptionPricer`;
- Heston CCF coefficient solver.

## 6.2 Rename the Heston CCF file and class parsimoniously

Preferred:

```text
OptionPricing/heston_ccf.py
```

with a concrete function or small concrete class.

Acceptable options:

```python
def heston_pricing_coefficients(...)
```

or:

```python
class HestonCcf:
    def coefficients(...)
```

Do not retain a generic `CcfSolver` abstraction when only Heston is used.

Avoid `Any` in the public numerical signature. Use NumPy arrays and `HestonRiskNeutralParameters`.

## 6.3 Dissolve `OptionPricing/types.py`

Remove generic and legacy containers that no longer have production call sites.

Move genuine value objects next to the code that owns them:

- `CoefficientTensor` near the Heston CCF implementation;
- `PreparedFixedCosBasis` near COS basis preparation;
- `FixedBasisPriceJacobian` near the COS pricer/Jacobian implementation.

Candidates to delete after call-site audit:

```text
VarianceScaledCosConfig
OptionPriceCubeConfig
OptionPriceCube
PricingStack
```

Delete `OptionPricing/types.py` when empty.

## 6.4 Delete legacy dense-cube panel path

Audit and, when confirmed unused by production code, delete:

```text
Python/OptionPricing/panel_generator.py
```

and:

```text
HestonOptionPriceCubeGenerator
OptionPriceCube
OptionPriceCubeConfig
PricingStack
```

The current production path is the generated market-data-style panel. Do not maintain a second dense price-cube architecture without a production consumer.

---

# 7. Clean-panel generation

Move the production clean-panel generation to:

```text
Python/OptionData/clean_panel.py
```

This module should contain only:

- panel column definitions if they are specific to clean-panel generation;
- option-type selection from log moneyness;
- generation of clean panel records from one Heston path;
- concise transformation logic from path states to strikes, prices, IVs and vegas.

Do not include:

- CSV writing;
- Parquet writing;
- metadata-sidecar writing;
- generic table reading;
- historical dense cube generation;
- physical-to-risk-neutral conversion logic.

Delete the old duplicate/legacy panel-generation path rather than merging two implementations into one larger file.

Update `Scripts/generation.py` to call the single production clean-panel function.

---

# 8. Split noise code by scenario

Replace the current large `OptionPricing/noisy_panel.py`.

Move noise-panel functionality under `OptionData`, because it transforms generated option data rather than implementing an option-pricing model.

Use:

```text
OptionData/noise_common.py
OptionData/noise_low_iid.py
OptionData/noise_spatial.py
OptionData/noise_persistent_factor.py
OptionData/add_noise.py
```

## 8.1 Responsibilities

### `noise_common.py`

Contain only shared concepts:

- common marginal scale;
- deterministic scenario seed;
- common price bounds;
- tick rounding and price-bound mechanics;
- conversion from noisy IV to observed price and back to observed IV;
- shared result dataclass if genuinely needed;
- shared concise terminology.

### `noise_low_iid.py`

Contain only the low-IID noise formula.

### `noise_spatial.py`

Contain only the spatial-correlation noise formula and its Cholesky/eigenvalue fallback.

### `noise_persistent_factor.py`

Contain only:

- persistent-factor parameter checks that are scientifically necessary;
- stationary variance/innovation variance conversion;
- factor recursion;
- residual-scale formula;
- factor-record construction.

### `add_noise.py`

Contain the small orchestration layer:

- select scenario;
- create RNG;
- call one scenario implementation;
- apply shared price mechanics;
- return noisy rows and optional factor records.

Do not recreate a plugin framework or registry. A clear `if/elif` dispatch over three fixed scenarios is acceptable.

## 8.2 Terminology

Avoid:

```text
contaminate
contaminated
contamination
```

Use:

```text
add_noise
noisy
noise
raw_noisy_iv
noise_scenario
```

Rename, where practical without breaking the scientific output contract:

```text
_contaminate_rows          -> apply_noise_to_rows
raw_contaminated_iv        -> raw_noisy_iv
```

Update internal variables, function names, comments, logs and tests.

For persisted panel columns, prefer renaming `raw_contaminated_iv` to `raw_noisy_iv` now, before production samples are frozen. Update readers and tests accordingly. Do not retain duplicate alias columns.

---

# 9. Apply a consistent I/O discipline

Retain I/O modules, but use them consistently.

## 9.1 `DGPSimulation/io.py`

Own only:

- save Heston path NPZ;
- load Heston path NPZ;
- path-format metadata.

## 9.2 `OptionData/io.py`

Own all panel and panel-adjacent persistence:

- read CSV/Parquet records;
- write CSV/Parquet records atomically;
- load an `OptionPanel`;
- read/write metadata sidecars;
- locate clean/noisy panel files;
- write persistent-factor CSV files.

Remove from computational modules:

- `parquet_available`;
- `write_panel`;
- `read_panel_csv`;
- `read_table`;
- `write_table`;
- `write_factor_file`.

A computational module may return rows or factor records; it should not decide how they are stored.

## 9.3 Run-level orchestration output

`Scripts/generation.py` may own:

- completion markers;
- run-level generation manifest;
- run metadata.

Alternatively, a small `Scripts/generation_io.py` may own them if this materially improves readability. Do not add another layer merely for symmetry.

---

# 10. Remove test-only and legacy production methods

Perform a repository-wide call-site audit.

Delete from production code when no non-test caller exists:

```text
price_one_variance_scaled_reference
price_matrix_variance_scaled_reference
variance_scaled_effective_width
VarianceScaledCosConfig
```

These are explicitly legacy/reference paths and should not remain in `CosOptionPricer`.

Where a test still requires a small reference calculation:

- implement the helper locally in the test file; or
- use the independent Fourier reference in `fft_pricer.py`.

Do not retain production methods solely to support one test.

Repeat the audit for other public or module-level methods:

1. identify functions called only from tests;
2. decide whether the test is still methodologically valuable;
3. move a required helper into the test;
4. otherwise delete the helper and refactor/remove the test.

Do not delete genuine numerical production functions simply because their call is indirect.

Include the audit findings in the final Codex report.

---

# 11. Make COS configuration validation parsimonious

Move `FixedCosBasisConfig` to:

```text
Python/OptionPricing/config.py
```

Retain only checks that prevent realistic scientific/configuration errors:

- non-empty aligned maturities and effective widths;
- positive finite maturities and widths;
- positive generation and estimation term counts;
- positive maturity tolerance;
- unique configured maturities;
- requested maturity is present.

Keep one straightforward panel/estimator compatibility check covering:

- maturity grid;
- effective widths;
- generation term count.

Do not compare estimation term count because it may legitimately differ.

Remove duplicated validation passes and defensive checks for internal states that are already constructed by the same validated code path.

Simplify tests accordingly:

- one valid lookup test;
- one unknown-maturity test;
- one mismatched configuration test;
- one panel compatibility test.

Do not retain an extensive boundary-condition matrix for a fixed production experiment.

---

# 12. Configuration naming and documentation

## 12.1 Naming convention

Use `config.py` consistently for numerical setting objects:

```text
DGPSimulation/config.py
OptionPricing/config.py
Estimation/ISCGMM/config.py
```

Use specific container names rather than generic `types.py`.

## 12.2 Add brief field documentation

Add concise comments or docstrings explaining the role of each setting.

### `ImpliedStateConfig`

Explain:

- `cos_basis`: fixed maturity-specific COS widths and term counts;
- `v_min`, `v_max`: admissible latent-variance interval;
- `tol`: inner scalar solver convergence tolerance;
- `max_iter`: maximum inner iterations/evaluations as used by the chosen method;
- `boundary_tol`: distance used to flag a solution near a bound;
- `warm_start_window`: local interval around the preceding date's estimate;
- `state_solver`: primary latent-state inversion algorithm;
- `fallback_solver`: scalar fallback after a derivative-based method fails;
- finite-difference steps: perturbation sizes for numerical Jacobian validation/solver;
- `minimum_black_vega`: contracts below this vega are unsafe for price-to-IV differentiation.

### `CcfQuadratureConfig`

Explain:

- `dimension`: dimension of the integration-node vector;
- `order`: Gauss-Hermite order per dimension;
- `scale`: scale of the quadrature nodes.

### `CgmmConfig`

Explain:

- `instrument_precision`: frequency scales used in the C-GMM instruments;
- `transition_rk_steps`: RK4 resolution when the numerical transition CCF is selected;
- `transition_cf_method`: analytic or RK4 transition CCF;
- `dt`: optional fixed observation spacing; otherwise inferred;
- `spacing_tolerance`: tolerance for checking an equally spaced panel.

### Powell configuration

Explain:

- `base_start`: central starting parameter vector;
- `natural_bounds`: bounds in interpretable Heston coordinates;
- `candidate_relative_perturbations`: alternative starts screened before Powell;
- `coarse_pass`: broad initial Powell search;
- `refinement_pass`: tighter Powell search starting from the coarse result;
- `penalty_value`: finite objective returned after an expected numerical failure;
- `progress_every`: objective evaluations between progress messages.

Keep comments short and written in simple English.

## 12.3 Configuration README

The JSON production configuration cannot contain comments.

Add:

```text
Python/Scripts/configs/README.md
```

Document:

- each top-level section;
- each estimation setting;
- the distinction between first-step C-GMM and the two Powell passes;
- the fixed COS generation/estimation term counts;
- noise scenario parameters.

Keep it concise.

---

# 13. Rename the two Powell stages

The two current optimisation runs are **not** first- and second-stage C-GMM estimators.

They are two Powell passes over the same first-step criterion:

```text
coarse Powell pass
        ↓
refinement Powell pass
```

Keep both numerical passes and their current settings.

Remove ambiguous names.

## 13.1 Python names

Rename:

```text
PowellStageConfig      -> PowellPassConfig
OptimizerConfig        -> PowellConfig
stage1                 -> coarse_pass
stage2                 -> refinement_pass
OptimizerStageResult   -> PowellPassResult
stage_results          -> powell_passes
stage_number           -> pass_name or pass_index
```

Prefer log messages:

```text
Powell coarse pass started
Powell coarse pass completed
Powell refinement pass started
Powell refinement pass completed
```

Do not write:

```text
Stage 1
Stage 2
second stage
two-stage optimisation
```

unless explicitly referring to econometric two-stage C-GMM, which is not implemented.

## 13.2 JSON names

Replace:

```json
"stages": [
  {...},
  {...}
]
```

with:

```json
"coarse_pass": {
  "max_evaluations": 120,
  "xtol": 0.02,
  "ftol": 0.002
},
"refinement_pass": {
  "max_evaluations": 300,
  "xtol": 0.0002,
  "ftol": 0.00002
}
```

Preserve the exact values.

## 13.3 Results and metadata

Change serialised result names consistently:

```text
stage_results   -> powell_passes
stage           -> pass_name
```

Expected values:

```text
coarse
refinement
```

This output schema has not yet been used for the final Monte Carlo run, so prefer a clean rename rather than compatibility aliases.

Keep:

```text
FirstStepEstimate
estimate_first_step
CgmmFirstStepCriterion
```

because they correctly refer to first-step C-GMM estimation.

---

# 14. Comments and readability

Add mathematical or conceptual comments only where they materially help a reader:

- physical-to-risk-neutral mapping;
- affine Heston CCF;
- branch-stable complex logarithm;
- fixed COS basis reuse;
- analytical initial-variance derivative;
- Andersen QE regime split;
- persistent-factor recursion and loading;
- first-step C-GMM versus Powell-pass distinction.

Do not add comments that merely restate Python syntax.

Use concise British English.

Prefer descriptive numerical names:

```text
affine_a
affine_b
candidate_variance
coarse_pass
refinement_pass
raw_noisy_iv
```

Avoid single-letter names outside equations or short local linear-algebra expressions.

---

# 15. Tests

## 15.1 Preserve numerical coverage

Retain and update tests for:

- Heston parameter mapping;
- Heston path determinism and dimensions;
- Andersen QE non-negative variance;
- fixed COS pricing;
- analytical initial-variance Jacobian;
- clean-panel generation;
- all three noise scenarios;
- panel I/O;
- implied-state recovery;
- C-GMM criterion;
- Powell estimation.

## 15.2 Remove structural tests for deleted architecture

Delete or rewrite tests that assert:

- abstract base-class inheritance;
- generic pricing stacks;
- dense option-price cube classes;
- legacy variance-scaled production methods;
- excessive COS validation cases;
- obsolete import re-exports;
- old `stage1`/`stage2` result names;
- old contamination terminology.

## 15.3 Add focused tests

Add tests for:

1. `HestonParameters.from_physical(...).to_risk_neutral()`;
2. clean-panel generation uses this mapping;
3. clean-panel and noisy-panel I/O use `OptionData/io.py`;
4. each noise module reproduces the previous deterministic output for a fixed seed;
5. the unified noise orchestrator dispatches the three scenarios;
6. Powell coarse and refinement passes preserve the previous sequence and exact settings;
7. result serialisation uses `powell_passes` and `pass_name`;
8. no production import references either removed `base.py`;
9. no production call references legacy variance-scaled COS methods;
10. production configuration loads the renamed `coarse_pass` and `refinement_pass`.

Tests 8 and 9 may be implemented as simple source/import audits only if they remain short and robust. Do not build a custom architecture-testing framework.

---

# 16. Behavioural-equivalence checks

Before editing, capture deterministic reference outputs for:

1. one short Heston path;
2. one short clean panel;
3. each of the three noise scenarios;
4. one fixed-basis COS price grid;
5. one analytical Jacobian grid;
6. one short first-step Powell estimation.

After editing, compare:

- arrays with appropriate numerical tolerances;
- panel rows and ordering;
- noise draws and seeds exactly;
- output columns, except the deliberate `raw_contaminated_iv` to `raw_noisy_iv` rename;
- parameter estimates and criterion values within numerical tolerance;
- Powell coarse/refinement evaluation sequence and settings.

Do not claim numerical equivalence without running these checks.

---

# 17. Suggested implementation order

## Phase 1 — inventory and reference capture

- inspect HEAD versus baseline;
- map all production and test call sites;
- capture deterministic outputs;
- identify every interface, legacy path and test-only production helper.

## Phase 2 — DGP and parameter model

- add `HestonParameters.from_physical`;
- split DGP config/path files;
- remove DGP base classes;
- delete benchmark;
- update imports and tests.

## Phase 3 — pricing

- remove `OptionPricing/base.py`;
- simplify the Heston CCF;
- move pricing value objects to their owning modules;
- remove dense cube path;
- remove variance-scaled legacy methods;
- simplify COS configuration and validation;
- add mathematical comments.

## Phase 4 — option-panel organisation

- move clean-panel generation to `OptionData`;
- centralise panel I/O;
- split the three noise models;
- replace contamination terminology;
- update generation orchestration.

## Phase 5 — estimation naming

- rename Powell stages to coarse/refinement passes;
- update Python config, JSON config, results, logging and tests;
- add concise setting documentation and configuration README.

## Phase 6 — verification

- run focused tests after each phase;
- run the full Python test suite;
- run deterministic equivalence checks;
- run one temporary sample-000 generation smoke test;
- run a short state-only and first-step estimation smoke test.

---

# 18. Files likely to change

Likely additions:

```text
Python/DGPSimulation/config.py
Python/DGPSimulation/path.py
Python/OptionPricing/config.py
Python/OptionPricing/heston_ccf.py
Python/OptionData/clean_panel.py
Python/OptionData/noise_common.py
Python/OptionData/noise_low_iid.py
Python/OptionData/noise_spatial.py
Python/OptionData/noise_persistent_factor.py
Python/OptionData/add_noise.py
Python/Scripts/configs/README.md
```

Likely deletions:

```text
Python/DGPSimulation/base.py
Python/DGPSimulation/types.py
Python/DGPSimulation/benchmark.py
Python/OptionPricing/base.py
Python/OptionPricing/types.py
Python/OptionPricing/panel_generator.py
Python/OptionPricing/clean_panel.py
Python/OptionPricing/noisy_panel.py
```

Likely modifications:

```text
Python/Models/Heston/parameters.py
Python/Models/Heston/affine_transform.py
Python/OptionPricing/cos_basis.py
Python/OptionPricing/cos_pricer.py
Python/OptionData/io.py
Python/Estimation/ISCGMM/config.py
Python/Estimation/ISCGMM/results.py
Python/Estimation/ISCGMM/estimate.py
Python/Scripts/experiment_config.py
Python/Scripts/generation.py
Python/Scripts/configs/heston_experiment_run_001.json
affected scripts and tests
```

Adjust this list based on actual call sites. Avoid unrelated changes.

---

# 19. Acceptance criteria

The task is complete only when:

1. `OptionPricing/base.py` and `DGPSimulation/base.py` are deleted;
2. no replacement abstract-interface framework is introduced;
3. `HestonParameters.from_physical` owns the \(P\)-to-\(Q\) construction;
4. clean-panel code no longer contains parameter-conversion or file-I/O logic;
5. Heston CCF code contains a brief affine-equation preamble;
6. affine locals are named `affine_a` and `affine_b`;
7. `PreparedFixedCosBasis` clearly documents reuse and the variance derivative;
8. the historical DGP benchmark is deleted;
9. DGP `config.py` and `path.py` replace vague `types.py`;
10. the legacy dense option-price cube path is deleted when confirmed unused;
11. the three noise models live in separate files;
12. contamination terminology is removed from current code and persisted columns;
13. panel and factor persistence is centralised in `OptionData/io.py`;
14. no production method exists solely for one test;
15. variance-scaled legacy COS methods and config are removed;
16. COS validation is materially smaller but retains essential scientific checks;
17. configuration field roles have concise documentation;
18. a configuration README exists;
19. both Powell passes remain numerically unchanged;
20. Powell terminology uses `coarse_pass` and `refinement_pass`;
21. no Powell code or result uses ambiguous `stage1`, `stage2` or `second stage`;
22. first-step C-GMM names remain intact;
23. deterministic simulation, pricing, panel and noise outputs are preserved;
24. the analytical Jacobian and criterion are numerically unchanged;
25. the current production state solver and fallback are unchanged;
26. focused tests and the full feasible test suite pass;
27. sample-000 smoke generation and short first-step estimation pass;
28. no production samples `001`–`099` are generated;
29. report files are unchanged;
30. the final diff contains no unrelated formatting or redesign.

---

# 20. Final Codex report

Codex must report:

1. starting and ending Git SHA;
2. implementation plan actually followed;
3. files added, moved, changed and deleted;
4. interfaces removed;
5. legacy/test-only methods removed;
6. call-site audit findings;
7. before/after package structure;
8. I/O ownership after the refactor;
9. noise-module split;
10. exact terminology renames;
11. confirmation that both Powell passes remain and their settings are unchanged;
12. confirmation that this is still first-step C-GMM only;
13. confirmation that the production state solver/fallback are unchanged;
14. deterministic equivalence checks and tolerances;
15. exact test commands and results;
16. exact smoke-test commands and results;
17. any checks not run;
18. any remaining structural inconsistencies;
19. confirmation that reports and samples `001`–`099` were not changed.
