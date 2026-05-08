# Heston Option Pricing Stack — Technical Design

## Purpose
Define a production-ready, computationally efficient pricing stack that can:
1. Solve characteristic-function coefficients (CCF) via an analytic Heston implementation with stable complex-log branching.
2. Price strike/maturity matrices via COS-FFT style workflows with coefficient reuse.
3. Build weekly option price panels from simulated Heston state paths.
4. Persist/reload priced panels with reproducibility metadata.

This design intentionally follows the interface-first naming and boundary style used in `Python/sim` (`PathSimulator`, `VarianceDrawer`) so pricing components remain model-agnostic where possible.

---

## Current state (gap summary)
- `Python/OptionPricing/cos_transform_pricer.py` currently computes only partial COS payoff terms and does not expose a full model/pricing API.
- No dedicated interface layer exists for CCF solving vs pricing vs orchestration.
- No panel-oriented save/load API exists for option prices.
- No comprehensive validation suite yet for:
  - Monte Carlo (MC) cross-checks.
  - Black–Scholes (BS) limit/edge-case consistency.

---

## Proposed architecture

### 1) Base interfaces / contracts
Create a `Python/OptionPricing/base.py` module with abstract contracts:

- `CcfSolver` (interface/base class)
  - Responsibility: return CF/CCF coefficients on a requested transform grid.
  - Method shape (illustrative):
    - `solve_coefficients(u_grid, maturity_grid, model_params, measure_params) -> coefficients`

- `OptionPricer` (interface/base class)
  - Responsibility: consume coefficients and return prices.
  - Method shape:
    - `price_matrix(state_matrix, strike_grid, maturity_grid, rate_grid, coefficients, pricing_config) -> price_matrix`

- `OptionPanelGenerator` (interface/base class)
  - Responsibility: orchestrate state paths + solver + pricer into indexed panel output.

- `OptionPanelStore` (interface/base class)
  - Responsibility: persistence boundary for save/load.

Naming rationale: avoid `I*` prefixes and keep style compatible with existing `sim` package naming.

### 2) Data types
Create `Python/OptionPricing/types.py` with dataclasses:

- `HestonPricingParamsQ` (risk-neutral parameters, includes variance risk premium convention used by solver).
- `CosPricingConfig` (`n_cos`, truncation strategy, bounds/cumulant policy, branch policy).
- `OptionPanelConfig` (strike convention, maturities, rate handling, return type).
- `OptionPanel` (tensor + explicit axis metadata + provenance fields).
- `CoefficientTensor` (or structured container with named slices and shape validation).

### 3) Heston analytic CCF solver
Create `Python/OptionPricing/heston_ccf_solver.py` implementing `CcfSolver`:

- `HestonAnalyticCcfSolver`:
  - Analytic Heston coefficients under Q.
  - Numerically stable complex square root/log branch logic.
  - Supports precomputation over all maturities in one call.

Implementation notes:
- Use branch handling consistent with "Little Heston Trap" / Lord-Kahl guidance.
- Include invariant checks (no NaN/Inf, continuity diagnostics over `u_grid`).

### 4) COS matrix pricer
Create `Python/OptionPricing/cos_pricer.py` implementing `OptionPricer`:

- `CosOptionPricer`:
  - Precompute static COS terms (`k_grid`, payoff coefficients, weights).
  - Vectorized pricing over strike grid for each maturity/state batch.
  - Accept optional precomputed coefficients for repeated evaluations.

Optimization strategy:
- Reuse coefficient slices for all strikes at a maturity.
- Avoid Python loops over strikes.
- Keep contiguous array layouts and minimal temporary allocations.

### 5) Panel orchestration
Create `Python/OptionPricing/panel_generator.py` implementing `OptionPanelGenerator`:

- `HestonOptionPanelGenerator`:
  - Input: weekly state paths (`logS_week`, `V_week`) + strike/maturity specs.
  - Flow:
    1. Build transform grid.
    2. Solve/load coefficients (cached by maturity/config hash).
    3. Price matrix for each observation date (or in batches).
    4. Return `OptionPanel` with explicit indexes and metadata.

### 6) Panel persistence
Create `Python/OptionPricing/io.py` implementing `OptionPanelStore`:

- `NpzOptionPanelStore`:
  - Save/load `OptionPanel` using NPZ + JSON metadata.
  - Metadata includes: format version, pricing config, model params, source simulation seed, and axis semantics.

Optional future path:
- Add chunked backend (`zarr`) behind same interface when panel size grows.

---

## Testing strategy (required additions)

## A) Unit tests
- Shape/index contracts for all public interfaces.
- Determinism where seeded and deterministic code paths are used.
- Input validation and error-path tests (bad grid shapes, negative maturities, etc.).

## B) Numerical tests — analytic/COS layer
- Regression fixtures from MATLAB reference outputs for selected parameter sets.
- Convergence checks in `n_cos` (monotone/stabilizing price error behavior).

## C) Monte Carlo cross-check tests (explicitly requested)
Use existing Heston simulator in `Python/sim` as benchmark engine:
- For selected strikes/maturities:
  1. Simulate under risk-neutral dynamics consistent with pricing params.
  2. Compute discounted payoff estimates and confidence intervals.
  3. Verify COS prices are within MC CI or agreed tolerance bands.

Test design:
- Keep one "fast" CI-sized suite and one "slow" high-precision suite (marked/optional).

## D) Black–Scholes edge/limit tests (explicitly requested)
Add tests for limiting behavior where Heston should collapse to BS-like pricing:
- `sigma_v -> 0` and stable variance path near constant level.
- `rho = 0` sanity checks.
- short-maturity behavior vs BS reference where appropriate.

Acceptance criteria:
- COS-Heston price differences vs BS baseline remain within tolerance under defined limit settings.

## E) Performance tests / benchmarks
- Measure:
  - coefficient solve time per maturity grid,
  - matrix pricing time per weekly observation,
  - full panel generation time for ~500 observations.
- Report with and without coefficient reuse cache.

---

## Suggested module layout

- `Python/OptionPricing/base.py`
- `Python/OptionPricing/types.py`
- `Python/OptionPricing/heston_ccf_solver.py`
- `Python/OptionPricing/cos_pricer.py`
- `Python/OptionPricing/panel_generator.py`
- `Python/OptionPricing/io.py`
- `Python/OptionPricing/benchmark.py`
- `Python/OptionPricing/tests/`
  - `test_ccf_solver.py`
  - `test_cos_pricer_shapes.py`
  - `test_cos_pricer_mc_parity.py`
  - `test_cos_pricer_bs_limits.py`
  - `test_panel_generator.py`
  - `test_option_panel_io.py`

---

## External references (implementation guidance)

1. Fang, F. & Oosterlee, C.W. (2008), COS method foundational paper.
2. Albrecher, H., Mayer, P., Schoutens, W., Tistaert, J. (2007), "The Little Heston Trap".
3. Lord, R. & Kahl, C. (2010 / preprint 2008), complex logarithms in Heston-like models.
4. Le Floc’h, F. (2017), stable Fourier-COS call pricing refinements.

These references should drive branch handling, truncation selection, and stability checks in implementation.

---

## PR sequencing plan

### PR 1 — API skeleton + types + tests for contracts
No pricing math yet; establish stable public boundaries.

### PR 2 — Heston analytic CCF solver + branch-stability tests
Include regression fixtures and continuity diagnostics.

### PR 3 — Vectorized COS pricer + strike/maturity matrix support
Include coefficient reuse path and performance benchmark.

### PR 4 — Panel generation + IO + MC/BS validation suite
Add end-to-end pipeline and persistence.

This phased approach minimizes integration risk and makes numerical issues easier to isolate.
