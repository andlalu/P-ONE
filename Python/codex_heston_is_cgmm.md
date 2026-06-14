# Codex Task: Minimal Heston Implied-State C-GMM Estimator

## Objective

Implement a clean, minimal, fast baseline Heston implied-state C-GMM estimator in the current Python codebase.

The target is **not** a wholesale port of the old MATLAB code. The target is a modern Python implementation that:

1. reuses the current simulation, option-pricing, characteristic-function / Riccati, and implied-volatility infrastructure wherever possible;
2. ports only the minimum methodological components needed from the old Boswijk--Laeven--Lalu MATLAB code;
3. implements a first working Heston implied-state C-GMM estimator for generated option-panel samples;
4. provides diagnostics and tests before extending to large Monte Carlo runs or second-step estimation.

The old MATLAB folder should be treated as a methodological reference, not as the architecture to reproduce.

---

## Source material to inspect

Inspect the following files first.

### Current report and model definition

- `Reports/18062026/draft_proposalJun26.tex`

Use this file as the authoritative source for the current Heston model, notation and parameter partition.

The current Heston parameter vector is:

\[
\theta = [\eta,\kappa,\bar V,\sigma_v,\rho,\eta_V].
\]

Under \(\mathbb P\), \(\eta\) enters the log-price drift. Under \(\mathbb Q\), \(\eta_V\) changes the variance drift, so that

\[
\kappa_{\mathbb Q}=\kappa-\eta_V,
\qquad
\bar V_{\mathbb Q}=\frac{\kappa\bar V}{\kappa-\eta_V}.
\]

The implementation must use this current model definition rather than the older MATLAB conventions.

### Current Python infrastructure

Inspect before adding new code:

- `Python/`
- `Python/DGPSimulation/`
- `Python/OptionPricing/`
- `Python/Scripts/configs/clean_generation_run_001.json`
- existing tests under `Python/**/tests/`

Look specifically for:

- Heston characteristic-function code;
- Riccati / affine transform code;
- COS option pricer code;
- implied-volatility solver interface;
- generated panel schema;
- noisy panel schema;
- conventions for config parsing, dataclasses, tests and package layout.

If characteristic-function or Riccati logic already exists for option pricing, reuse or refactor it. Do **not** create a parallel C-GMM-only Riccati solver unless reuse is impossible.

### Old MATLAB code: methodological reference

Inspect these files:

- `BLL_CODE/Code ReVIEW/Implying States/impl_state_func_fut.m`
- `BLL_CODE/Code ReVIEW/Derivative Pricing/ccf.m`
- `BLL_CODE/Code ReVIEW/Derivative Pricing/price_fft_fut.m`
- `BLL_CODE/Code ReVIEW/criterion.m`
- `BLL_CODE/Code ReVIEW/Cont CCF criterion/funcH0.m`
- `BLL_CODE/Code ReVIEW/Cont CCF criterion/critstep2.m`
- `BLL_CODE/Code ReVIEW/Cont CCF criterion/estimate.m`
- `BLL_CODE/Code ReVIEW/Cont CCF criterion/estimate2.m`
- `BLL_CODE/Code ReVIEW/Cont CCF criterion/boundtofree.m`
- `BLL_CODE/Code ReVIEW/Cont CCF criterion/freetobound.m`

Use these files to recover the estimator design:

- implied-state inversion from option-implied volatilities;
- date-by-date warm-starting of the latent state;
- conditional characteristic-function residuals;
- continuum GMM criterion;
- Gaussian quadrature / integration over CCF arguments;
- analytic integration over exponential instruments;
- first-step and optional second-step weighting;
- parameter constraints.

Do not port:

- old simulation code;
- old plotting/export code;
- old `.mat` workflows;
- old FFT option pricer unless needed for benchmarking;
- old `fminsearch` wrappers;
- duplicate experimental routines;
- old sine-based parameter transformation unless used only as a reference.

---

## Relevant methodology recovered from the MATLAB code

### Implied-state inversion

The old routine `impl_state_func_fut.m` implies latent states by minimizing Black--Scholes implied-volatility errors. The relevant structure is:

- one state for the stochastic-volatility-only case;
- starting value equal to the long-run variance;
- lower and upper bounds for the latent variance;
- first date initialized from long-run variance;
- subsequent dates initialized from the previous implied state;
- objective uses model-implied Black--Scholes IV minus observed IV.

For the Heston implementation, implement the scalar implied variance path

\[
\widehat V_i^\theta
=
\arg\min_{v\in[v_{\min},v_{\max}]}
\sum_{j=1}^{q_i}
\left(
\sigma_{i,j}^{\mathrm{obs}}
-
\sigma_{i,j}^{\mathrm{model}}
(S_i,v,K_i^{(j)},\tau_i^{(j)};\theta)
\right)^2.
\]

Default objective: IV-space least squares.

Optional later extensions:

- price-space errors;
- vega-scaled price errors;
- clean-IV state recovery for diagnostics.

### CCF residual

The old criterion uses a conditional characteristic-function residual of the form

\[
e_i(s;\theta)
=
\exp\left(i s^\top X_{i+1}^{\theta}\right)
-
\phi(s,X_i^\theta,\Delta;\theta),
\]

where

\[
\phi(s,X_i^\theta,\Delta;\theta)
=
\mathbb E_\theta
\left[
\exp\left(i s^\top X_{i+1}^{\theta}\right)
\mid X_i^\theta
\right].
\]

For the Heston case, use a two-dimensional estimation state, for example

\[
X_i^\theta = (r_i,\widehat V_i^\theta),
\qquad
r_i = \log S_i-\log S_{i-1}.
\]

Keep the state convention explicit and documented. The state convention must be consistent across the implied-state module, CCF module and criterion module.

### Continuum moment function

The continuum moment condition is indexed by \(\tau=(r,s)\):

\[
\widehat h_n(\tau;\theta)
=
\frac{1}{n-1}
\sum_{i=1}^{n-1}
\exp\left(i r^\top X_i^\theta\right)
e_i(s;\theta).
\]

At the true parameter value,

\[
\mathbb E_{\theta_0}
\left[
e_i(s;\theta_0)\mid X_i^{\theta_0}
\right]=0,
\]

so multiplying by instruments measurable at \(t_i\) preserves the moment restriction.

### First-step C-GMM criterion

The first-step criterion is the integrated squared norm

\[
Q_{1,n}(\theta)
=
\int
\left|
\widehat h_n(\tau;\theta)
\right|^2
\pi(\tau)\,d\tau.
\]

Following the BLL computational design, integrate out the instrument coordinate analytically under Gaussian instrument weights. This gives a Gaussian kernel

\[
K_r(X_i^\theta-X_j^\theta)
=
\exp\left(
-\frac12
(X_i^\theta-X_j^\theta)^\top
\Sigma_r
(X_i^\theta-X_j^\theta)
\right).
\]

The criterion becomes

\[
Q_{1,n}(\theta)
\approx
\frac{1}{(n-1)^2}
\sum_{k=1}^{K_s}
\omega_k
\sum_{i=1}^{n-1}
\sum_{j=1}^{n-1}
K_r(X_i^\theta-X_j^\theta)
e_i(s_k;\theta)\overline{e_j(s_k;\theta)}.
\]

Implement this with vectorized NumPy operations. For the current sample size, an \((n-1)\times(n-1)\) kernel matrix is acceptable.

### Optional second-step C-GMM

Do not implement this first unless the first-step estimator is stable.

The second step should:

1. compute a first-step estimate \(\widehat\theta_1\);
2. compute moment residuals at \(\widehat\theta_1\);
3. estimate the covariance/operator approximation;
4. apply a regularized inverse weighting object.

The old MATLAB code uses a regularized object of the form

\[
(\alpha I+C_m^2)^{-1}.
\]

In Python:

- do not form explicit matrix inverses;
- use linear solves, Cholesky, eigendecomposition or SVD;
- expose the regularization parameter \(\alpha\);
- return condition diagnostics.

---

## Implementation requirements

### 1. Add a minimal estimator package

Use the repo’s existing package conventions. A possible layout is:

```text
Python/Estimation/ISCGMM/
    __init__.py
    types.py
    heston_params.py
    implied_state.py
    heston_transition_cf.py
    quadrature.py
    cgmm_moments.py
    cgmm_criterion.py
    estimate.py
    diagnostics.py
    tests/
```

Adapt names to the existing code style.

### 2. Data input

The estimator should operate on existing generated panel outputs.

Required fields, or their equivalents:

- date/sample index;
- spot or log spot;
- strike;
- maturity;
- option type;
- clean price;
- clean IV;
- observed price;
- observed IV;
- vega if available;
- true variance if available for diagnostics.

Do not introduce a new incompatible data format unless necessary. If the existing generated panel schema is awkward, write a thin adapter.

### 3. Implied-state module

Implement:

```python
imply_heston_variance_path(theta, panel, config) -> ImpliedStateResult
```

Requirements:

- scalar bounded optimization per date;
- warm start from previous date;
- default \(V_{\min}=1e-8\), \(V_{\max}=1.0\), configurable;
- default objective: IV residuals;
- use current Heston pricing and IV infrastructure;
- avoid re-building static contract metadata on every objective call;
- return diagnostics:
  - objective values by date;
  - boundary hit flags;
  - failed inversion flags;
  - number of function evaluations if available;
  - implied variance path.

Performance guidance:

- for scalar \(V\), prefer `scipy.optimize.minimize_scalar(..., bounds=..., method="bounded")` or a similarly efficient scalar bounded method;
- do not use a general multidimensional least-squares solver for the one-state Heston case;
- cache any maturity/strike/type arrays that do not change within a date;
- use previous \(V\) as warm start only if the chosen optimizer supports it; otherwise use a narrow bracket around previous \(V\) as an optional refinement.

### 4. Heston transition CCF

Implement or expose:

```python
heston_p_transition_cf(s_nodes, x_prev, dt, theta) -> complex ndarray
```

where `s_nodes` are CCF evaluation points and `x_prev` contains the previous estimation states.

Before implementing, search for the Heston characteristic-function / Riccati code used by COS pricing. Reuse/refactor it so pricing and C-GMM share the same affine transform machinery.

For Heston-only C-GMM, do not implement a full generic AJD framework unless it already exists.

Sanity checks:

\[
\phi(0,X_i,\Delta;\theta)=1.
\]

Where applicable:

\[
\phi(-s,X_i,\Delta;\theta)=\overline{\phi(s,X_i,\Delta;\theta)}.
\]

### 5. Quadrature

Implement a small quadrature module for \(s\)-integration.

Initial acceptable default:

- tensor Gauss-Hermite or Gauss-Legendre nodes in two dimensions;
- configurable number of nodes per dimension;
- default low order suitable for tests;
- no dependency on the old MATLAB `nwspgr` unless a sparse-grid implementation already exists.

Expose:

```python
make_cgmm_quadrature(dim=2, order=..., scheme=...) -> nodes, weights
```

The first version should prioritize clarity and determinism.

### 6. First-step criterion

Implement:

```python
cgmm_first_step_criterion(theta_free_or_struct, panel, config) -> float
```

or a class-based equivalent.

The criterion should:

1. transform free parameters to constrained Heston parameters;
2. imply \(V_i^\theta\);
3. construct \(X_i^\theta=(r_i,\widehat V_i^\theta)\);
4. compute CCF residuals \(e_i(s_k;\theta)\);
5. construct Gaussian instrument kernel \(K_r\);
6. compute the integrated quadratic criterion;
7. return a real scalar;
8. return optional diagnostics in a separate method, not in the scalar optimizer callback.

Be careful with complex arithmetic:

- use conjugates explicitly;
- take the real part only at the end;
- guard against small negative numerical real parts if mathematically non-negative.

### 7. Parameter transformation and constraints

Do not copy the old sine transform.

Use a computationally simple, non-periodic parameterization:

\[
\kappa=\exp(a_\kappa),\qquad
\bar V=\exp(a_{\bar V}),\qquad
\sigma_v=\exp(a_\sigma),\qquad
\rho=\tanh(a_\rho).
\]

For the risk-neutral variance drift, enforce

\[
\kappa_{\mathbb Q}=\kappa-\eta_V>0
\]

by parameterizing

\[
\kappa_{\mathbb Q}=\exp(a_{\kappa Q}),
\qquad
\eta_V=\kappa-\kappa_{\mathbb Q}.
\]

For \(\eta\), use either an unconstrained value or a broad configurable bound.

Requirements:

- provide `to_free(theta)` and `from_free(free_theta)`;
- avoid periodic transforms;
- make transformations vectorized and cheap;
- provide tests for round-trip accuracy.

### 8. Optimizer

Initial optimizer:

- `scipy.optimize.minimize`;
- `method="L-BFGS-B"` or `method="Powell"`;
- deterministic objective;
- single-start first;
- configurable multi-start later.

Do not attempt gradients until the objective is stable.

### 9. Diagnostics

Implement diagnostic outputs:

- parameter vector;
- criterion value;
- state-inversion success rate;
- boundary-hit rate;
- implied \(V\) summary;
- true-vs-implied \(V\) comparison when true \(V\) is available;
- number of panel dates used;
- number of option contracts used;
- timing for:
  - state inversion;
  - CCF residuals;
  - kernel construction;
  - criterion evaluation.

### 10. Tests

Add focused tests.

#### Implied-state tests

- On clean generated panels and true parameters, implied variance should be close to true variance.
- Warm-start logic should be used or at least not regress date ordering.
- Boundary hits should be reported.

#### CCF tests

- \(\phi(0)=1\).
- finite complex output for representative \(s\)-nodes.
- conjugacy check where applicable.

#### Criterion tests

- deterministic criterion value for same input;
- finite criterion for true parameters;
- true parameters should usually have lower criterion than a strongly perturbed parameter on clean data;
- criterion should not mutate input panel data.

#### Parameter transform tests

- round-trip `theta -> free -> theta`;
- constraints enforced:
  - \(\kappa>0\);
  - \(\bar V>0\);
  - \(\sigma_v>0\);
  - \(\rho\in(-1,1)\);
  - \(\kappa-\eta_V>0\).

#### Second-step tests, if implemented

- no explicit inverse;
- regularized solve succeeds;
- diagnostics include eigenvalue or condition-number information.

---

## Performance constraints

The implementation must be fast enough to be usable inside Monte Carlo estimation loops.

Important constraints:

1. Do not re-price static metadata repeatedly.
2. Do not re-create quadrature nodes on every objective call.
3. Do not re-parse panel files inside the objective.
4. Do not use multidimensional least squares for scalar \(V\)-inversion.
5. Reuse/refactor Heston characteristic-function code from option pricing.
6. Avoid explicit matrix inverses.
7. Keep second-step optional.
8. Add profiling diagnostics before optimizing prematurely.

The main expected bottleneck is:

\[
\theta
\mapsto
\widehat V^\theta
\mapsto
X^\theta
\mapsto
Q_n(\theta),
\]

especially repeated option pricing inside implied-state inversion.

---

## Deliverables

### Required deliverable 1: first-step estimator

A working first-step Heston implied-state C-GMM estimator for one generated sample and one panel scenario.

Must include:

- Python modules;
- config object or JSON config;
- command-line or script entry point;
- tests;
- diagnostics;
- README or docstring explaining how to run on one sample.

### Required deliverable 2: validation run

Run a small validation:

1. load one clean or low-noise generated sample;
2. evaluate the criterion at the true parameter vector;
3. evaluate the criterion at at least two perturbed parameter vectors;
4. report implied-state diagnostics;
5. store a small JSON or text summary.

### Optional deliverable 3: second-step estimator

Only after first-step works:

- implement regularized two-step weighting;
- expose as `method="two_step"`;
- include numerical stability diagnostics.

---

## Non-goals for this task

Do not implement:

- particle filtering;
- full Monte Carlo batch estimation;
- asymptotic standard errors;
- bootstrap confidence intervals;
- old SVHJ model;
- Hawkes jumps;
- old FFT pricer as production code;
- plotting/reporting beyond minimal diagnostics;
- new generated data schema unless necessary.

---

## Suggested Codex work sequence

1. Inspect current Python infrastructure.
2. Identify or expose reusable Heston characteristic-function / Riccati functionality.
3. Implement Heston parameter transform.
4. Implement implied variance path inversion using existing pricing/IV code.
5. Implement quadrature nodes and weights.
6. Implement CCF residual computation.
7. Implement first-step criterion.
8. Add tests.
9. Add a one-sample runner and diagnostics.
10. Only then consider second-step weighting.

---

## Acceptance criteria

The task is acceptable when:

- tests pass;
- first-step criterion can be evaluated for one generated sample;
- implied variance recovery works on clean panels at true parameters;
- criterion is finite and deterministic;
- implementation reuses existing pricing/CF code or clearly documents why a small new helper was needed;
- no large legacy MATLAB structure is copied into Python;
- no explicit matrix inverse is used in production estimator code;
- code is modular enough to be called later in Monte Carlo loops.
