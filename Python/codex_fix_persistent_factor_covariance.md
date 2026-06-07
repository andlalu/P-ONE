# Codex task: verify and fix persistent-factor noise covariance semantics

## Goal

Investigate and, if needed, fix the implementation of the `persistent_factor` contaminated option-panel noise design so that the configuration semantics match the report specification.

The key issue is:

> `q_diag` must mean the diagonal of the **innovation covariance matrix** \(Q\), not a vector of innovation standard deviations and not the stationary factor standard deviations.

The user-facing configuration may also provide `stationary_factor_std`. These are the unconditional standard deviations of the persistent common factors. If both `stationary_factor_std` and `q_diag` are present, the code must either validate that they are consistent or use one as the source of truth in a clearly documented way.

Do not implement or modify IS estimation, PF estimation, AWS Batch, Docker, Terraform, Student-t noise, mixture noise, or stale quote logic in this task.

---

## Report specification to implement

For Design C, the common factor process is

\[
f_i=A f_{i-1}+\xi_i,
\qquad
\xi_i\sim\mathcal{N}(0,Q),
\]

where

\[
A=\mathrm{diag}(a_0,a_m,a_\tau)
\]

and \(Q\) is the innovation covariance matrix.

The report calibration is

\[
A=\mathrm{diag}(0.85,0.65,0.65).
\]

The stationary standard deviations of the three factor components are

\[
(\bar s_{f,0},\bar s_{f,m},\bar s_{f,\tau})
=
(0.0007,\ 0.0025,\ 0.0008).
\]

Therefore, if \(Q\) is diagonal,

\[
Q
=
\mathrm{diag}
\left(
(1-a_0^2)\bar s_{f,0}^2,\,
(1-a_m^2)\bar s_{f,m}^2,\,
(1-a_\tau^2)\bar s_{f,\tau}^2
\right).
\]

Numerically this is

\[
Q=
\mathrm{diag}
\left(
1.35975\cdot10^{-7},\,
3.609375\cdot10^{-6},\,
3.696\cdot10^{-7}
\right).
\]

These values are **variances**, not standard deviations.

---

## Desired config semantics

The config should support the following persistent-factor block:

```json
"persistent_factor": {
  "enabled": true,
  "alpha_0": 0.0015,
  "alpha_m": 2.0,
  "alpha_tau": 0.05,
  "tau_min": 0.01984126984126984,
  "a_diag": [0.85, 0.65, 0.65],
  "stationary_factor_std": [0.0007, 0.0025, 0.0008],
  "q_diag": [1.35975e-7, 3.609375e-6, 3.696e-7],
  "residual_policy": "match_total_marginal_scale",
  "factor_initialization": "zero"
}
```

### Required interpretation

- `a_diag`: diagonal entries of \(A\).
- `stationary_factor_std`: stationary/unconditional standard deviations of \(f_i\).
- `q_diag`: diagonal entries of the innovation covariance matrix \(Q\).
- `q_diag` entries are **variances**.
- when drawing innovations,
  \[
  \xi_{k,i}\sim\mathcal{N}(0,q_k)
  \]
  so the implementation should use `sqrt(q_diag[k])` as the normal draw multiplier.

### Source-of-truth policy

Prefer the following behavior:

1. If `stationary_factor_std` is present and `q_diag` is absent:
   - compute `q_diag` from
     \[
     q_k=(1-a_k^2)\bar s_{f,k}^2.
     \]

2. If both `stationary_factor_std` and `q_diag` are present:
   - compute the implied `q_diag` from `a_diag` and `stationary_factor_std`;
   - validate that the supplied `q_diag` agrees within a small numerical tolerance;
   - raise a clear `ValueError` if inconsistent.

3. If `q_diag` is present and `stationary_factor_std` is absent:
   - interpret `q_diag` as innovation covariance diagonal;
   - compute the implied stationary factor standard deviations as diagnostics:
     \[
     \bar s_{f,k} = \sqrt{\frac{q_k}{1-a_k^2}}.
     \]

4. Reject invalid values:
   - `a_diag` must have length 3;
   - `stationary_factor_std`, if present, must have length 3 and non-negative entries;
   - `q_diag`, if present, must have length 3 and non-negative entries;
   - all absolute values of `a_diag` must be strictly below 1 for stationarity.

---

## Residual scale policy for Design C

The design goal is that Design C has the same unconditional pointwise implied-volatility observation error scale as Designs A and B.

Let

\[
\Omega_f
=
\mathrm{diag}
\left(
\bar s_{f,0}^2,\bar s_{f,m}^2,\bar s_{f,\tau}^2
\right).
\]

For each option contract \(j\),

\[
b_i^{(j)}
=
(1,m_i^{(j)},\tau_i^{(j)})^\top.
\]

The factor contribution has unconditional variance

\[
\left(b_i^{(j)}\right)^\top\Omega_f b_i^{(j)}.
\]

If the marginal Design A/B scale is \(s_{\sigma,i}^{(j)}\), the idiosyncratic residual scale should be

\[
s_{\sigma,u,i}^{(j)}
=
\sqrt{
\max\left[
\left(s_{\sigma,i}^{(j)}\right)^2
-
\left(b_i^{(j)}\right)^\top\Omega_f b_i^{(j)},
0
\right]
}.
\]

Implement `residual_policy = "match_total_marginal_scale"` using this formula.

This replaces any older interpretation such as:

```json
"residual_scale_multiplier": 0.50
```

Do not keep `residual_scale_multiplier` in the canonical config unless it is needed for backwards compatibility. If backwards compatibility is retained, document that it is deprecated and not used when `residual_policy = "match_total_marginal_scale"`.

---

## Expected implementation work

1. Locate the contaminated panel noise implementation.
   - Search for `persistent_factor`, `q_diag`, `a_diag`, `factor_initialization`, and `residual_scale_multiplier`.
   - Determine whether `q_diag` is currently treated as:
     - innovation covariance diagonal,
     - innovation standard deviation,
     - stationary factor standard deviation,
     - or something else.

2. Update the parser / config dataclass / validation logic.
   - Make the semantics above explicit.
   - Add helper functions if useful:
     - `compute_q_diag_from_stationary_std(a_diag, stationary_factor_std)`
     - `compute_stationary_std_from_q_diag(a_diag, q_diag)`
     - `validate_persistent_factor_config(...)`.

3. Update the factor innovation draw.
   - If \(q_k\) is the innovation variance, the innovation draw must be:
     ```python
     xi = rng.normal(loc=0.0, scale=np.sqrt(q_diag), size=3)
     ```
   - Do not use `q_diag` directly as the normal scale.

4. Update residual scale computation.
   - Implement the `match_total_marginal_scale` policy.
   - Ensure the factor variance is computed from `stationary_factor_std`, not from `q_diag` directly unless stationary standard deviations are inferred from it.

5. Update the canonical run config.
   - Use:
     ```json
     "a_diag": [0.85, 0.65, 0.65],
     "stationary_factor_std": [0.0007, 0.0025, 0.0008],
     "q_diag": [1.35975e-7, 3.609375e-6, 3.696e-7],
     "residual_policy": "match_total_marginal_scale",
     "factor_initialization": "zero"
     ```
   - Remove or deprecate:
     ```json
     "residual_scale_multiplier": 0.50
     ```

6. Update documentation / README comments if the repo has them.
   - State clearly that `q_diag` is the innovation covariance diagonal.
   - State clearly that `stationary_factor_std` is the unconditional factor standard deviation.

---

## Tests to add

Add focused tests. Do not rely only on end-to-end panel generation.

### 1. `q_diag` computation test

Given:

```python
a_diag = np.array([0.85, 0.65, 0.65])
stationary_factor_std = np.array([0.0007, 0.0025, 0.0008])
```

the computed `q_diag` must be approximately:

```python
np.array([1.35975e-7, 3.609375e-6, 3.696e-7])
```

### 2. `q_diag` consistency validation test

If both fields are supplied and inconsistent, validation must fail with a clear error message.

Example inconsistent config:

```python
q_diag = np.array([0.0007, 0.0025, 0.0008])
```

This should fail because those are stationary standard deviations, not innovation variances.

### 3. Innovation draw scale test

Mock or seed the RNG and verify that the code uses `sqrt(q_diag)` as the innovation standard deviation.

This can be tested indirectly by simulating a long factor path and checking that the empirical stationary standard deviations are close to `stationary_factor_std`.

### 4. Residual scale matching test

For the configured 5-by-3 grid:

```python
m_grid = [-0.15, -0.075, 0.0, 0.075, 0.15]
tau_grid = [1/12, 1/4, 1/2]
```

and Design A scale:

```python
alpha_0 = 0.0015
alpha_m = 2.0
alpha_tau = 0.05
tau_min = 1/252 * 5
```

verify for each grid point:

\[
\sqrt{
(b^\top\Omega_f b) + s_{\sigma,u}^2
}
\approx
s_{\sigma}(m,\tau).
\]

### 5. Backwards compatibility test, if applicable

If legacy configs are still supported, add one test documenting the behavior. If not supported, add one test confirming that deprecated ambiguous configs fail with a useful error message.

---

## Acceptance criteria

- The persistent-factor code no longer has ambiguous `q_diag` semantics.
- The canonical config uses `stationary_factor_std` and `q_diag` consistently.
- Innovation draws use `sqrt(q_diag)`.
- Design C residual noise is calibrated so the unconditional pointwise scale matches Design A/B.
- Tests cover the covariance/standard-deviation distinction.
- Existing clean panel and noisy panel generation tests continue to pass.
- No changes are made to IS, PF, AWS Batch, Docker, or Terraform.
