# Chapter 3 — Contaminated option panel generation

## 1. Goal

Extend the clean Chapter 3 option-panel generation workflow to create three contaminated panel versions for every clean sample.

This task is based on Section 1.3.1, **Option Prices with Observation Errors**, in:

```text
Reports/20052026/draft_proposalMay26.tex
```

The LaTeX section is the conceptual source of truth for the three observation-error designs. The implementation should follow the notation and modelling choices in that section.

The immediate deliverable is:

```text
100 clean samples
  -> 100 clean option panels
  -> 3 contaminated panel versions per clean panel
```

The contaminated panel versions are:

```text
Design A: Low I.I.D. noise
Design B: Cross-sectionally correlated noise
Design C: Temporally persistent noise
```

Do **not** implement IS estimation, particle filtering, AWS Batch, Docker, Terraform, Step Functions, or Student-t / mixture / stale-quote noise in this task.

---

## 2. Existing context

The repo already contains the P-measure Heston path simulator under:

```text
Python/sim/
```

The existing clean generation workflow should remain the input to this task. The contaminated-panel module must not re-simulate paths or re-price the Heston model from scratch if clean panels already exist.

Expected clean panel inputs should contain, or be extended to contain, at least:

```text
sample_id
week_index
t
S
logS
V
maturity_years
tau
K / strike
F / forward
log_moneyness
option_type
clean_model_price
clean_model_iv
clean_model_vega
pricing_method
```

If current column names differ, adapt minimally and document the mapping.

---

## 3. Core modelling principle

All three contaminated panel designs are **IV-first**.

For every design, first construct a raw contaminated implied volatility:

```text
tilde_sigma
```

Then map it to a raw contaminated price using Black-Scholes:

```text
tilde_price = BS(S, K, tau, r, q, tilde_sigma, option_type)
```

Then apply market quote mechanics:

```text
raw contaminated IV
  -> Black-Scholes price
  -> tick rounding
  -> no-arbitrage capping / projection
  -> final observed price
  -> final observed IV by Black-Scholes inversion
```

This means the vega approximation:

```text
delta_price ≈ vega * delta_sigma
```

must **not** be used to generate Design A or Design B. Vega should still be stored for diagnostics and possible vega-scaled estimation later.

---

## 4. Notation to implement

For option contract `j` at observation date `i`, use:

```text
m_i_j   = log(K_i_j / F_i_j)
tau_i_j = T_i_j - t_i
```

In code:

```text
log_moneyness
tau
```

The marginal IV error scale is:

```text
s_sigma_i_j = alpha_0 * (
    1
    + alpha_m * abs(log_moneyness_i_j)
    + alpha_tau / sqrt(max(tau_i_j, tau_min))
)
```

The lower bound for raw contaminated implied volatility is:

```text
sigma_min
```

The raw contaminated IV should be capped below at `sigma_min`:

```text
tilde_sigma = max(sigma_min, clean_model_iv + noise)
```

Use normal distributions only.

---

## 5. Design A: Low I.I.D. noise

### Definition

For each contract independently:

```text
tilde_sigma_i_j_A =
    max(
        sigma_min,
        clean_model_iv_i_j + s_sigma_i_j * Z_i_j
    )

Z_i_j ~ iid N(0, 1)
```

### Interpretation

This is the low-noise benchmark. It represents small quote observation errors whose scale varies with moneyness and maturity, without imposing dependence across contracts.

### Implementation details

- Generate independent standard normal shocks per row.
- Use deterministic seeding per `(sample_id, scenario_name)`.
- Write one contaminated panel per sample.

Scenario name:

```text
low_iid
```

---

## 6. Design B: Cross-sectionally correlated noise

### Definition

At each observation date `i`, let:

```text
epsilon_i_B ~ N(0, D_i R_i D_i)
```

where:

```text
D_i = diag(s_sigma_i_1, ..., s_sigma_i_q)
```

and the correlation matrix is:

```text
R_i_jk = exp(
    - abs(m_i_j - m_i_k) / ell_m
    - abs(tau_i_j - tau_i_k) / ell_tau
)
```

Then for each contract:

```text
tilde_sigma_i_j_B =
    max(
        sigma_min,
        clean_model_iv_i_j + epsilon_i_j_B
    )
```

### Interpretation

This design allows nearby strikes and maturities to share quote or surface construction errors. It is intended to reflect spatial dependence in the option observation error, following the idea in the report text.

### Implementation details

- Group each clean panel by `week_index`.
- For each weekly panel, construct `D_i` and `R_i`.
- Draw one multivariate normal vector per `week_index`.
- The dimension is the number of contracts at that week.
- Add a small diagonal jitter to `R_i` if needed for numerical stability:

```text
R_i += correlation_jitter * I
```

- Use Cholesky if possible; if Cholesky fails, increase the jitter up to a configured maximum and record that the fallback was used.

Scenario name:

```text
spatial_corr
```

---

## 7. Design C: Temporally persistent noise

### Definition

Let:

```text
f_i = A f_{i-1} + xi_i
xi_i ~ N(0, Q)
```

where:

```text
f_i = (f_0_i, f_1_i, f_2_i)'
```

The matrix `A` is diagonal with entries inside the unit circle:

```text
A = diag(a_0, a_1, a_2)
```

The matrix `Q` is the covariance matrix of the factor innovations. In the first implementation, support diagonal `Q` via factor innovation standard deviations:

```text
q_diag = [q_0, q_1, q_2]
Q = diag(q_0^2, q_1^2, q_2^2)
```

For each contract, define the loading:

```text
b_i_j = (1, m_i_j, tau_i_j)'
```

and generate:

```text
tilde_sigma_i_j_C =
    max(
        sigma_min,
        clean_model_iv_i_j + b_i_j' f_i + nu_i_j
    )

nu_i_j ~ N(0, s_sigma_u_i_j^2)
```

For the residual scale, use the same marginal scale function with a separate multiplier:

```text
s_sigma_u_i_j = residual_scale_multiplier * s_sigma_i_j
```

### Interpretation

This design introduces a persistent common component in the observation error. The common component affects the level, moneyness slope and maturity slope of the contaminated implied-volatility panel. It is not a stochastic process for the clean surface; it is a stochastic process for the observation error.

### Implementation details

- Draw the factor process once per sample and scenario, across all `week_index` values.
- Sort the panel by `week_index` before generating factors.
- Initialize `f_0` either at zero or from the stationary distribution. Default to zero for reproducibility and simplicity.
- Add independent residual noise per contract using `s_sigma_u_i_j`.
- Use deterministic seeding per `(sample_id, scenario_name)`.
- Store factor draws by week in a sidecar file, e.g.:

```text
noise_factors/persistent_factor/sample_000.csv
```

Scenario name:

```text
persistent_factor
```

---

## 8. Black-Scholes price mapping

For all scenarios, after generating `tilde_sigma`, compute:

```text
tilde_model_price = black_scholes_price(
    S,
    K,
    tau,
    r,
    q,
    tilde_sigma,
    option_type
)
```

Requirements:

- Support calls and puts.
- The current clean panel uses OTM put and call wings; preserve `option_type`.
- The BS function must handle small maturities and small volatilities robustly.
- If `tau <= 0`, fail clearly unless the clean panel explicitly contains expiry-date contracts. The expected panel has strictly positive maturities.

---

## 9. Tick rounding and no-arbitrage capping

After computing the raw contaminated price, apply tick rounding:

```text
rounded_price = tick_size * round(tilde_model_price / tick_size)
```

Then project / cap into no-arbitrage bounds.

For a call:

```text
lower = max(S * exp(-q * tau) - K * exp(-r * tau), 0)
upper = S * exp(-q * tau)
```

For a put:

```text
lower = max(K * exp(-r * tau) - S * exp(-q * tau), 0)
upper = K * exp(-r * tau)
```

The final observed price is:

```text
observed_price = min(
    upper - price_epsilon,
    max(lower + price_epsilon, rounded_price)
)
```

Do not re-sample invalid draws. Use capping/projection. Record capping indicators.

Required output flags:

```text
was_price_capped
cap_direction      # "none", "lower", "upper"
raw_price_before_rounding
price_after_rounding
observed_price
```

---

## 10. Implied-vol inversion

After obtaining `observed_price`, compute:

```text
observed_iv = implied_volatility(
    observed_price,
    S,
    K,
    tau,
    r,
    q,
    option_type
)
```

Requirements:

- Prefer an implementation of Jaeckel's "Let's Be Rational" if already available or easy to add.
- If using an external package such as `py_vollib` / `lets_be_rational`, isolate it behind a local wrapper.
- If no robust dependency is available, implement a fallback Brent solver behind the same interface.
- The rest of the code must call only the local wrapper, not the third-party package directly.

Suggested wrapper location:

```text
Python/pricing/implied_vol.py
```

Suggested interface:

```python
def implied_volatility_from_price(
    price: float,
    S: float,
    K: float,
    tau: float,
    r: float,
    q: float,
    option_type: str,
) -> float:
    ...
```

Also expose vectorized helpers if useful, but keep the scalar function correct and tested first.

---

## 11. Output schema

For each contaminated panel row, store all original clean-panel columns plus:

```text
noise_scenario              # low_iid / spatial_corr / persistent_factor
raw_contaminated_iv          # tilde_sigma
raw_price_before_rounding
price_after_rounding
observed_price
observed_iv
was_price_capped
cap_direction
noise_draw                   # for A: scalar shock; for B/C can be contract-level realized epsilon
noise_seed
```

For Design C, also store or write separately:

```text
factor_0
factor_1
factor_2
```

per `week_index`, either as repeated columns or a sidecar file. Prefer a sidecar factor file to avoid repeating factors across rows.

Recommended output layout:

```text
outputs/ch3/run_001/
    panels_clean/
        sample_000.parquet
        ...

    panels_observed/
        low_iid/
            sample_000.parquet
            ...
        spatial_corr/
            sample_000.parquet
            ...
        persistent_factor/
            sample_000.parquet
            ...

    noise_factors/
        persistent_factor/
            sample_000.csv
            ...

    config/
        noise_scenarios.json
        manifest_noisy_panels.csv
```

---

## 12. Configuration

Add a config file, or extend the existing run config, with:

```json
{
  "noise": {
    "base_seed": 9900000,
    "sigma_min": 0.0001,
    "price_epsilon": 1e-10,
    "tick_size": 0.01,
    "scenarios": {
      "low_iid": {
        "enabled": true,
        "alpha_0": 0.0025,
        "alpha_m": 1.0,
        "alpha_tau": 0.05,
        "tau_min": 0.01984126984126984
      },
      "spatial_corr": {
        "enabled": true,
        "alpha_0": 0.0050,
        "alpha_m": 1.0,
        "alpha_tau": 0.05,
        "tau_min": 0.01984126984126984,
        "ell_m": 0.10,
        "ell_tau": 0.25,
        "correlation_jitter": 1e-10,
        "max_correlation_jitter": 1e-6
      },
      "persistent_factor": {
        "enabled": true,
        "alpha_0": 0.0030,
        "alpha_m": 1.0,
        "alpha_tau": 0.05,
        "tau_min": 0.01984126984126984,
        "a_diag": [0.85, 0.80, 0.75],
        "q_diag": [0.0015, 0.0010, 0.0010],
        "residual_scale_multiplier": 0.50,
        "factor_initialization": "zero"
      }
    }
  }
}
```

The numerical values are placeholders. They should be kept configurable and not hard-coded.

---

## 13. CLI

Add a command-line entry point:

```bash
PYTHONPATH=Python python -m experiments.generate_ch3_noisy_panels   --config Python/experiments/configs/ch3_run_001.json   --run-root outputs/ch3/run_001   --workers 8
```

Requirements:

- Use multiprocessing across `sample_id`.
- Do not use threads for CPU-bound work.
- Use worker count from CLI.
- Support:

```text
--sample-start
--sample-end
--scenarios low_iid,spatial_corr,persistent_factor
--skip-existing
```

The script should process clean panel files from:

```text
<run-root>/panels_clean/sample_XXX.parquet
```

and write contaminated panels under:

```text
<run-root>/panels_observed/<scenario>/sample_XXX.parquet
```

---

## 14. Manifest

Write:

```text
<run-root>/config/manifest_noisy_panels.csv
```

Columns:

```text
sample_id
noise_scenario
seed
input_clean_panel
output_observed_panel
factor_file
status
elapsed_seconds
n_rows
n_capped_lower
n_capped_upper
n_capped_total
error_message
```

---

## 15. Validation command

Add or extend:

```bash
PYTHONPATH=Python python -m experiments.validate_ch3_noisy_panels   --run-root outputs/ch3/run_001
```

Checks:

- For every clean sample and enabled scenario, there is an observed panel file.
- Observed row count equals clean row count.
- Required columns exist.
- `observed_price` is finite.
- `observed_iv` is finite.
- `observed_iv >= sigma_min`.
- Price bounds are satisfied.
- Cap flags are consistent with observed prices.
- For Design B, generated covariance matrices are numerically valid.
- For Design C, factor sidecar files exist and contain the expected number of weeks.

---

## 16. Tests

Add tests with a miniature panel:

- 2 samples
- 2 weeks
- 2 maturities
- 3 moneyness points
- both puts and calls if available

Tests:

1. Design A creates output and is deterministic.
2. Design B creates output and has non-zero cross-sectional covariance.
3. Design C creates output and factor files.
4. Tick rounding works.
5. Capping works and records flags.
6. IV inversion works after capping.
7. `--skip-existing` avoids recomputation.
8. Validation passes.

---

## 17. Non-goals

Do not implement:

- Student-t noise
- mixture noise
- stale quote logic
- price-space vega approximation noise
- IS estimation
- PF estimation
- AWS Batch
- Docker
- Terraform
- Step Functions

---

## 18. Acceptance criteria

The task is complete when:

- The repo can generate contaminated panel versions A/B/C from existing clean panels.
- The implementation is deterministic by `(sample_id, scenario_name, base_seed)`.
- The panel generator is stand-alone and reusable.
- The same IV solver wrapper is used for clean and contaminated panels where applicable.
- The output format is compatible with later IS and PF estimation.
- Tests pass for the miniature run.
- Existing simulator tests continue to pass.
