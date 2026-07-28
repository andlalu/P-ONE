# Experiment configuration

`heston_experiment_run_001.json` is the production configuration for path
simulation, clean option panels, noise scenarios and first-step IS-CGMM
estimation.

## Production runner

One process owns one complete sample locally and on AWS:

```bash
python Python/Scripts/run_heston_sample.py \
  --config Python/Scripts/configs/heston_experiment_run_001.json \
  --sample-id 0 \
  --output-root outputs/run_001
```

The run root contains one authoritative `run.json`. Each `sample_NNN`
directory contains only `path.npz`, `panels.parquet`, `record.json` and
`sample.log`. The Parquet file stores `clean`, `low_iid`, `spatial_corr` and
`persistent_factor` rows together; `record.json` stores embedded validation
and all four first-step estimates.

Use `--resume` to verify recorded hashes and continue the first incomplete
stage. Use `--overwrite` to replace only the requested sample directory.
AWS wrappers should keep S3 transfer outside the Python runner and use the
array index as `--sample-id`.

## Top-level sections

- `run` sets the run identifier, sample count, base seed, output directory,
  panel format and worker count.
- `dgp` contains the physical-measure Heston parameters.
- `simulation` sets the daily step, weekly sampling interval, number of weeks,
  burn-in, initial spot and daily-path retention.
- `q_measure.eta_v` maps physical variance dynamics to risk-neutral pricing
  dynamics.
- `panel` sets the log-moneyness grid and the option type used exactly at the
  money.
- `cos` fixes the maturity grid, calibrated effective widths and term counts.
  `generation_n_cos` prices generated panels; `estimation_n_cos` prices inside
  state inversion and may legitimately differ.
- `noise` contains common rounding/bound controls and the three scenarios.
- `estimation` contains implied-state, quadrature, C-GMM and Powell settings.

## Estimation settings

`implied_state` sets the admissible variance interval, scalar-solver tolerance
and iteration limit, boundary flag, warm-start window, primary solver and
fallback, finite-difference steps and minimum Black vega. Production uses
`bounded_brent` with `golden_section` fallback.

`quadrature` gives the Gauss-Hermite dimension, order and node scale. `cgmm`
sets instrument frequency scales, transition-CF method and RK4 resolution,
optional fixed spacing, and the equal-spacing tolerance.

`optimizer` screens the central Heston start and at most two relative
perturbations within the natural parameter bounds. It then runs:

1. `coarse_pass`, a broad Powell search;
2. `refinement_pass`, a tighter Powell search starting from the coarse result.

Both passes minimise the same first-step C-GMM criterion. They are optimisation
passes, not first- and second-stage C-GMM estimators. `penalty_value` is the
finite value returned after an expected numerical failure, and
`progress_every` controls logging frequency.

## Noise scenarios

All scenarios use `alpha_0`, `alpha_m`, `alpha_tau` and `tau_min` for their
marginal scale.

- `low_iid` adds independent Gaussian IV noise.
- `spatial_corr` adds within-date correlation using moneyness and maturity
  length scales, with bounded Cholesky jitter and an eigenvalue fallback.
- `persistent_factor` uses three persistent factor coefficients, stationary
  factor standard deviations and the equivalent innovation variances. Its
  residual policy preserves the configured total marginal scale.

`base_seed` fixes scenario seeds. `sigma_min`, `price_epsilon` and `tick_size`
control the final noisy-IV floor, strict price bounds and price rounding.
