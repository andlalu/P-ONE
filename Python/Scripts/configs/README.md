# Experiment configuration

`heston_experiment_run_002.json` is the 500-sample production configuration for
path simulation, clean option panels, four noise scenarios and first-step
IS-CGMM estimation. `heston_experiment_run_001.json` remains the frozen
historical 100-sample configuration.

## Production runners

The single-sample runner remains the numerical entry point for local debugging:

```bash
python Python/Scripts/run_heston_sample.py \
  --config Python/Scripts/configs/heston_experiment_run_002.json \
  --sample-id 0 \
  --output-root outputs/run_002
```

The range runner is a thin orchestration layer around the same `run_sample(...)`
implementation:

```bash
python Python/Scripts/run_heston_samples.py \
  --config Python/Scripts/configs/heston_experiment_run_002.json \
  --output-root /data/p-one/outputs/run_002 \
  --sample-start 0 \
  --sample-end 8 \
  --sample-workers 8 \
  --scenario clean \
  --resume \
  --s3-uri s3://<bucket>/p-one/run_002
```

`sample-start` is inclusive and `sample-end` is exclusive, so `[0, 8)` owns
samples `000` through `007`. Processes run in parallel across samples only;
each process owns one complete sample directory, while scenarios remain
sequential within that sample. `--generation-only` still creates and validates
all five panel variants. `--scenario clean` reuses or creates those shared
panels and estimates only the clean scenario.

The run root contains one authoritative `run.json`. Each `sample_NNN`
directory contains only `path.npz`, `panels.parquet`, `record.json` and
`sample.log`. The Parquet file stores `clean`, `low_iid`, `spatial_corr`,
`persistent_factor` and `variance_linked_factor` rows together; `record.json`
stores embedded validation and all five first-step estimates.

Use `--resume` to verify recorded hashes and continue the first incomplete
stage. Use `--overwrite` to replace only the requested sample directory.

On EC2, persistent EBS is the active filesystem and S3 is an optional durable
replica. The parent range process uploads `run.json` and synchronises each
sample only after its worker returns. It never mounts S3, archives directories
or uses `aws s3 sync --delete`. `--restore-from-s3` requires `--resume` and
restores `run.json` plus only the requested sample range before workers launch.
Use an EBS volume configured with `DeleteOnTermination=false`.

The EC2 wrapper builds `lets_be_rational` once and then invokes the range
runner. For the first eight-sample clean pilot:

```bash
SAMPLE_START=0 \
SAMPLE_END=8 \
SAMPLE_WORKERS=8 \
SCENARIO=clean \
CONFIG=Python/Scripts/configs/heston_experiment_run_002.json \
OUTPUT_ROOT=/data/p-one/outputs/run_002 \
S3_URI=s3://<bucket>/p-one/run_002 \
bash scripts/run_heston_samples_ec2.sh
```

Set `GENERATION_ONLY=1` for generation and validation, or
`RESTORE_FROM_S3=1` to restore the requested range before resuming. The
instance IAM role supplies AWS credentials.

The full sample-000 production integration test is deliberately excluded from
the ordinary fast suite. Run it explicitly with:

```bash
PYTHONPATH=Python python -m pytest \
  Python/Scripts/tests/test_sample_000_clean_production.py \
  -m production_integration \
  -s \
  -vv
```

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
- `noise` contains common IV-floor/price-bound controls and the four noisy
  scenarios.
- `estimation` contains implied-state, quadrature, C-GMM and Powell settings.

## Estimation settings

`implied_state` sets the admissible variance interval, scalar-solver tolerance
and iteration limit, boundary flag, warm-start window, primary solver and
fallback, finite-difference steps and minimum Black vega. Production uses
`bounded_brent` with `golden_section` fallback.

The fixed production COS widths remain `(0.75, 1.25, 2.0)`, with 576 terms
for both generation and estimation.

`quadrature` gives the Gauss-Hermite dimension, order and node scale. The
scale may be one positive number for every coordinate or one positive number
per coordinate, ordered as return then variance. Production uses the even
order-four rule, which has no `(0, 0)` tensor node, with coordinate scales
`(1.0, 4.0)`. `cgmm` sets instrument frequency scales, transition-CF method
and RK4 resolution, optional fixed spacing, and the equal-spacing tolerance.

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
- `variance_linked_factor` is Design C plus a common factor linked to the
  standardised latent variance state; `latent_variance_share` sets its target
  variance share.

`base_seed` fixes scenario seeds. `sigma_min` and `price_epsilon` control the
final noisy-IV floor and strict no-arbitrage price bounds.
