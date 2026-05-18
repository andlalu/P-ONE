# Codex task — clean DGP path and option-panel generation

## Goal

Implement the first production-quality clean data-generation stage:

1. Generate **100 Heston paths** under the proposed `P`-measure DGP.
2. Store **daily path values** and weekly retained observations.
3. For each path, generate **526 weekly clean option panels**.
4. Store for each option point:
   - option type
   - strike
   - forward
   - maturity
   - clean Heston model price
   - clean Black/Black-Scholes implied volatility
5. Make the panel generation routine **stand-alone and reusable**, because later tasks will contaminate prices with noise and then re-imply observed IVs.
6. Keep the workflow compatible with running locally for 1–2 paths and on EC2 for 100 paths using multiprocessing.

Do not implement noisy option panels, IS estimation, or PF estimation in this task. However, design the generated data layout and IV-solver interface so those steps can reuse the clean panel generation and implied-volatility routines.

---

## Existing repo context

The existing simulator lives under `Python/DGPSimilation`.

Use the existing simulation classes and helpers:
- `HestonParamsP`
- `HestonSimConfig`
- `HestonPath`
- `HestonPathSimulator`
- `AndersenQeVarianceDrawer`
- `save_heston_path_npz`
- `load_heston_path_npz`

The current simulation convention is:
- daily internal step: `delta = 1/252`
- weekly observation every `m_week = 5` daily steps
- retained weekly observations: `t_week = 525`
- weekly state arrays therefore have length `526`
- burn-in: `burnin_days = 756`
- default `s0 = 100`

For this task, set `return_daily = true`, because the daily path values should be stored.

---

## Design principles

1. **Separation of concerns**
   - Path simulation stays in `DGPSimulation`.
   - Option pricing stays in `OptionPricin`.
   - Implied-volatility inversion stays in `ImpliedVolatility`.
   - Experiment orchestration stays in `Scripts`.

2. **Reusable clean panel generation**
   - The clean panel generator should be callable from scripts and from future noise-contamination routines.
   - Do not tie it to AWS.
   - Do not hard-code the run configuration inside the numerical modules.

3. **Parallelism**
   - One sample is one independent unit of work.
   - Use multiprocessing across `sample_id`.
   - Do not use threading for CPU-bound generation.
   - Set numerical backend thread environment variables to avoid oversubscription.

4. **Determinism**
   - `seed = base_seed + sample_id`
   - Output filenames depend on `sample_id`.
   - Re-running the same sample with the same config should reproduce identical path outputs.

5. **No fake prices**
   - If no usable Heston COS pricer exists, add a clean pricing interface and fail explicitly with `NotImplementedError`.
   - Do not generate placeholder option prices.

---

## New or updated package layout

Add or update:

```text
Python/
    ImpliedVolatility/
        __init__.py
        black_iv.py
        black_price.py

    OptionPricing/
        __init__.py
        contracts.py
        heston_cos_pricer.py
        clean_panel.py

    Scripts/
        __init__.py
        configs/
            clean_generation_run_001.json
        generation.py
        generation_local.py
        generation_ec2.py
        validate_generation.py

scripts/
    generation_local.sh
    generation_ec2.sh
    sync_run_to_s3.sh
```

If the repo already has equivalent folders, follow the existing structure - the above structure is just a rough outline.

---

## Configuration file

Create:

```text
Python/Scripts/configs/clean_generation_run_001.json
```

with this shape:

```json
{
  "run_id": "run_001",
  "n_samples": 100,
  "base_seed": 1234500,
  "output_root": "outputs/generation/run_001",
  "dgp": {
    "eta": 1.5,
    "kappa": 3.0,
    "vbar": 0.04,
    "sigma_v": 0.4,
    "rho": -0.7,
    "r": 0.02,
    "q": 0.0
  },
  "simulation": {
    "delta": 0.003968253968253968,
    "m_week": 5,
    "t_week": 525,
    "burnin_days": 756,
    "s0": 100.0,
    "return_daily": true
  },
  "q_measure": {
    "eta_v": 0.0
  },
  "panel": {
    "maturities_years": [0.08333333333333333, 0.25, 0.5],
    "log_moneyness": [-0.15, -0.075, 0.0, 0.075, 0.15],
    "option_type_rule": "otm_by_forward_moneyness",
    "atm_option_type": "call",
    "pricing_method": "COS",
    "iv_method": "lets_be_rational",
    "cos": {
      "n_terms": 256,
      "truncation_L": 10
    }
  },
  "parallel": {
    "workers": null
  }
}
```

Parameter values are placeholders for now. Read them from the config; do not hard-code them.

---

## OTM option-panel convention

The generated option panel should use OTM options:

Let

```text
F = S * exp((r - q) * tau)
K = F * exp(log_moneyness)
```

Use:

```text
log_moneyness < 0  -> put
log_moneyness >= 0 -> call
```

For `log_moneyness == 0`, use the `atm_option_type` in the config. Default: `"call"`.

The point is to avoid storing ITM options dominated by intrinsic value.

---

## Clean panel schema

Store one panel file per sample:

```text
outputs/generation/run_001/panels_clean/sample_000.parquet
...
outputs/generation/run_001/panels_clean/sample_099.parquet
```

Preferred format: Parquet.

If Parquet dependencies are not available, add the dependency cleanly if the repo already uses dependency files. Otherwise fall back to CSV only with an explicit warning.

Required columns:

```text
run_id
sample_id
week_index
t
S
logS
V
r
q
maturity_years
expiry_time
forward
log_moneyness
strike
option_type
is_otm
pricing_method
model_price
model_iv
iv_method
```

Expected row count per sample:

```text
526 * n_maturities * n_log_moneyness
```

With the default config:

```text
526 * 3 * 5 = 7,890 rows per sample
```

For 100 samples:

```text
789,000 rows
```

This is small enough for per-sample Parquet files.

---

## Heston COS pricer interface

Create a clean interface, even if an implementation already exists somewhere else in the repo.

Suggested public API:

```python
@dataclass(frozen=True)
class EuropeanOptionContract:
    maturity_years: float
    log_moneyness: float
    strike: float
    option_type: str  # "call" or "put"

@dataclass(frozen=True)
class HestonQParams:
    kappa_q: float
    vbar_q: float
    sigma_v: float
    rho: float
    r: float
    q: float

class HestonCosPricer:
    def price(
        self,
        *,
        S: float,
        V: float,
        tau: float,
        K: float,
        option_type: str,
        params_q: HestonQParams,
    ) -> float:
        ...
```

Use the `P -> Q` mapping consistently with the report:

```text
kappa_q = kappa - eta_v
vbar_q = kappa * vbar / (kappa - eta_v)
```

If `eta_v` is not yet finalized, keep it explicit in the config under `q_measure`. Do not silently assume `eta_v = 0` unless the config explicitly says so.

---

## Implied-volatility solver

Create a reusable implied-volatility module under `Python/iv`.

Use a forward/discounted Black formulation if possible, because the panel is generated in terms of `F`, `K`, `tau`, and discount factor.

Suggested API:

```python
def black76_price(
    *,
    forward: float,
    strike: float,
    tau: float,
    vol: float,
    discount_factor: float,
    option_type: str,
) -> float:
    ...

def implied_vol_black76(
    *,
    price: float,
    forward: float,
    strike: float,
    tau: float,
    discount_factor: float,
    option_type: str,
    on_bounds: str = "raise",
) -> float:
    ...
```

The `implied_vol_black76` implementation should use Peter Jaeckel's "Let's Be Rational" approach if an acceptable dependency is available.

Preferred implementation options:
1. Use an existing dependency wrapping Jaeckel's algorithm, such as `py_vollib` / `lets_be_rational`, if this can be added cleanly to the repo.
2. If a newer dependency such as `fast-vollib` is selected, justify it in a code comment and keep the wrapper interface stable.
3. If no dependency is added, create the wrapper API and raise `NotImplementedError`.

Do not implement a slow generic Brent solver as the main production solver unless explicitly requested. It may be useful as a test fallback, but not as the primary implementation.

### Bounds and noisy-panel compatibility

The IV solver must be immediately reusable when price noise is added later.

For each input price, check no-arbitrage bounds:

Call:
```text
lower = discount_factor * max(forward - strike, 0)
upper = discount_factor * forward
```

Put:
```text
lower = discount_factor * max(strike - forward, 0)
upper = discount_factor * strike
```

Behavior controlled by `on_bounds`:

```text
"raise" -> raise ValueError if outside bounds
"nan"   -> return np.nan if outside bounds
"clip"  -> clip price to [lower + eps, upper - eps] and solve
```

Use `"raise"` for clean panels. Future noisy panels can use `"nan"` or `"clip"` depending on the experiment design.

---

## Path generation

Each sample should be generated as:

```text
seed = base_seed + sample_id
return_daily = true
```

Store:

```text
outputs/generation/run_001/paths/sample_000.npz
...
outputs/generation/run_001/paths/sample_099.npz
```

The saved path must include:
- `t_week`
- `logS_week`
- `V_week`
- `dlogS_week`
- `logS_daily`
- `V_daily`

If the existing `save_heston_path_npz` does not preserve daily arrays correctly, minimally fix it and add a test.

---

## Local run script

Add:

```text
Python/Scripts/generation_local.py
```

Purpose:
- run only 1 or 2 samples
- use 1 or 2 workers
- write to `outputs/generation/local_run`
- validate outputs

CLI:

```bash
PYTHONPATH=Python python -m Scripts.generation_local \
  --config Python/Scripts/configs/clean_generation_run_001.json \
  --n-samples 2 \
  --workers 2 \
  --output-root outputs/generation/local_run
```

It should call the same underlying worker logic used by EC2.

---

## EC2 master script

Add:

```text
Python/Scripts/generation_ec2.py
```

Purpose:
- generate 100 samples on EC2
- use multiprocessing across `sample_id`
- make use of multiple vCPUs without committing to 8, 16, or 32 vCPUs

CLI:

```bash
PYTHONPATH=Python python -m Scripts.generation_ec2 \
  --config Python/Scripts/configs/clean_generation_run_001.json \
  --workers 15
```

Arguments:
- `--config`
- `--workers`
- `--output-root`
- `--sample-start`
- `--sample-end`
- `--skip-existing`
- `--paths-only`
- `--panels-only`

Default worker logic:
```python
workers = min(max(os.cpu_count() - 1, 1), n_samples)
```

Set these environment variables before worker launch:
```text
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

Reason: avoid oversubscription when each worker uses NumPy/SciPy internally.

---

## Per-sample worker

Implement one per-sample function used by both local and EC2 scripts:

```python
def generate_one_sample(sample_id: int, config: GenerationConfig) -> SampleGenerationResult:
    ...
```

Responsibilities:
1. Generate or load path.
2. Save daily + weekly path NPZ.
3. Generate clean panel from weekly path states.
4. Save clean panel Parquet.
5. Return status, elapsed time, paths to artifacts, and error if failed.

Do not let a failure in one sample kill all workers without reporting which sample failed.

---

## Manifests

Write:

```text
outputs/generation/run_001/config/manifest_generation.csv
```

Columns:

```text
run_id
sample_id
seed
path_file
panel_file
status
elapsed_seconds
error
```

The parent process should also write:

```text
outputs/generation/run_001/config/run_config_resolved.json
```

which includes all defaults after resolution.

---

## Validation command

Add:

```text
Python/Scripts/validate_generation.py
```

CLI:

```bash
PYTHONPATH=Python python -m Scripts.validate_generation \
  --run-root outputs/generation/run_001 \
  --expected-samples 100
```

Checks:
1. Exactly expected number of path files.
2. Exactly expected number of panel files.
3. Each path loads.
4. Path weekly arrays have length 526.
5. `dlogS_week` has length 525.
6. Daily arrays exist.
7. No negative variance.
8. No non-finite values.
9. Each panel has expected row count.
10. Model prices are finite and non-negative.
11. Model IVs are finite and positive.
12. Option type follows the OTM rule.

---

## S3 sync script

Add:

```text
scripts/sync_run_to_s3.sh
```

Usage:

```bash
bash scripts/sync_run_to_s3.sh outputs/generation/run_001 s3://<bucket>/p-one/generation/run_001
```

Implementation:
```bash
aws s3 sync "$LOCAL_RUN_ROOT" "$S3_PREFIX"
```

Do not use `--delete` by default.

---

## Shell wrappers

Add:

```text
scripts/generation_local.sh
scripts/generation_ec2.sh
```

These should:
- activate `.venv` if present
- export `PYTHONPATH=$PWD/Python`
- export numerical thread env vars
- call the appropriate Python module
- tee logs to `outputs/.../logs/`

---

## Tests

Add tests using a miniature config:

```text
n_samples = 2
t_week = 4
burnin_days = 2
maturities = [0.25]
log_moneyness = [-0.05, 0.0, 0.05]
workers = 2
```

Tests:
1. local runner completes
2. path files exist
3. daily arrays are saved
4. panel files exist if pricer exists
5. manifest exists
6. validation succeeds
7. `--skip-existing` does not recompute already completed samples
8. IV solver handles clean model prices
9. IV solver bound behavior works for `"raise"`, `"nan"`, and `"clip"`


---

## Acceptance criteria

1. Running the local test script for 2 samples works.
2. Running the EC2 script for 100 samples can use multiple processes.
3. Daily path values are stored.
4. Weekly clean panels are stored.
5. OTM put/call convention is implemented.
6. Clean model prices and model IVs are both stored.
7. The IV solver is exposed through a reusable wrapper suitable for later noisy panels.
8. No noisy-panel logic is implemented yet.
9. No IS/PF estimation logic is implemented yet.
10. No AWS Batch, Step Functions, Docker, Terraform, or CloudFormation is added in this task.
