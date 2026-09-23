# P-ONE local production generation specification — run_002, 500 samples, five panel types

## Objective

Execute the complete **local data-generation stage** for the new P-ONE Monte Carlo design and validate it fully in the same task.

The desired final local dataset consists of **500 independent samples**, using the project’s existing zero-based sample convention:

```text
sample IDs 000 through 499
sample range [0, 500)
```

For every sample, regenerate everything from scratch:

1. Heston physical-measure path;
2. clean option panel;
3. Design A — `low_iid`;
4. Design B — `spatial_corr`;
5. Design C — `persistent_factor`;
6. Design D — `variance_linked_factor`, i.e. temporally persistent noise linked to latent variance.

Each sample’s combined panel must therefore contain **five scenario/panel types**:

```text
clean
low_iid
spatial_corr
persistent_factor
variance_linked_factor
```

Do not reuse any `run_001` paths, panels, records, estimates, or other sample artefacts. The point of this task is to avoid cross-run migration complexity and produce one homogeneous fresh `run_002` dataset.

This task includes:
- preflight checks;
- focused tests;
- a short production pilot;
- generation of all 500 samples locally;
- an independent post-generation validation/audit of all 500 samples;
- generation of a concise validation summary.

This task does **not** include:
- S3/AWS upload;
- EC2 launch;
- IS-CGMM estimation;
- report-table generation;
- report rewriting.

After successful validation, stop. The next task will upload the validated dataset to S3 and launch the AWS estimation run.

---

# 1. Use the current local repository state as authoritative

The repository may contain Design-D / no-rounding changes made in the immediately preceding Codex task that are not yet visible on the remote `master` branch.

Therefore:

1. inspect the **current local checkout**;
2. do not replace it with remote `master`;
3. do not reset or discard approved local changes;
4. inspect the current diff/status before running production.

The production generation must use the currently approved implementation containing:
- Design D;
- no tick rounding;
- five total panel types;
- the `run_002` / 500-sample configuration.

If something required by that approved implementation is unexpectedly missing, make only the smallest necessary correction and report it explicitly.

Do not reopen the scientific design or refactor the code.

---

# 2. Reproducibility preflight

Production output must be tied to a stable Git state.

Before generating samples:

```bash
git status --short
git rev-parse HEAD
```

Preferred state: clean working tree at the approved code revision.

Do **not** use `--allow-dirty` for the production generation merely to bypass reproducibility checks.

If the tree is dirty:

- inspect the diff;
- if the only tracked changes are exactly the already-approved Design-D / no-rounding implementation from the immediately preceding task, preserve them and ensure the production run is tied to a stable revision before proceeding;
- do not silently include unrelated changes;
- if unrelated user changes are present, stop rather than modifying or committing them.

Do not commit generated output files.

Record the final Git SHA used for generation in the completion report.

---

# 3. Production configuration

Use the current production configuration for the new run, expected to be:

```text
Python/Scripts/configs/heston_experiment_run_002.json
```

Verify before generation that it resolves to:

```text
run_id       = run_002
n_samples    = 500
base_seed    = 1234500
panel_format = parquet
```

Verify that the noise scenario order is exactly:

```python
(
    "low_iid",
    "spatial_corr",
    "persistent_factor",
    "variance_linked_factor",
)
```

and therefore the full panel order is exactly:

```python
(
    "clean",
    "low_iid",
    "spatial_corr",
    "persistent_factor",
    "variance_linked_factor",
)
```

Verify:
- Design D uses the already-approved `latent_variance_share = 0.5`;
- tick rounding is absent from the active run-002 generation path;
- the existing Heston DGP, simulation grid, option grid, COS basis and seeds remain the approved values;
- `run_001` is not used.

Do not alter numerical design parameters during this task.

---

# 4. Expected production dimensions

For each sample:

```text
weekly dates          = 526
contracts per date    = 15
rows per scenario     = 526 * 15 = 7,890
scenarios per sample  = 5
rows per sample       = 39,450
```

Across all 500 samples:

```text
samples               = 500
rows per scenario     = 3,945,000
total panel rows      = 19,725,000
```

These are hard validation expectations unless the current approved production configuration explicitly differs. If it does, stop and explain the discrepancy rather than silently adjusting the expected design.

---

# 5. Fresh output root

Use the run-002 output root configured by the current production config, expected conceptually to be:

```text
outputs/run_002
```

The goal is a **fresh full regeneration**.

Before deleting or replacing anything:
- resolve the exact configured output path;
- confirm it is a generated `run_002` output directory and not a source directory.

If a previous incomplete/test `run_002` output exists there, remove that run-002 output only so the production generation starts clean.

Do not touch:
- `run_001`;
- source code;
- unrelated output directories;
- S3.

The final expected layout is:

```text
run_002/
    run.json
    sample_000/
        path.npz
        panels.parquet
        record.json
        sample.log
    ...
    sample_499/
        path.npz
        panels.parquet
        record.json
        sample.log
```

No estimation files/results are required beyond the normal empty/pending estimation structure in `record.json`.

---

# 6. Local environment and resource preflight

Before starting the full generation:

1. activate/use the project’s existing Python environment;
2. ensure `PYTHONPATH=Python` or the equivalent project setup;
3. verify required packages are importable;
4. verify the implied-volatility / `lets_be_rational` dependency is built/usable;
5. inspect logical CPU count and available memory;
6. inspect available disk space at the output filesystem.

Do not assume an arbitrary high worker count.

Use the existing range runner’s process-level parallelism, which already constrains numerical-library threads to one.

Choose a conservative worker count from the local machine after the pilot. Do not use so many workers that the machine begins swapping or becomes unstable.

For a typical workstation, a value in roughly the 4–8 range is likely sensible, but detect the actual machine and decide from evidence.

---

# 7. Run focused tests before production

Before any production generation, run the tests that exercise the changed generation path.

At minimum:

```bash
PYTHONPATH=Python python -m pytest \
  Python/OptionData/tests/test_noise.py \
  Python/OptionData/tests/test_panel_io.py \
  Python/Scripts/tests/test_production_configuration.py \
  Python/Scripts/tests/test_sample_run.py \
  Python/Scripts/tests/test_run_heston_samples.py \
  -q
```

Also run any new Design-D-specific test file introduced by the previous task if it is not already included above.

If these tests fail, fix only defects directly related to the approved Design-D/no-rounding production path. Rerun until green.

Do not start the 500-sample run with failing focused tests.

---

# 8. Short production pilot

Use the actual production config and final output root.

Generate a small initial block, e.g. samples `[0,4)`, using `--generation-only`.

This pilot is **part of the production run**, not throwaway data.

Example shape:

```bash
PYTHONPATH=Python python Python/Scripts/run_heston_samples.py \
  --config Python/Scripts/configs/heston_experiment_run_002.json \
  --output-root <resolved-run-002-root> \
  --sample-start 0 \
  --sample-end 4 \
  --sample-workers <conservative-worker-count> \
  --generation-only
```

Do not pass:
- `--scenario`;
- S3 arguments;
- AWS arguments;
- estimation options.

Inspect the pilot before continuing.

For every pilot sample confirm:
- status is generated/validated according to current runner semantics;
- internal validation passed;
- path and panel artefacts exist;
- exactly five scenario values are present;
- each scenario contains 7,890 rows;
- total panel rows = 39,450;
- no estimations were run;
- no tick rounding occurred;
- Design D fields/diagnostics are populated as intended;
- clean fields are exact clean model values.

If the pilot exposes a defect, fix it minimally, delete/restart the fresh run-002 output so all production samples are generated by one consistent code revision/configuration, rerun the focused tests, and repeat the pilot.

Do not mix samples produced under different code SHAs or configuration hashes.

---

# 9. Generate all 500 samples

Once the pilot passes, complete the entire range `[0,500)` using the same Git SHA and config.

Because `[0,4)` already exists and was validated, use the runner’s supported `--resume` semantics so the pilot is reused and the remaining samples are generated.

Conceptually:

```bash
PYTHONPATH=Python python Python/Scripts/run_heston_samples.py \
  --config Python/Scripts/configs/heston_experiment_run_002.json \
  --output-root <resolved-run-002-root> \
  --sample-start 0 \
  --sample-end 500 \
  --sample-workers <chosen-worker-count> \
  --generation-only \
  --resume
```

If the current CLI semantics differ after the previous code task, inspect the runner and use the equivalent supported command. Do not invent a separate generation framework.

Monitor the run to completion.

Do not start estimations.

Do not upload to S3 in this task.

If individual samples fail:
- retain logs;
- identify the cause;
- fix only genuine implementation/runtime issues;
- use the existing resumability mechanism to regenerate/retry incomplete samples;
- do not silently exclude failed samples.

The goal is **500/500 valid samples**.

---

# 10. First validation layer: built-in per-sample validation

The production runner already performs sample/path/panel validation. Treat that as the first validation layer.

At completion, inspect all 500 `record.json` files and require:

- correct `run_id`;
- correct sample ID;
- common configuration hash;
- common Git SHA;
- expected path seed;
- expected A/B/C/D scenario seeds;
- path validation passed;
- combined-panel validation passed;
- five-scenario order valid;
- per-scenario validation passed;
- no recorded hard generation error;
- generation status complete according to current runner semantics.

Do not accept 499/500.

---

# 11. Second validation layer: independent full-run audit

Perform a separate **read-only audit over all 500 completed outputs** in this same task.

This audit may be implemented as:
- a temporary Python script outside tracked source; or
- a small explicit validation utility if the current repository already has an appropriate place for it.

Do not add unnecessary permanent framework code merely for this audit.

Write the audit results into the output root, preferably:

```text
generation_validation.json
generation_validation.txt
```

These are run artefacts, not source files.

The audit must cover all 500 samples, not a random subset.

## 11.1 Directory and artefact integrity

Require:
- exactly sample directories `sample_000` through `sample_499`;
- no missing IDs;
- every sample contains:
  - `path.npz`
  - `panels.parquet`
  - `record.json`
  - `sample.log`
- no temporary `.tmp` artefacts left behind;
- `run.json` exists and matches the current config/hash/SHA.

Where `record.json` stores SHA-256 checksums for path/panel artefacts, recompute and verify them.

---

## 11.2 Panel dimensions and order

For every sample:
- scenario order exactly:
  - clean
  - low_iid
  - spatial_corr
  - persistent_factor
  - variance_linked_factor
- exactly 7,890 rows per scenario;
- exactly 39,450 rows total;
- exactly 526 dates per scenario;
- exactly 15 contracts per date;
- deterministic contract ordering is preserved.

Across the whole run require:
- 500 samples;
- 3,945,000 rows per scenario;
- 19,725,000 rows in total.

---

## 11.3 Parquet metadata

For every `panels.parquet`, inspect metadata and require consistency for:
- run ID;
- sample ID;
- configuration hash;
- Git SHA;
- format version;
- scenario order;
- COS-basis metadata.

All samples must be produced under one common configuration hash and Git SHA.

---

## 11.4 Path sanity

Across all samples verify:
- finite path values;
- non-negative variance;
- correct weekly/daily lengths;
- correct deterministic seed `base_seed + sample_id`;
- no grossly impossible storage/read failures.

Aggregate and report:
- minimum/maximum spot over all paths;
- minimum/maximum variance over all paths;
- distribution/summary of sample-level mean variance;
- distribution/summary of sample-level variance standard deviation.

These are diagnostics, not additional DGP calibration tests.

---

## 11.5 Clean panel exactness

For every clean scenario:
- `estimation_price == model_price` exactly as intended by current implementation;
- `estimation_iv == model_iv` exactly;
- no noise seed;
- no observation-error factor diagnostics populated;
- no capping.

Report any mismatch as a validation failure.

---

## 11.6 Noisy-panel basic validity

For A/B/C/D:
- all relevant price/IV fields finite;
- observed/estimation IV non-negative and above the configured floor as required;
- observed prices lie inside the enforced no-arbitrage interval;
- cap flags and directions are internally consistent;
- deterministic noise seed is correct;
- scenario-specific required factor fields are present;
- fields that should be nullable for that scenario are null.

Aggregate and report, by scenario:
- number of rows;
- mean signed raw IV error;
- MAE;
- RMSE;
- median absolute error;
- 95% absolute-error quantile;
- maximum absolute error;
- lower-bound cap count;
- upper-bound cap count;
- total cap count.

Use **raw IV error before any price projection/inversion effect** where the current schema exposes it, so the scenario noise calibration is assessed directly.

---

# 12. Explicit no-tick-rounding validation

The new production design must have no tick rounding.

Validate this in the actual produced data/code path, not merely by trusting the config.

At minimum:
- verify active run-002 config has no operative tick-rounding setting;
- verify active generation helper does not round raw noisy prices to a price grid;
- verify the format-2/new-schema panel contains the new raw-price representation expected by the approved code;
- verify there is no active `rounding_consistency` requirement;
- demonstrate on generated noisy rows that raw noisy prices are not mechanically constrained to cent/tick multiples.

This is a sanity check only; do not deliberately perturb data.

The output validation report must state clearly:

```text
tick rounding applied: no
no-arbitrage projection/capping retained: yes
```

---

# 13. Design-D validation

Validate Design D over the complete generated run.

## 13.1 Configuration/theoretical calibration

Using the production grid and configuration:
- recompute the Design-D common loading `gamma_V`;
- verify `latent_variance_share = 0.5`;
- verify the equal-weight mean across the 15 unique contract locations of

\[
\gamma_V^2/s_j^2
\]

equals 0.5 to tight numerical tolerance;
- verify the Design-D residual variance is non-negative at every contract location;
- verify the total marginal target identity at every contract location.

Report:
- `gamma_V` in volatility units;
- `gamma_V` in IV basis points;
- min/mean/max theoretical latent-variance variance share across the 15 contract locations;
- min/mean/max theoretical total-error correlation with latent variance implied by the construction.

The expected common loading under the approved design is approximately 13.55 IV bp, but use the actual current configuration as the source of truth.

## 13.2 Generated-factor identity

For every sample/date, or equivalently over all unique date-level D-factor rows, verify:

\[
f_{V,i}
=
\gamma_V
\frac{V_i-\bar V}
{\sqrt{\bar V\sigma_v^2/(2\kappa)}}
\]

to tight numerical tolerance.

Require the factor to be constant across the 15 contracts at a date.

## 13.3 Empirical sanity diagnostics

Across the complete run, report empirical diagnostics for Design D:
- correlation between the variance-linked factor and \(V_i\);
- correlation between total raw Design-D IV error and \(V_i\), preferably by contract location and then min/mean/max across locations;
- empirical variance contribution of the variance-linked component relative to the configured total marginal variance.

Do not require finite-sample empirical quantities to equal theoretical targets exactly. They are sanity diagnostics.

A materially wrong sign or a gross departure from the intended state linkage should be treated as a failure requiring investigation.

---

# 14. Cross-design calibration check

The purpose of A/B/C/D is to vary dependence structure while keeping the overall marginal noise magnitude broadly comparable.

Using all 500 samples, report the raw-IV-error scale by scenario and by contract location.

Check that Design D has not mechanically become a substantially higher-noise design despite the 50% state-linked variance share.

Do not demand identical realised sample RMSEs across independently generated designs. Assess consistency with the common target scale.

For A/B/C/D, produce a compact aggregate table with:
- mean raw error;
- raw-error SD/RMSE;
- MAE;
- 95% absolute quantile.

Optionally include min/max across contract locations if useful.

---

# 15. Determinism spot check

After the full run is complete, perform a non-destructive determinism check on a very small number of samples, e.g. IDs:

```text
000
249
499
```

Regenerate them into a temporary directory using the same config/SHA and compare:
- path arrays;
- scenario rows;
- noise draws;
- Design-D factor values

against the production artefacts.

Require exact equality where the existing deterministic implementation guarantees exact equality.

Delete the temporary validation output afterwards.

This is a validation-only operation, not a replacement of production samples.

---

# 16. Run-level summary artefacts

Create two concise run artefacts in the run-002 output root:

```text
generation_validation.json
generation_validation.txt
```

They should include at least:

```text
run_id
git_sha
configuration_hash
sample_count
sample_id_min
sample_id_max
scenario_order
rows_per_scenario
rows_per_sample
total_rows
samples_with_validation_passed
samples_with_validation_failed
missing_sample_ids
checksum_failures
tick_rounding_applied
total cap counts by scenario
noise summary by scenario
Design-D gamma_V
Design-D theoretical variance-share min/mean/max
Design-D theoretical correlation min/mean/max
Design-D empirical state-error correlation min/mean/max
determinism spot-check result
overall_passed
```

`overall_passed` must only be true if all hard requirements pass.

Do not create report-ready dissertation figures/tables in this task.

---

# 17. Optional second built-in resume validation

After the independent audit passes, it is acceptable and useful to invoke the existing range runner once more over `[0,500)` with:

```text
--generation-only --resume
```

to confirm that all completed samples are recognised as valid/resumable and no regeneration is attempted.

This is a validation step within the same task, not a separate user prompt.

If it is computationally wasteful under the exact current runner semantics, skip it and explain why the independent audit is sufficient.

---

# 18. Failure policy

Do not declare success unless all 500 samples are valid.

If validation finds:
- missing sample;
- corrupt Parquet/path;
- wrong scenario set/order;
- wrong configuration hash/SHA;
- failed built-in validation;
- Design-D calibration/identity failure;
- evidence of tick rounding;
- checksum mismatch;

then:
1. identify the affected sample(s);
2. correct the underlying generation/runtime issue if required;
3. regenerate only affected samples using the same committed code/config;
4. rerun the complete audit;
5. finish only when `overall_passed = true`.

Do not silently remove problematic samples from the Monte Carlo population.

---

# 19. Do not upload or estimate yet

This task must stop after local generation and validation.

Do **not**:
- run `aws s3 sync`;
- call the project’s S3 persistence functions against the production bucket;
- start EC2;
- start IS-CGMM estimation locally;
- generate estimator outputs.

The next Codex task will take the validated local run and:
1. upload it to a new S3 production prefix;
2. provision/restore AWS compute;
3. execute and babysit the 500 × 5 estimation workload.

---

# 20. Final response required from Codex

When everything is complete, provide a concise but complete production report containing:

1. Git SHA used;
2. configuration file/path and configuration hash;
3. resolved local output root;
4. worker count and approximate generation runtime;
5. confirmation that samples are IDs 000–499;
6. confirmation of exactly 500 valid samples;
7. confirmation of five scenario types in the required order;
8. exact total row count;
9. focused-test results;
10. built-in validation result;
11. independent audit result;
12. checksum/integrity result;
13. no-tick-rounding confirmation;
14. no-arbitrage-capping confirmation;
15. noise summary for A/B/C/D;
16. Design-D `gamma_V` and theoretical variance-share/correlation diagnostics;
17. empirical Design-D state-link diagnostics;
18. determinism spot-check result;
19. paths to `generation_validation.json` and `generation_validation.txt`;
20. confirmation that no S3 upload, AWS launch, or estimation was performed.

If any requirement failed, say so explicitly rather than presenting the run as production-ready.
