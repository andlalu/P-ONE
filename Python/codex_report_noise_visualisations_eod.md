# Codex Task 3: EOD Report Update — Noise Panel Visualisations and Text for `Reports/18062026`

## Purpose

Create a focused report update for:

```text
Reports/18062026/draft_proposalJun26.tex
```

The EOD objective is **not** to run IS-CGMM estimation, not to run optimizer tests, and not to build AWS/local-run infrastructure. The objective is to improve the current report by adding:

1. publication-quality visualisations showing how the clean Heston implied-volatility panel is distorted by the three observation-error designs;
2. figure-generation code with BLL-style greyscale aesthetics;
3. LaTeX figure environments, captions, and cross-references;
4. concise explanatory text around the figures;
5. a short computational-status paragraph describing the current IS-CGMM implementation status.

This task is meant to produce a credible report update by end of day.

---

## Current report context

The report file already contains a section:

```text
Reports/18062026/draft_proposalJun26.tex
```

with a subsection:

```latex
\subsubsection{Option Prices with Observation Errors}
```

This section defines:

- observed option-price vectors;
- clean/model prices;
- clean implied volatilities;
- the interpretation of the observed-clean difference as option observation error;
- three implied-volatility-space noise specifications.

The current text has explicit placeholders:

```latex
\emph{Example sampling noise for specification A: [to be added].}
\emph{Example sampling noise for specification B: [to be added].}
\emph{Example sampling noise for specification C: [to be added].}
```

Replace or supplement these placeholders with figure references and explanatory text.

The current observation-error section already states that the noisy panels are generated in implied-volatility space, transformed back to prices through Black-Scholes, rounded/capped, and then inverted back to observed implied volatilities. Preserve this logic.

---

## Non-goals

Do **not** implement or run:

1. IS-CGMM estimation;
2. optimizer work;
3. local full-panel AWS-like estimation runs;
4. AWS/S3/Batch infrastructure;
5. particle filtering;
6. second-stage GMM;
7. new model specifications;
8. new estimation results.

This is a report-update and visualisation task only.

---

## Files to inspect first

Inspect these files before making changes:

```text
Reports/18062026/draft_proposalJun26.tex
Python/Scripts/generation.py
Python/Scripts/configs/clean_generation_run_001.json
Python/OptionPricing/noisy_panel.py
Python/OptionPricing/clean_panel.py
Python/DGPSimulation/heston_simulator.py
Python/Estimation/ISCGMM/panel.py
```

Also inspect any existing plotting or report scripts if present:

```text
Python/Scripts/
Reports/18062026/
Reports/18062026/figures/
Reports/18062026/scripts/
```

If no plotting scripts exist, create a new self-contained script as specified below.

---

# 1. Required script

Create a new self-contained plotting script:

```text
Python/Scripts/make_noise_visualisations.py
```

or, if the repository convention prefers report-local scripts:

```text
Reports/18062026/make_noise_visualisations.py
```

Prefer the first path if existing Python scripts live under `Python/Scripts`.

The script must be callable from the repository root with:

```bash
PYTHONPATH=Python python Python/Scripts/make_noise_visualisations.py \
  --config Python/Scripts/configs/clean_generation_run_001.json \
  --sample-id 0 \
  --date-index 100 \
  --output-dir Reports/18062026/figures/noise \
  --generate-if-missing
```

If `date-index=100` is unavailable, the script should choose a sensible valid interior date and record it in stdout.

The script should:

1. load or generate the clean panel for one sample;
2. load or generate the noisy panels for:
   - `low_iid`;
   - `spatial_corr`;
   - `persistent_factor`;
3. load the persistent-factor time series if available;
4. produce greyscale figures suitable for direct inclusion in the LaTeX report;
5. write all figures to:

```text
Reports/18062026/figures/noise/
```

Use `.pdf` as the main figure format. Optionally also write `.png` for quick preview.

---

# 2. Required figure set

Implement a small but high-value figure set. Avoid overproduction. The goal is to make the report visually convincing, not to create every possible diagnostic.

## Figure 1: Clean and contaminated implied-volatility panels at one date

Suggested filename:

```text
Reports/18062026/figures/noise/noise_surface_comparison_sample000_dateXXX.pdf
```

Purpose:

Show the clean implied-volatility panel and the corresponding observed panels under Designs A, B, and C at one representative date.

Preferred layout:

```text
2 × 2 panel

top-left:     clean IV surface
top-right:    Design A observed IV surface
bottom-left:  Design B observed IV surface
bottom-right: Design C observed IV surface
```

Because each weekly panel has only:

```text
3 maturities × 5 log-forward moneyness levels = 15 contracts
```

avoid glossy 3D surfaces as the main representation. Prefer one of:

1. heatmap with moneyness on x-axis and maturity on y-axis;
2. line plot by moneyness, with one line per maturity;
3. combined heatmap plus overlaid contour.

Recommended for academic print:

- use line plots by maturity, one subplot per design;
- x-axis: log-forward moneyness;
- y-axis: implied volatility;
- three maturity curves:
  - short maturity: black solid;
  - medium maturity: dark grey dashed;
  - long maturity: medium grey dotted;
- show clean surface in light grey in noisy panels as a reference;
- avoid rainbow colours.

This is often more readable than a 3D surface for only 15 grid points.

## Figure 2: Noise perturbation by design at one date

Suggested filename:

```text
Reports/18062026/figures/noise/noise_difference_comparison_sample000_dateXXX.pdf
```

Purpose:

Show the actual implied-volatility perturbation:

```math
Delta sigma_{i,j}^{d}
=
sigma_{i,j}^{obs,d}
-
sigma_{i,j}^{0},
d in {A,B,C}.
```

Preferred layout:

```text
1 × 3 panel

Design A: Delta IV heatmap
Design B: Delta IV heatmap
Design C: Delta IV heatmap
```

Axes:

- x-axis: log-forward moneyness;
- y-axis: maturity;
- cell label or contour optional.

Use a symmetric greyscale scale centered at zero, e.g.:

```text
vlim = max(abs(all_differences))
vmin = -vlim
vmax = +vlim
```

Use a diverging greyscale convention if possible:

- negative perturbation: light grey;
- zero: mid grey / white boundary;
- positive perturbation: dark grey.

If a true diverging greyscale is awkward, use contour lines and signed labels.

## Figure 3: Marginal noise scale and spatial correlation structure

Suggested filename:

```text
Reports/18062026/figures/noise/noise_design_cross_section_structure.pdf
```

Purpose:

Show the deterministic components of the noise design:

1. marginal IV noise scale `s_sigma(m,tau)`;
2. Design B spatial correlation matrix `R_i`.

Preferred layout:

```text
1 × 2 panel

left:  heatmap of s_sigma(m, tau)
right: heatmap of R_i ordered by maturity then moneyness
```

This figure directly supports the formal definitions in the text.

## Figure 4: Persistent-factor time dependence in Design C

Suggested filename:

```text
Reports/18062026/figures/noise/noise_persistent_factor_timeseries_sample000.pdf
```

Purpose:

Show how Design C differs fundamentally from A/B because it introduces serially persistent common components.

Preferred layout:

```text
2 × 1 panel

top:    time series of f_0, f_1, f_2
bottom: implied-volatility perturbation over time for three representative contracts
```

Representative contracts:

- ATM short maturity;
- low-moneyness short/medium maturity;
- high-moneyness long maturity.

Use black/dark grey/medium grey lines, distinct line styles.

If loading the factor time series is difficult, produce only the bottom panel from observed-clean IV differences for selected contracts, and explain in the caption that the plotted series reflect the realized factor-driven perturbation.

---

# 3. Required aesthetics

Use BLL-style greyscale / print-friendly formatting.

General plotting rules:

```python
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})
```

Do not use seaborn.

Use:

```text
black
dim grey
grey
light grey
white
```

Preferred line styles:

```text
solid
dashed
dotted
dash-dot
```

Recommended figure sizes:

```text
single-column: 6.0 × 3.5 inches
full-width:    7.2 × 4.8 inches
```

Save with:

```python
fig.savefig(path, bbox_inches="tight")
```

Do not use colour palettes that will look poor in black-and-white print.

---

# 4. Data assumptions

The script should be robust to both `.parquet` and `.csv`, following the existing loader conventions where possible.

Expected generated panel layout:

```text
outputs/generation/run_001/
  panels_clean/sample_000.parquet
  panels_observed/low_iid/sample_000.parquet
  panels_observed/spatial_corr/sample_000.parquet
  panels_observed/persistent_factor/sample_000.parquet
  noise_factors/persistent_factor/sample_000.csv
```

Do not require all noisy panels to exist at script startup if `--generate-if-missing` is passed. Reuse existing generation functions rather than duplicating generation logic.

The script should use columns such as:

```text
week_index
t
S
logS
maturity_years
log_moneyness
strike
option_type
model_iv or clean_iv
observed_iv
model_price or clean_price
observed_price
r
q
V
```

Be defensive about exact column names. If a required column is absent, raise a clear error.

---

# 5. Computation constraints

This script may generate one sample and noisy panels, but it must **not** run estimation.

Do not import or call:

```text
Estimation.ISCGMM.optimizer
Estimation.ISCGMM.estimate
CgmmFirstStepCriterion
```

This is intentionally independent from estimation.

---

# 6. Required LaTeX updates

Modify:

```text
Reports/18062026/draft_proposalJun26.tex
```

Add or replace text around the three placeholders.

## 6.1 Replace placeholder after Design A

Replace:

```latex
\emph{Example sampling noise for specification A: [to be added].}
```

with a sentence like:

```latex
Figure~\ref{fig:noise_surface_comparison} illustrates the resulting contaminated implied-volatility panels for a representative simulated date. The low i.i.d. design perturbs individual contracts independently while preserving the broad shape of the clean surface.
```

## 6.2 Replace placeholder after Design B

Replace:

```latex
\emph{Example sampling noise for specification B: [to be added].}
```

with a sentence like:

```latex
The cross-sectionally correlated design creates locally coherent distortions across nearby moneyness and maturity points. This is visible both in the contaminated panel in Figure~\ref{fig:noise_surface_comparison} and in the signed perturbation plots in Figure~\ref{fig:noise_difference_comparison}.
```

## 6.3 Replace placeholder after Design C

Replace:

```latex
\emph{Example sampling noise for specification C: [to be added].}
```

with a sentence like:

```latex
Unlike Designs A and B, the persistent-factor design induces serially dependent quote distortions. Figure~\ref{fig:noise_persistent_factor} displays the realized common factors and their effect on representative contracts over time.
```

## 6.4 Add figure environments

Add figure environments near the observation-error subsection. Use exact filenames produced by the script.

### Figure surface comparison

```latex
\begin{figure}[H]
    \centering
    \includegraphics[width=0.95\textwidth]{figures/noise/noise_surface_comparison_sample000_dateXXX.pdf}
    \caption{Clean and contaminated implied-volatility panels for a representative simulated date. The clean panel is generated exactly from the Heston pricing model. Designs A--C introduce increasingly structured quote distortions in implied-volatility space before prices are reconstructed and projected into admissible price intervals.}
    \label{fig:noise_surface_comparison}
\end{figure}
```

### Figure differences

```latex
\begin{figure}[H]
    \centering
    \includegraphics[width=0.95\textwidth]{figures/noise/noise_difference_comparison_sample000_dateXXX.pdf}
    \caption{Signed implied-volatility perturbations relative to the clean model-implied panel for Designs A--C. The plots show observed minus clean implied volatility across log-forward moneyness and maturity.}
    \label{fig:noise_difference_comparison}
\end{figure}
```

### Figure design structure

```latex
\begin{figure}[H]
    \centering
    \includegraphics[width=0.95\textwidth]{figures/noise/noise_design_cross_section_structure.pdf}
    \caption{Cross-sectional structure of the observation-error designs. The left panel shows the marginal implied-volatility error scale. The right panel shows the spatial correlation matrix used in Design B after ordering contracts by maturity and log-forward moneyness.}
    \label{fig:noise_cross_section_structure}
\end{figure}
```

### Figure persistent factor

```latex
\begin{figure}[H]
    \centering
    \includegraphics[width=0.95\textwidth]{figures/noise/noise_persistent_factor_timeseries_sample000.pdf}
    \caption{Temporal dependence in Design C. The upper panel shows the realized level, moneyness-slope and maturity-slope observation-error factors. The lower panel shows the resulting serially persistent implied-volatility perturbations for selected representative contracts.}
    \label{fig:noise_persistent_factor}
\end{figure}
```

Replace `dateXXX` with the actual date index selected by the script.

---

# 7. Add computational status paragraph

Add a short paragraph near the end of the main text, or before the appendix, titled for example:

```latex
\paragraph{Current computational implementation.}
```

Suggested content:

```latex
The simulation engine and option-panel generation code have been implemented for the Heston design. The first-stage implied-state C-GMM estimator has also been implemented in Python. It uses scalar date-by-date implied-variance inversion from option-implied volatilities, a Heston \(\mathbb{P}\)-transition characteristic-function criterion, and a first-step C-GMM objective with analytically integrated Gaussian instruments. The transition characteristic function can be evaluated either by an analytic Heston Riccati solution or by a Runge--Kutta reference solver. The next computational step is to introduce a bounded staged derivative-free optimizer and then use the resulting estimator in the Monte Carlo study.
```

Adjust wording to fit the surrounding report style.

Do not claim that Monte Carlo estimation results are available.

---

# 8. Optional table of figure-generation settings

If useful, add a small LaTeX table or textual note in the appendix listing:

```text
sample_id
date_index
maturities
log-forward moneyness grid
noise parameter values
```

This is optional. Do not spend too much effort if time is limited.

---

# 9. Script diagnostics / manifest

The plotting script should write a small manifest:

```text
Reports/18062026/figures/noise/noise_visualisation_manifest.json
```

Containing:

```json
{
  "sample_id": 0,
  "date_index": 100,
  "scenario_files": {
    "clean": "...",
    "low_iid": "...",
    "spatial_corr": "...",
    "persistent_factor": "..."
  },
  "figures": [
    "...pdf"
  ],
  "created_by": "Python/Scripts/make_noise_visualisations.py"
}
```

This helps keep the report reproducible.

---

# 10. Acceptance criteria

The task is complete when:

1. a self-contained plotting script exists;
2. the script produces all required figures under `Reports/18062026/figures/noise/`;
3. figures are greyscale, print-friendly, and suitable for the report;
4. `draft_proposalJun26.tex` includes the figure environments and captions;
5. the three placeholder sentences are replaced or supplemented;
6. a computational-status paragraph is added;
7. no estimation run is required;
8. no generated outputs are committed except report figures and possibly the manifest if desired;
9. the report compiles or, if compilation is not available, the LaTeX syntax is manually checked.

---

# 11. Suggested commands

Generate figures:

```bash
PYTHONPATH=Python python Python/Scripts/make_noise_visualisations.py \
  --config Python/Scripts/configs/clean_generation_run_001.json \
  --sample-id 0 \
  --date-index 100 \
  --output-dir Reports/18062026/figures/noise \
  --generate-if-missing
```

If the report is compiled manually:

```bash
cd Reports/18062026
pdflatex draft_proposalJun26.tex
bibtex draft_proposalJun26
pdflatex draft_proposalJun26.tex
pdflatex draft_proposalJun26.tex
```

If bibliography compilation is not configured, at least run a single `pdflatex` pass or perform a syntax check.

---

# 12. Final Codex response

At completion, report:

1. files changed;
2. figures created;
3. selected sample/date;
4. whether panels were generated or loaded;
5. where figures were inserted in the LaTeX file;
6. whether the report compiled;
7. any remaining limitations.
