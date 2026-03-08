# Chapter 3 — Heston path generation (P-measure)  
**Daily simulation / Weekly observations (10 years = 525 weeklies)**

## 1. Scope (this file)
This document specifies **only** the simulation of Heston state paths under **P**:
- state: \(X_t = [\log(S_t), V_t]\)
- internal step: daily \(\delta = 1/252\)
- observation step: weekly (every \(m=5\) internal steps)
- sample length: \(T_{\text{week}} = 525\) weekly observations (≈ 10 years)
- measurement errors in option prices: **out of scope** here
- option pricing: **out of scope** here (will be a separate module/spec)

The goal is to implement clean, reusable code that can later be paired with a COS-based pricer to form option panels.

---

## 2. Model under \(\mathbb{P}\) (DGP)
Let \(y_t=\log(S_t)\). Under \(\mathbb{P}\) the DGP is:

\[
\frac{dS_t}{S_t} = (r-q+\eta V_t)\,dt + \sqrt{V_t}\,dW^{\mathbb{P}}_{t,S},
\]
\[
dV_t = \kappa(\bar V - V_t)\,dt + \sigma_v \sqrt{V_t}\,dW^{\mathbb{P}}_{t,V}, 
\quad dW^{\mathbb{P}}_{t,S} dW^{\mathbb{P}}_{t,V} = \rho\,dt.
\]

Equivalently,
\[
dy_t = (r-q)\,dt + \left(\eta-\tfrac12\right)V_t\,dt + \sqrt{V_t}\,dW^{\mathbb{P}}_{t,S}.
\]

Parameters (P):
- \(\eta,\kappa,\bar V,\sigma_v,\rho\)
- \(r,q\) (constants for now; can be generalized later)

---

## 3. Simulation design (daily internal / weekly output)

### 3.1 Time grids
- daily grid: \(t_n = n\delta\), \(\delta = 1/252\)
- weekly grid: record every \(m=5\) days, \(\Delta = m\delta\)

### 3.2 Lengths
- weekly observations recorded: `T_week = 525`
- daily steps recorded (excluding burn-in): `N_days = m * T_week = 2625`
- burn-in (daily steps): configurable; default **756 days** (≈ 3 years; conservative)

### 3.3 Initialization
- \(S_0 = 100\) by default (configurable)
- \(V_0 = \bar V\) by default (configurable)

---

## 4. Variance stepping: Andersen QE scheme (default choice)

We simulate \(V_{n+1}\) from \(V_n\) using Andersen’s Quadratic-Exponential (QE) scheme applied to the CIR process.

### 4.1 Conditional mean/variance of CIR (exact moments)
For a step size \(\delta\), define:
- \(e = \exp(-\kappa \delta)\)
- conditional mean:
  \[
  m = \bar V + (V_n-\bar V)e
  \]
- conditional variance:
  \[
  s^2 = V_n \frac{\sigma_v^2 e(1-e)}{\kappa} + \bar V \frac{\sigma_v^2 (1-e)^2}{2\kappa}
  \]
- \(\psi = s^2/m^2\)

### 4.2 QE branching rule
Use threshold \(\psi_c = 1.5\).

#### Case A: \(\psi \le \psi_c\)  (quadratic Gaussian form)
Draw \(Z \sim \mathcal{N}(0,1)\) and set:
- \(b^2 = \frac{2}{\psi}-1 + \sqrt{\frac{2}{\psi}}\sqrt{\frac{2}{\psi}-1}\)
- \(a = \frac{m}{1+b^2}\)
- \(V_{n+1} = a(b+Z)^2\)

#### Case B: \(\psi > \psi_c\) (point mass at 0 + exponential)
Draw \(U \sim \mathcal{U}(0,1)\), set:
- \(p = \frac{\psi-1}{\psi+1}\)
- \(\beta = \frac{1-p}{m}\)

Then:
- if \(U \le p\): \(V_{n+1}=0\)
- else:
  - \(U_2 = \frac{U-p}{1-p}\) (uniform on (0,1))
  - \(V_{n+1} = -\frac{1}{\beta}\ln(\max(U_2,\varepsilon))\), with \(\varepsilon\) ~ 1e-16

### 4.3 Positivity
The QE update must produce \(V_{n+1}\ge 0\) by construction.

---

## 5. Log-price update (correlation-preserving, uses variance increment)

After computing \(V_{n+1}\), define the integrated variance proxy (trapezoid):
\[
\widehat I_n = \frac{\delta}{2}(V_n+V_{n+1}).
\]

Draw \(Z^\perp \sim \mathcal{N}(0,1)\), independent of the QE draws.

Use the correlation-preserving decomposition:
\[
\Delta y_n = y_{n+1}-y_n
= (r-q)\delta
+ \left(\eta-\tfrac12\right)\widehat I_n
+ \frac{\rho}{\sigma_v}\left[(V_{n+1}-V_n) - \kappa(\bar V \delta - \widehat I_n)\right]
+ \sqrt{1-\rho^2}\sqrt{\widehat I_n}\,Z^\perp.
\]

Then:
- \(y_{n+1} = y_n + \Delta y_n\)
- \(S_{n+1} = \exp(y_{n+1})\)

Notes:
- This update preserves the leverage channel via the variance increment term, without needing an explicit \(dW_V\) draw.
- Requires \(\widehat I_n \ge 0\); with QE and trapezoid, \(\widehat I_n\) is nonnegative.

---

## 6. Output (weekly observations)

After burn-in:
- record weekly indices \(n = n_0, n_0+m, n_0+2m, \dots, n_0+mT_{\text{week}}\)
- output arrays of length `T_week+1`:
  - `t_week[i]` (in years, i.e., \(i\Delta\))
  - `logS_week[i]`
  - `V_week[i]`
- output weekly log-returns of length `T_week`:
  - `dlogS_week[i] = logS_week[i+1] - logS_week[i]`

Optional (debug only):
- ability to also return daily series (off by default to reduce memory).

---

## 7. Software design (separation of concerns / SOLID)

### 7.1 Modules (recommended)
- `Python/sim/`  
  - `HestonPathSimulator` (or `HestonSimulator`)
  - `HestonParamsP`, `HestonSimConfig`, `HestonPath`
  - `VarianceStepperAndersenQE`

- `Python/pricing/`  
  (not implemented in this task; future)
  - COS-based Heston pricer, option contract types, panel generator

### 7.2 Interfaces (dependency injection)
- `IVarianceStepper`  
  - `double step(double V_n, double delta, RNG& rng)`  
- `HestonPathSimulator` accepts:
  - `HestonParamsP`
  - `HestonSimConfig`
  - `IVarianceStepper&` (default: Andersen QE)

### 7.3 RNG requirements
- deterministic engine (e.g., `mt19937_64` in C++)
- controlled seeding via config
- avoid global RNG state

---

## 8. Acceptance tests (must be implemented)

### 8.1 Determinism
- same seed + same config → identical weekly outputs (bitwise or tight tolerance)

### 8.2 Positivity
- verify \(V_n \ge 0\) for all n across:
  - typical parameter sets
  - stress sets (high vol-of-vol, low kappa, high rho magnitude)

### 8.3 QE conditional moments (unit test on the stepper)
For fixed \((V_n, \delta)\), compute analytic \(m, s^2\) (Section 4.1).
Generate many draws of \(V_{n+1}\) and assert:
- sample mean close to \(m\)
- sample variance close to \(s^2\)

### 8.4 Correlation sanity (simulation-level)
Over many replications:
- `Corr(dy_daily, dV_daily)` is monotone in \(\rho\) and matches the sign of \(\rho\)

---

## 9. Micro-benchmark harness (to confirm cost profile)
Add a benchmark (executable or script) that reports:
- ns/step for `VarianceStepperAndersenQE::step`
- ns/step for full daily update (V + logS)
- total time to generate one replication:
  - burn-in days + 2625 daily steps
- allow varying `burnin_days` to confirm it is negligible relative to later pricing/PF work

---

## 10. Implementation checklist (Codex task order)
1. Implement `HestonParamsP`, `HestonSimConfig`, `HestonPath`
2. Implement `VarianceStepperAndersenQE` per Section 4
3. Implement `HestonPathSimulator` per Sections 5–6
4. Add tests 8.1–8.4
5. Add micro-benchmark (Section 9)
