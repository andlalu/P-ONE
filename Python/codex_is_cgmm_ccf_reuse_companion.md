# Companion Note: Reuse the Existing Analytic Heston CCF Solver for IS-CGMM Transition CCF

## Purpose

This note is intended to be used **alongside**:

```text
codex_is_cgmm_task0_followup.md
```

It refines the transition-characteristic-function part of that instruction set.

The key point is:

> The current option-pricing code already contains an analytic Heston CCF solver. Do not implement a separate analytic Heston CCF from scratch unless reuse/refactoring is impossible. First inspect and refactor the existing solver so the same analytic Riccati machinery can serve both option pricing and IS-CGMM transition moments.

Current relevant file:

```text
Python/OptionPricing/heston_ccf_solver.py
```

Current relevant IS-CGMM file:

```text
Python/Estimation/ISCGMM/heston_transition_cf.py
```

---

## Background

The current option-pricing solver is named:

```python
HestonAnalyticCcfSolver
```

and implements an analytic Heston affine-coefficient formula for option pricing. The IS-CGMM implementation currently has a separate RK4-based Riccati solver for the physical-measure transition CCF.

This duplication should be reviewed. The two objects are related, but not identical.

---

## Distinguish the two CCFs

### 1. Option-pricing CCF

The option-pricing CCF is used by the COS pricer under the risk-neutral measure.

It is essentially a transform of the future log price or log-forward variable, with terminal variance-transform coefficient equal to zero. In affine notation it computes coefficients:

```math
E^Q[exp(i u X_{t+tau}) | X_t, V_t]
=
exp(A_Q(tau,u)+B_Q(tau,u)V_t+i u X_t).
```

In the current code this is implemented by `HestonAnalyticCcfSolver.solve_coefficients(...)`.

### 2. IS-CGMM transition CCF

The IS-CGMM criterion needs the conditional CCF of the next estimation state under the physical measure:

```math
X_i^theta=(R_i,V_i^theta),
```

where `R_i` is a log-return increment and `V_i^theta` is the implied Heston variance state.

The required transform is:

```math
phi_P(u,w;v,Delta)
=
E^P_theta[
exp(i u R_{t,t+Delta}+i w V_{t+Delta})
| V_t=v
].
```

This is not exactly the same object as the option-pricing CCF because:

1. it is under `P`, not `Q`;
2. the return drift contains `eta`;
3. the terminal transform includes `i w V_{t+Delta}`, so the Riccati equation has initial condition `B(0)=i w`, not `B(0)=0`;
4. the estimator conditions on the current variance state `V_t`, not on an option-pricing forward log variable.

Therefore, do not simply call the option-pricing solver unchanged. Instead, refactor the common analytic Riccati logic.

---

## Required implementation strategy

### Step 1: Inspect the existing analytic solver

Inspect:

```text
Python/OptionPricing/heston_ccf_solver.py
```

Confirm that it implements the closed-form Heston Riccati coefficients for option pricing.

The solver currently contains the Heston square-root discriminant, `g`-ratio, logarithmic term, and affine `A/B` coefficients. This is the right place to start.

### Step 2: Extract a reusable analytic Riccati helper

Create a small internal helper that solves the scalar constant-coefficient Riccati equation:

```math
B'(tau)=a+bB(tau)+cB(tau)^2,
B(0)=B_0,
```

and also returns:

```math
int_0^Delta B(tau) dtau.
```

Suggested location:

```text
Python/OptionPricing/heston_riccati.py
```

or, if the project style prefers fewer files:

```text
Python/OptionPricing/heston_ccf_solver.py
```

as internal helper functions.

The helper should be Heston-specific or at most scalar-Riccati-specific. Do not build a generic affine-jump-diffusion framework.

Suggested function signature:

```python
def solve_constant_riccati(
    *,
    a: np.ndarray,
    b: np.ndarray,
    c: float,
    b0: np.ndarray,
    tau: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return B(tau) and integral_0_tau B(s) ds."""
```

It must support complex NumPy arrays.

### Step 3: Re-express option-pricing solver using this helper if safe

If the refactor is low risk, modify `HestonAnalyticCcfSolver.solve_coefficients(...)` to use the shared helper.

If that introduces too much risk, leave the option-pricing solver unchanged and use the helper only for the IS-CGMM transition CCF. But document why the option-pricing solver was not refactored.

### Step 4: Implement analytic IS-CGMM transition CCF using the shared helper

For the physical-measure transition CCF:

```math
phi_P(u,w;v,Delta)
=
exp(A(Delta;u,w)+B(Delta;u,w)v).
```

The Riccati coefficients should be:

```math
B' = a+bB+cB^2,  B(0)=i w,
```

with:

```math
a = i u (eta-1/2)-1/2 u^2,
b = -kappa+rho sigma_v i u,
c = 1/2 sigma_v^2.
```

The affine `A` coefficient is:

```math
A(Delta)
=
i u(r-q)Delta
+
kappa * vbar * int_0^Delta B(s) ds.
```

Check signs carefully against the current RK implementation in `heston_transition_cf.py`.

### Step 5: Keep RK4 as a benchmark/fallback

Do not remove the RK implementation.

Expose a method option:

```python
heston_p_transition_cf(..., method="analytic")
```

Accepted values:

```text
analytic
rk4
```

The default may be `"analytic"` only if the tests described below pass. Otherwise default to `"rk4"` and document analytic as experimental.

Also expose this through:

```python
CgmmConfig
```

and the validation CLI:

```bash
--transition-cf-method analytic
--transition-cf-method rk4
```

---

## Analytic Riccati formula guidance

For:

```math
B' = a+bB+cB^2,  B(0)=B_0,
```

with `c != 0`, define:

```math
d=sqrt(b^2-4ac),
r_1=(-b+d)/(2c),
r_2=(-b-d)/(2c).
```

Then one possible form is:

```math
y_0=(B_0-r_1)/(B_0-r_2),
y_Delta=y_0 exp(d Delta),
B(Delta)=(r_1-y_Delta r_2)/(1-y_Delta).
```

The integral can be written as:

```math
int_0^Delta B(s) ds
=
r_1 Delta
+
((r_1-r_2)/d) [log(1-y_0)-log(1-y_Delta)].
```

Complex-log branch handling is the main numerical issue. If this form creates branch discontinuities, investigate an equivalent Heston-CCF “little trap” / stable-log representation or fall back to RK4.

Do not force a fragile analytic implementation into production.

---

## Required references / documentation to consult

Codex should consult the existing implementation first. If external references are available, use these as guidance for branch and formulation issues:

1. Heston (1993), “A Closed-Form Solution for Options with Stochastic Volatility...”
2. Albrecher, Mayer, Schoutens and Tistaert, “The Little Heston Trap.”
3. Kahl, Jäckel and Lord work on Heston implementation / branch-cut issues.
4. Lord and Kahl material on Heston characteristic-function implementation.
5. Cui, del Bano Rollin and Germano, “Full and Fast Calibration of the Heston Stochastic Volatility Model.”
6. del Bano Rollin, Ferreiro-Castilla and Utzet, “A New Look at the Heston Characteristic Function.”

The implementation should not add paper references to code comments unless useful. A short comment explaining the branch/log choice is enough.

---

## Tests required

### 1. Existing option-pricing tests must still pass

Any refactor of `HestonAnalyticCcfSolver` must preserve COS pricing behavior.

If there are no existing option-pricing regression tests, add a minimal test that checks the solver output shape and finite values for representative parameters.

### 2. Analytic transition CCF vs RK4

Add tests comparing analytic and RK4 transition CCF values for:

- weekly `Delta=5/252`;
- several `V_t` values;
- several `(u,w)` nodes, including moderate positive and negative arguments;
- the current dissertation parameter values:
  ```text
  eta=5, kappa=7, vbar=0.0225, sigma_v=0.4, rho=-0.5
  ```

The tolerance can be relaxed if RK4 is not extremely accurate, but the discrepancy should be small for `rk_steps=128` or higher.

### 3. Moment-derivative test

Keep or add the deterministic first-moment derivative test from the main instruction file.

For the transition CCF:

```math
E[V_{t+Delta}|V_t=v] = vbar + (v-vbar) exp(-kappa Delta)
```

and:

```math
E[R_{t,t+Delta}|V_t=v]
=
(r-q)Delta
+
(eta-1/2)
[
vbar Delta + (v-vbar)(1-exp(-kappa Delta))/kappa
].
```

Use finite differences of the CCF around zero to recover these moments.

### 4. Transition method passed through C-GMM criterion

Add a test that:

- constructs a tiny panel;
- evaluates the first-step criterion with `transition_cf_method="rk4"`;
- evaluates it with `transition_cf_method="analytic"`;
- verifies both are finite and close enough.

### 5. Slow Monte Carlo test remains gated

If the Monte Carlo test remains in the suite, gate it with `RUN_SLOW_TESTS=1` or a pytest slow marker.

---

## Acceptance criteria

The task is complete when:

1. The code first tries to reuse/refactor the existing analytic Heston CCF machinery rather than creating an unrelated second analytic implementation.
2. The IS-CGMM transition CCF supports explicit method selection:
   - `analytic`;
   - `rk4`.
3. The analytic method either:
   - passes analytic-vs-RK4 and moment-derivative tests and can be made default; or
   - is documented as experimental while RK4 remains default.
4. The validation CLI exposes the method selection.
5. Fast tests do not run the heavy Monte Carlo check by default.
6. The main Task 0 follow-up file remains valid and this companion file is implemented consistently with it.

---

## Suggested final Codex response

Report:

1. Whether the existing option-pricing analytic CCF solver was reused or refactored.
2. What common helper was introduced, if any.
3. Whether analytic IS-CGMM transition CCF is now available.
4. Which method is default and why.
5. Test results:
   - fast tests;
   - slow tests, if run.
6. Remaining concerns about branch cuts, complex logarithms, or numerical stability.
7. Example validation command using `--transition-cf-method analytic`.
