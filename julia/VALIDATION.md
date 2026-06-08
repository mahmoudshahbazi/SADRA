# SADRA.jl Validation

This document records the validation of `SADRA.jl` against two independent
references:

1. the **AIMMS** reference implementation of SADRA (the paper's reference), for
   the controlled 1354-bus PEGASE case; and
2. [PowerModelsACDC.jl](https://github.com/Electa-Git/PowerModelsACDC.jl)
   (PMACDC), for the uncontrolled 3120-bus case.

## Summary

| Case | Controls | Reference | SADRA objective | Reference | Status |
|------|----------|-----------|----------------:|----------:|--------|
| 1354-bus PEGASE (2 DC grids) | full (5 VSC + PST + CTT) | AIMMS | 74,037.87 | 74,037 | LOCALLY_SOLVED |
| 3120-bus | none (cost-min) | PMACDC | 2,143,038 | 2,142,975 | LOCALLY_SOLVED |

Both cases use the solver setting `tol = 1e-6` (see
[Known limitations](#known-limitations)).

---

## Controlled case: 1354-bus PEGASE

The modified PEGASE system from the SADRA paper (1354 AC buses, two DC
subgrids: one point-to-point HVDC link and one 3-node meshed MTDC grid; five
VSCs; one controlled phase-shifting transformer (PST) and one controlled
tap-changing transformer (CTT)). Run with:

```julia
using SADRA, Ipopt, JuMP
opt = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6)
res = solve_sadra_fubm("test/data/sadra_case1354pegase_2MTDC_ctrls.m", opt)
```

### Control actions (paper Table IV) — all reproduced

| Element | Control | Target | SADRA result |
|---------|---------|--------|--------------|
| VSC1 | theta = 0, Vac | Vac = 1.08 | matched |
| VSC2 | droop, Vac | Vac = 1.075 | matched |
| VSC3 | Vdc, Vac | Vdc = 1.01, Vac = 1.06 | matched |
| VSC4 | Pf, Vac | Pf = -5.0 pu, Vac = 1.07 | Pf = -5.0022, matched |
| VSC5 | Pf, free | Pf = -4.5 pu | Pf = -4.5020, matched |
| PST  | from-side P | Pf = -1.73 pu | -1.7300 |
| CTT  | to-side V | Vt = 1.1 | 1.1000 |

The VSC4/VSC5 active-power values land at the setpoint **plus the converter
loss** (e.g. -5.0022 = -5.0 + dummy-generator power), exactly matching the
AIMMS `Pf = PfShiftControl + Pg` relation. DC line reactive power is ~1e-17 pu;
each DC grid converges to a single common voltage phase angle. Objective
74,037.87 matches the AIMMS reference (74,037) to displayed precision.

### How the controlled case is built

All control data is read **from the case file**, not hardcoded:

- control **type** from the `CONV_A` column (1 = type I, 2 = type II,
  4 = type III), with the type-I theta=0 / Pf split decided by whether a `PF`
  setpoint is present;
- **Vac** from `VT_SET`, **droop Vref** from `VF_SET`, **droop slope** from
  `KDP`, **Pf / Pref** from `PF`;
- **loss coefficients** from `ALPHA1/2/3`;
- **DC-bus voltage bounds** from the bus matrix.

Controlled transformers are detected structurally (a non-VSC transformer branch
carrying a `PF` or `VT_SET` setpoint), with no hardcoded element list.

> **Note on the case file.** The distributed
> `sadra_case1354pegase_2MTDC_ctrls.m` has been **corrected to the AIMMS
> reference values** and therefore differs from the original FUBM distribution.
> The original file's control columns (`VT_SET`, `ALPHA1/2/3`, droop columns)
> did not match the values that produced the paper's results; they have been
> updated, and a header comment in the `.m` documents the changes.

---

## Uncontrolled case: 3120-bus

```julia
using SADRA, Ipopt, JuMP
opt = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6)
res = solve_sadra_opf("test/data/case3120sp_acdc.m", opt;
                      setting = Dict("sadra_controls" => false))
```

| Quantity | SADRA | PMACDC (matched impedance) | PMACDC (original) |
|----------|------:|---------------------------:|------------------:|
| Objective | 2,143,038 | 2,142,975 | 2,142,635 |
| Status | LOCALLY_SOLVED | LOCALLY_SOLVED | NUMERICAL_ERROR |
| Objective gap vs SADRA | — | 0.003% | 0.019% |

"Matched impedance" means PMACDC was run with the converter transformer flag
enabled, so it uses the same `rtf + j*xtf` series impedance SADRA applies to
the VSC branch. With this match the objective gap drops to 0.003% and PMACDC
converges cleanly (it returns NUMERICAL_ERROR with the flag off, while SADRA
solves the same system without difficulty).

### Why a small gap remains (0.003%)

Even with matched series impedance, PMACDC retains an intermediate filter bus
between the transformer and the converter — one extra node per converter.
SADRA lumps the converter station into a single branch. The two networks are
therefore not topologically identical, giving marginally different voltage
profiles and losses. A difference at the 0.003% level is well within the
tolerance expected for two locally-optimal solutions of a non-convex problem.

### DC grid flow non-uniqueness

On meshed DC grids the per-branch DC power distribution is not unique at the
optimum. For `case3120sp_acdc` the DC buses 1-2-3 form a triangle, so multiple
DC line flow patterns achieve the same objective and (near-)identical bus
voltages. SADRA and PMACDC may report different DC line flows (up to ~12 MW on
individual lines) while agreeing on objective and voltages. Both are valid
optima.

---

## Core mechanism checks

Verified on both cases:

- **Common DC phase angle.** All DC buses in a connected DC grid converge to an
  identical voltage phase angle — the artefact of relaxing the DC grid into AC
  phasors, and what keeps the DC network free of reactive power.
- **Reactive-free DC lines.** DC lines are purely resistive (`x = 0`);
  combined with the common phase angle, DC-line reactive flow is ~0 (≤1e-17 pu).
- **Physical converter losses.** Dummy-generator active power equals the
  converter loss model `P_loss = a + b*I + c*I^2`, matching the AIMMS
  dummy-generator powers exactly on the 1354-bus case.

---

## Known limitations

These are stated plainly so the 1.0 release is not over-claimed.

- **Convergence tolerance.** A setting of `tol = 1e-6` is recommended. At
  Ipopt's default `tol = 1e-8` the uncontrolled 3120-bus case fails with
  `Restoration Failed`, because the unscaled objective (~2e6) is poorly
  conditioned for the interior-point method. `tol = 1e-6` is a case-independent
  engineering tolerance that resolves this; an objective scaling factor was
  considered but rejected as case-specific.
- **Corrected case file.** The 1354-bus `.m` has been corrected to the AIMMS
  reference and differs from the original FUBM distribution (see note above).
- **Validation scope.** Validation to date covers the two cases above (one
  controlled vs AIMMS, one uncontrolled vs PMACDC). Behaviour on arbitrary
  cases for which no independent reference was cross-checked is not guaranteed.
- **Droop voltage reference** is read from the `VF_SET` column; cases without
  it set will default to 1.0 p.u.
- **Not modelled:** LCC (line-commutated) converters; storage and unit
  commitment; security-constrained OPF (SCOPF).
- **Residual.** A persistent bound residual of ~1.6e-7 appears at the solution
  of the controlled case; it is within the `tol = 1e-6` criterion and does not
  affect the reported values, but its exact origin was not traced.

---

## Reproducing the PMACDC comparison

PowerModelsACDC is **not a dependency of SADRA** — install it separately, in
its own environment, to reproduce the 3120-bus comparison. The validated PMACDC
environment uses PowerModelsACDC 0.7.4, PowerModels 0.19.10, JuMP 1.15.1,
MathOptInterface 1.25.0, Ipopt 1.x. In that environment the converter-loss
preprocessing `process_additional_data!` runs correctly and the OPF entry point
is `run_acdcopf`. SADRA itself requires only PowerModels, JuMP and Ipopt.
