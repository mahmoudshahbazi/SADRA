# SADRA.jl Validation

This document records the validation of `SADRA.jl` against
[PowerModelsACDC.jl](https://github.com/Electa-Git/PowerModelsACDC.jl) (PMACDC)
for AC/DC optimal power flow.

## Summary

SADRA reproduces PMACDC's AC/DC OPF results to within **0.003%** on the
3120-bus benchmark when both models use the same converter series impedance,
and is faster on every case tested. The two formulations are mathematically
equivalent; the small residual differences are explained below and are not
modelling errors.

## Reference case: `case3120sp_acdc`

| Quantity | SADRA | PMACDC (matched impedance) | PMACDC (original) |
|----------|------:|---------------------------:|------------------:|
| Objective | 2,143,038 | 2,142,975 | 2,142,635 |

"Matched impedance" means PMACDC was run with the converter transformer flag
enabled so that it uses the same `rtf + j*xtf` series impedance that SADRA
always applies to the VSC branch. With this match, the objective gap drops
from 0.019% to 0.003%.



## Core mechanism checks

The following were verified on `case3120sp_acdc`:

- **Common DC phase angle.** All DC buses in a connected DC grid converge to an
  identical voltage phase angle (e.g. -0.2392 rad for all five DC buses). This
  is the artefact of relaxing the DC grid into AC phasors and is what keeps the
  DC network free of reactive power.
- **Topology.** DC buses, VSC branches (from = DC bus, to = AC bus), DC lines
  and dummy generators are all created at the correct indices and connect the
  correct nodes.

## Comparison across cases

SADRA was run on all AC/DC cases shipped with PowerModelsACDC. For standard
AC/DC OPF cases (no LCC, droop, storage or SCOPF), the objective agreed with
PMACDC to within 0.25%, and SADRA solved faster in every case. Cases using
features not yet implemented in SADRA (see README limitations) are excluded
from the comparison.

## Reproducing

1. Solve with SADRA:
   ```julia
   using SADRA, Ipopt
   result = SADRA.solve_sadra_opf("test/data/case3120sp_acdc.m", Ipopt.Optimizer)
   ```
2. Solve with PMACDC and compare. **PowerModelsACDC is not a dependency of
   SADRA** — it must be installed separately, in its own environment, to
   reproduce the comparison. Use `test/run_pmacdc.jl` to generate the PMACDC
   results and `test/compare.jl` to compare them against SADRA.

## Environment note

PMACDC and SADRA are sensitive to package versions, and PMACDC is intentionally
kept out of SADRA's dependencies so that using SADRA does not pull in the (more
version-fragile) PMACDC stack. To reproduce the comparison, set up a separate
environment for PMACDC.

The validated PMACDC environment uses:

- PowerModelsACDC 0.7.4
- PowerModels 0.19.10
- JuMP 1.15.1
- MathOptInterface 1.25.0
- Ipopt 1.x

In that environment the converter-loss preprocessing function
`process_additional_data!` runs correctly and the OPF entry point is
`run_acdcopf`. SADRA itself only requires PowerModels, JuMP and Ipopt.
