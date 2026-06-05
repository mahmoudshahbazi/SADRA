# SADRA.jl

A [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) implementation
of the **SADRA** universal AC/DC branch model for optimal power flow in hybrid
AC/DC networks.

SADRA models VSC-based HVDC converters and DC grids using only standard AC
power-flow equations, by representing each converter as a standard branch (with
variable tap ratio and phase shift) plus a dummy generator that accounts for
converter losses. This lets a hybrid AC/DC OPF be solved with the existing AC
machinery in PowerModels, with no DC-specific solver code.

The method is described in:

> M. Shahbazi et al., "SADRA: A Universal Branch Model for Steady-State
> Analysis of Hybrid AC/DC Networks," *IEEE Transactions on Power Systems*,
> vol. 40, no. 4, 2025. DOI: 10.1109/TPWRS.2024.3514815

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/SADRA")   # or Pkg.add once registered
```

Dependencies: PowerModels, JuMP, and a nonlinear solver such as Ipopt.

## Usage

```julia
using SADRA, Ipopt

# Solve AC/DC OPF from a PowerModelsACDC-format .m case file
result = solve_sadra_opf("case3120sp_acdc.m", Ipopt.Optimizer)

println(result["termination_status"])
println(result["objective"])
```

The input is a standard PowerModelsACDC `.m` case file containing `busdc`,
`convdc` and `branchdc` sections. SADRA transforms these into an augmented
pure-AC network internally, solves, and returns a PowerModels-style result
dictionary.

## How it works

`sadra_transform!` augments the parsed data dictionary:

- each **DC bus** becomes an AC bus with the DC voltage limits;
- each **DC line** becomes a purely resistive AC branch (`x = 0`);
- each **VSC converter** becomes an AC branch from the DC bus to the AC bus,
  with a variable tap ratio `ma` and phase shift `phi` and the converter
  transformer impedance `rtf + j*xtf`;
- a **dummy generator** is added at each DC bus to represent converter losses
  (`P_loss = a + b*I + c*I^2`) and to balance reactive power.

The OPF then uses standard PowerModels variables and power balance, with a
small set of SADRA-specific constraints for the converter Ohm's law, current
definition and loss equation.

## Validation

SADRA reproduces PowerModelsACDC's AC/DC OPF to within 0.003% on the 3120-bus
benchmark (with matched converter impedance) and is faster on every case
tested. See [VALIDATION.md](VALIDATION.md) for the full comparison, numbers and
reproduction steps. PowerModelsACDC is required only to reproduce the
comparison and is not a dependency of SADRA.

## Current limitations (v0.1)

This is an initial release covering standard VSC-based AC/DC OPF. Not yet
supported:

- **Converter control modes** — droop control and fixed P/V/Q setpoints
  (paper eq. 22-24). The OPF currently leaves converter set-points free to
  minimise cost, matching PowerModelsACDC's `run_acdcopf`. Control-mode
  constraints are present in the code but disabled by default.
- **LCC (line-commutated converters)** — only VSC is modelled.
- **Back-to-back VSC links** — `case5_b2bdc` is not yet handled.
- **Storage and unit commitment** — cases with storage components are not
  supported.
- **Security-constrained OPF (SCOPF)**.

Two PMACDC test cases (`case5_2grids_uc_hvdc`, `case67acdc_scopf`) currently
return locally-infeasible and are under investigation.

## Status

v0.1 — research preview. The core VSC/DC-grid OPF is validated; control modes
and the features listed above are planned for subsequent releases.

## License

MIT License — see [LICENSE](LICENSE).

## Acknowledgement

Built on PowerModels.jl and validated against PowerModelsACDC.jl.
