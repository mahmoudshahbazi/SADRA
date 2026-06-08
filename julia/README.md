# SADRA.jl

A [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) implementation
of the **SADRA** universal AC/DC branch model for optimal power flow in hybrid
AC/DC networks, including full VSC converter control modes and controlled
transformers.

SADRA models VSC-based HVDC converters and DC grids using only standard AC
power-flow equations, by representing each converter as a standard branch (with
variable tap ratio and phase shift) plus a dummy generator that accounts for
converter losses. This lets a hybrid AC/DC OPF be solved with the existing AC
machinery in PowerModels, with no DC-specific solver code.

The method is described in:

> M. Shahbazi, "An Efficient Universal AC/DC Branch Model for Optimal Power
> Flow Studies in Hybrid AC/DC Systems," *IEEE Transactions on Power Systems*,
> vol. 40, no. 4, pp. 3211-3221, July 2025. DOI: 10.1109/TPWRS.2024.3514815

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/SADRA")   # or Pkg.add once registered
```

Dependencies: PowerModels, JuMP, and a nonlinear solver such as Ipopt. (XLSX is
used only by the validation/comparison helpers.)

## Usage

SADRA accepts two input formats and offers an optional, per-run control switch.

### PowerModelsACDC-format cases

```julia
using SADRA, Ipopt, JuMP

opt = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6)

# Uncontrolled AC/DC OPF (converter setpoints free, minimise cost)
result = solve_sadra_opf("case3120sp_acdc.m", opt;
                         setting = Dict("sadra_controls" => false))

println(result["termination_status"])   # LOCALLY_SOLVED
println(result["objective"])            # 2,143,038
```

A PowerModelsACDC `.m` case contains `busdc`, `convdc` and `branchdc` sections.
`solve_sadra_opf` transforms these into an augmented pure-AC network, solves,
and returns a PowerModels-style result dictionary.

### MATPOWER-FUBM-format cases (with control actions)

The SADRA paper's 1354-bus PEGASE case is distributed in FUBM format, where the
DC data is encoded inline in extended branch columns. Use `solve_sadra_fubm`:

```julia
using SADRA, Ipopt, JuMP

opt = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6)

# Controlled AC/DC OPF (controls applied by default)
result = solve_sadra_fubm("test/data/sadra_case1354pegase_2MTDC_ctrls.m", opt)
```

### Control modes

Converter and transformer control actions (paper eq. 22-24, Table I/IV) are
applied automatically when present in the case. They can be disabled with
`setting = Dict("sadra_controls" => false)` to recover the cost-minimising
uncontrolled OPF, and restricted to specific converters (by AC bus number) for
debugging with `setting = Dict("sadra_control_buses" => [7282])`.

Supported control actions:

- **VSC active-power setpoint** (`Pf`), with the converter loss correctly
  included (`Pf = setpoint + P_dummy`).
- **VSC DC-voltage control** (`Vdc`).
- **VSC droop control** (`Pf = Pref - kdp*(Vdc - Vref) + P_dummy`).
- **VSC phase pinning** (`theta = 0`, type I with no power setpoint).
- **VSC AC-voltage control** (`Vac`).
- **Controlled phase-shifting transformer (PST)** — from-side active power.
- **Controlled tap-changing transformer (CTT)** — to-side voltage.

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
definition, loss equation, and (optionally) the control actions above.

The `julia/` package reads everything from the case file: the FUBM ingest
derives the control type from the `CONV_A` column, the setpoints from
`VT_SET`/`VF_SET`/`PF`/`KDP`, the loss coefficients from `ALPHA1/2/3`, and the
DC-bus bounds from the bus matrix. No case-specific values are hardcoded.

## Validation

SADRA has been validated against two independent references:

- **Controlled case (1354-bus PEGASE, 2 DC grids):** reproduces the AIMMS
  reference (paper Table IV-VI) on all five VSC control modes and both
  controlled transformers, with objective 74,037.87.
- **Uncontrolled case (3120-bus):** reproduces PowerModelsACDC's AC/DC OPF to
  within 0.003% (with matched converter impedance), objective 2,143,038, and
  solves faster on every case tested.

See [VALIDATION.md](VALIDATION.md) for the full numbers, the reproduction
steps, and the limitations.

## Limitations

See [VALIDATION.md](VALIDATION.md#known-limitations) for the full list. In
brief: a convergence tolerance of `tol = 1e-6` is recommended (the unscaled
3120-bus objective fails to converge at Ipopt's default `tol = 1e-8`); the
distributed 1354-bus case file has been corrected to the AIMMS reference and
therefore differs from the original FUBM distribution; LCC converters, storage,
unit commitment and SCOPF are not modelled; and validation to date covers the
two cases above.

## License

MIT License — see [LICENSE](LICENSE).

## Acknowledgement

Built on PowerModels.jl and validated against PowerModelsACDC.jl and the SADRA
AIMMS reference implementation. Thanks to Abraham Alvarez-Bustos for the FUBM
work that SADRA builds on.

## Contact

Dr Mahmoud Shahbazi — <mahmoud.shahbazi@durham.ac.uk>
