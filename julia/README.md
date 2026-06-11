# SADRA.jl

A [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) implementation
of the **SADRA** universal AC/DC branch model for optimal power flow and
(from v1.1) power flow in hybrid AC/DC networks, including full VSC converter
control modes and controlled transformers.

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

### Power flow (v1.1)

```julia
# PowerModelsACDC-format case
result = solve_sadra_pf("test/data/case5_acdc.m", opt)

# MATPOWER-FUBM-format case
result = solve_sadra_fubm_pf("test/data/sadra_case1354pegase_2MTDC_ctrls.m", opt)
```

The power flow is posed, as in PowerModels and PowerModelsACDC, as a
feasibility problem (no objective) solved by the same nonlinear solver:
generator active power is fixed at all generators except the slack, voltage
magnitude is fixed at PV and slack buses, and the reference angle is set.
Converter controls are **mandatory** in the power flow -- they are the PF
specification itself, with each converter control pinning one converter
degree of freedom. Each DC island must contain exactly one DC-voltage-
controlled (or droop) converter acting as the DC slack; this is enforced by
a determinacy check that errors on under-determined cases rather than
returning an arbitrary feasible point.

Unlike the OPF, the power flow **reports** limit violations rather than
enforcing them: thermal limits, angle-difference limits, voltage bounds and
converter current ratings are not imposed. Converter variables use wide
finite bounds (`ma` in [0.5, 2], `phi` in [-pi, pi], converter current
bounded below by 0 only) for solver conditioning; a power flow that is
infeasible because `ma` would need to leave [0.5, 2] indicates a genuinely
extreme operating point.

Note on conventions: on the PowerModelsACDC-format path, the converter
active-power setpoint `P_g` is pinned at the grid side (PCC), matching the
MatACDC/PowerModelsACDC convention. On the FUBM-format path the setpoint
keeps the DC-side convention (`Pf = setpoint + P_dummy`) used by the AIMMS
reference.

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

- **Power flow (v1.1):** validated internally by reproducing SADRA OPF
  states to machine precision (max bus-voltage deviation 9e-16 on
  case5_acdc, 1e-9 on a 10-converter meshed-DC case; 2e-5 on the
  3-zone/2-grid case24, fully attributed to active OPF inequality bounds
  that the unbounded PF legitimately drifts off), and externally against
  PowerModelsACDC's AC/DC power flow (max bus-voltage difference 1.2e-4 on
  case5_acdc and 4.6e-4 on the 10-converter case, attributed to the
  converter-station representation and the loss-measurement convention).
  Known comparison caveats: PowerModelsACDC's PF can leave a Vac-controlled
  converter's AC degree of freedom unpinned (observed on
  case24_3zones_acdc bus 204), and at a voltage-controlled bus that also
  hosts free-Q generators the reactive split between converter and
  generators is non-unique (the bus state is unique). The script
  `test_pf_case.jl` reproduces all of these checks.

See [VALIDATION.md](VALIDATION.md) for the full numbers, the reproduction
steps, and the limitations.

## Changes vs v1.0 (important)

The bipolar DC branch resistance is now divided by `dcpol` on the
PowerModelsACDC-format path, as required by the MatACDC convention
(`P = dcpol * (1/r) * v_f * (v_f - v_t)`). In v1.0 the effective DC line
resistance on this path was a factor `dcpol` (typically 2) too large, so
**v1.0 results on PowerModelsACDC-format cases are not reproduced by
v1.1**. The FUBM-format path carries no `dcpol` key and is unaffected; the
1354-bus PEGASE reference objective is unchanged. The PowerModelsACDC-path
active-power setpoint convention also moved from the DC side to the grid
side (see Power flow notes above), affecting controlled OPF and PF on that
path.

## Limitations

See [VALIDATION.md](VALIDATION.md#known-limitations) for the full list. In
brief: a convergence tolerance of `tol = 1e-6` is recommended (the unscaled
3120-bus objective fails to converge at Ipopt's default `tol = 1e-8`); the
distributed 1354-bus case file has been corrected to the AIMMS reference and
therefore differs from the original FUBM distribution; LCC converters, storage,
unit commitment and SCOPF are not modelled; validation to date covers the
cases above; droop control on the PowerModelsACDC-format path is not yet
exercised by a power-flow test case; the FUBM-format power flow requires a
convention for converters with no AC-side control (`type_ac = 0`), still to
be settled; and the converter station is represented by a single branch
(transformer impedance), so the filter shunt and phase-reactor impedance of
the PowerModelsACDC station model are not yet represented (the ~1e-4
power-flow agreement floor; planned for v1.2).

## License

MIT License — see [LICENSE](LICENSE).

## Acknowledgement

Built on PowerModels.jl and validated against PowerModelsACDC.jl and the SADRA
AIMMS reference implementation. Thanks to Abraham Alvarez-Bustos for the FUBM
work that SADRA builds on.

## Contact

Dr Mahmoud Shahbazi — <mahmoud.shahbazi@durham.ac.uk>
