# sadra_opf.jl
# Build function and top-level solve entry points for SADRA AC/DC OPF.

import PowerModels as _PM
using JuMP

function build_sadra_opf(pm::_PM.AbstractACPModel)

    # Read control flag from solve setting (default: controls ON).
    # Set setting["sadra_controls"]=false to recover the v0.1 uncontrolled OPF.
    apply_controls = get(pm.setting, "sadra_controls", true)
    # Optional: restrict controls to a subset of VSC AC buses, for debugging.
    # setting["sadra_control_buses"] = [7282] enables only that converter's
    # controls; nothing/absent => all converters.
    control_buses = get(pm.setting, "sadra_control_buses", nothing)

    # Populate SADRA ref lookups (must run after PM has built :branch, :gen, :bus)
    sadra_build_ref!(pm)

    # --- Variables ---
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    variable_sadra_converter(pm)
    variable_xfmr_control(pm)

    # --- Objective ---
    _PM.objective_min_fuel_cost(pm)

    # --- Reference bus and voltage model ---
    _PM.constraint_model_voltage(pm)
    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    # --- Power balance at all buses (AC and DC) ---
    for i in _PM.ids(pm, :bus)
        _PM.constraint_power_balance(pm, i)
    end

    # Identify SADRA branches
    vsc_branch_ids = Set(keys(_PM.ref(pm, :sadra_vsc_branch)))
    dc_branch_ids  = Set(keys(_PM.ref(pm, :sadra_dc_branch)))
    xfmr_ctrl_ids  = Set(keys(_PM.ref(pm, :sadra_ctrl_xfmr)))
    sadra_ids      = union(vsc_branch_ids, dc_branch_ids, xfmr_ctrl_ids)

    # --- Standard AC branches ---
    for i in _PM.ids(pm, :branch)
        i in sadra_ids && continue
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i)
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    # --- Controlled transformers (PST/CTT): variable-tap Ohm's law + control ---
    for i in xfmr_ctrl_ids
        constraint_xfmr_ohms(pm, i)
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
        if apply_controls
            kind = _PM.ref(pm, :branch, i)["sadra_xfmr_ctrl"]
            if kind == "pst"
                constraint_pst_pf(pm, i)
            elseif kind == "ctt"
                constraint_ctt_vt(pm, i)
            end
        end
    end

    # --- VSC branches: SADRA model ---
    for i in vsc_branch_ids
        constraint_vsc_ohms_from(pm, i)
        constraint_vsc_ohms_to(pm, i)
        constraint_vsc_current(pm, i)
        constraint_vsc_losses(pm, i)
        #constraint_vsc_zero_qf(pm, i)
        constraint_vsc_qg_link(pm, i)      # <-- replacing above and more stable,; essentially guarantees no reactive power on the dc lines.
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)

        # --- Optional control actions (paper eq 22-24, Table I/IV) ---
        # Each VSC gets one DC-side and one AC-side control constraint.
        #   type_dc: 1=P_f setpoint (eq22), 2=Vdc pin, 3=droop (eq24)
        #   type_ac: 1=Q setpoint (eq23),   2=Vac pin (done via bus bounds in data.jl)
        # Disable entirely with setting["sadra_controls"]=false.
        if apply_controls
            branch  = _PM.ref(pm, :branch, i)
            ac_bus  = branch["sadra_ac_bus_i"]
            # subset filter for debugging
            if control_buses !== nothing && !(ac_bus in control_buses)
                continue
            end
            type_dc = branch["sadra_type_dc"]
            type_ac = branch["sadra_type_ac"]

            if type_dc == 1
                constraint_vsc_p_setpoint(pm, i)
            elseif type_dc == 2
                constraint_vsc_vdc_setpoint(pm, i)
            elseif type_dc == 3
                constraint_vsc_droop(pm, i)
            elseif type_dc == 0
                # type I with theta=0 (paper Table IV, VSC1): pin phase shift.
                JuMP.fix(_PM.var(pm, :phi, i), 0.0; force=true)
            end

            # type_ac==1: Q setpoint. type_ac==2: Vac already pinned via bus
            # bounds in sadra_transform!, so no extra constraint here.
            if type_ac == 1
                constraint_vsc_q_setpoint(pm, i)
            end
        end
    end

    # --- DC branches: purely resistive (x=0); standard Ohm's law gives q=0 ---
    for i in dc_branch_ids
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end
end

# ====================================================================
# Top-level solve entry points
# ====================================================================

"""
    solve_sadra_opf(file::String, optimizer; setting=Dict())

Parse a PowerModelsACDC-format `.m` case file, apply the SADRA transform,
and solve the AC/DC OPF using the SADRA universal branch model.
"""
function solve_sadra_opf(file::String, optimizer; setting=Dict{String,Any}())
    data = _PM.parse_file(file)
    return solve_sadra_opf(data, optimizer; setting=setting)
end

"""
    solve_sadra_fubm(file::String, optimizer; setting=Dict())

Parse a MATPOWER-FUBM `.m` case (DC data inline in extended branch columns,
control assignment from paper Table IV), apply the SADRA transform, and solve.
Use this for the 1354-bus PEGASE case.
"""
function solve_sadra_fubm(file::String, optimizer; setting=Dict{String,Any}())
    data = parse_fubm_acdc(file)
    return solve_sadra_opf(data, optimizer; setting=setting)
end

"""
    solve_sadra_opf(data::Dict, optimizer; setting=Dict())

As above, but takes an already-parsed PowerModels data dictionary.
"""
function solve_sadra_opf(data::Dict{String,Any}, optimizer; setting=Dict{String,Any}())
    sadra_transform!(data)
    return _PM.solve_model(
        data,
        _PM.ACPPowerModel,
        optimizer,
        build_sadra_opf;
        setting = setting
    )
end
