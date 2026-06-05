# sadra_opf.jl
# Build function and top-level solve entry points for SADRA AC/DC OPF.

import PowerModels as _PM
using JuMP

function build_sadra_opf(pm::_PM.AbstractACPModel)

    # Populate SADRA ref lookups (must run after PM has built :branch, :gen, :bus)
    sadra_build_ref!(pm)

    # --- Variables ---
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    variable_sadra_converter(pm)

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
    sadra_ids      = union(vsc_branch_ids, dc_branch_ids)

    # --- Standard AC branches ---
    for i in _PM.ids(pm, :branch)
        i in sadra_ids && continue
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i)
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
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

        # Optional control setpoints (eq 22-23). Disabled by default in v0.1:
        # standard OPF leaves converter power free to minimise cost, matching
        # PowerModelsACDC's run_acdcopf. To enable fixed-setpoint control,
        # uncomment the block below.
        #
        # branch = _PM.ref(pm, :branch, i)
        # if branch["sadra_type_dc"] == 1
        #     constraint_vsc_p_setpoint(pm, i)
        # end
        # if branch["sadra_type_ac"] == 1
        #     constraint_vsc_q_setpoint(pm, i)
        # end
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
