# sadra_pf.jl
# Power flow (PF) problem for the SADRA universal branch model.
#
# Same physics as build_sadra_opf, but posed as a feasibility problem
# (no objective) with a square specification:
#   - AC slack (ref bus): vm and va fixed; its gens absorb the residual.
#   - PV buses: vm fixed to bus["vm"], gen pg fixed to gen["pg"], qg free.
#   - PQ buses: load fixed (data), vm/va free.
#   - VSC controls are MANDATORY (they ARE the PF specification):
#       type_dc 1/2/3/0 and type_ac 1/2 each pin one converter DOF.
#   - Dummy loss generators are NOT fixed: their pg comes from the loss
#     equation and qg from the qg-link. Fixing them would over-determine.
#
# Setpoint precedence (validated against the 1354 PEGASE case, bus 7282):
# if a bus is both PV and Vac-pinned by a converter (type_ac==2) or a CTT,
# the converter/transformer setpoint wins -> the gen-V fix is SKIPPED there.
# Gen P is still fixed; gen Q stays free. One vm pin per bus, never two.
#
# NOTE vs textbook PF: variable bounds are kept (the Vac pin lives in the
# bus bounds set by sadra_transform!), but thermal limits and angle-
# difference constraints are dropped — a PF state may violate ratings,
# and PF must report it, not exclude it.

import PowerModels as _PM
using JuMP

# ------------------------------------------------------------------
# Determinacy guard: every DC island must carry exactly one voltage
# anchor (one Vdc pin, type_dc==2) OR at least one droop (type_dc==3).
# Zero anchors  -> island voltage floats (under-determined; Ipopt will
#                  return an arbitrary feasible level and report success).
# Two+ Vdc pins -> over-determined (errors unless values happen to agree).
# ------------------------------------------------------------------
function sadra_check_pf_determinacy(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)
    vsc = _PM.ref(pm, nw, :sadra_vsc_branch)
    dcbr = _PM.ref(pm, nw, :sadra_dc_branch)

    # union-find over DC buses
    parent = Dict{Int,Int}()
    find(x) = (parent[x] == x ? x : (parent[x] = find(parent[x])))
    function union!(a, b)
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end
    for (_, br) in vsc
        parent[br["f_bus"]] = get(parent, br["f_bus"], br["f_bus"])
    end
    for (_, br) in dcbr
        parent[br["f_bus"]] = get(parent, br["f_bus"], br["f_bus"])
        parent[br["t_bus"]] = get(parent, br["t_bus"], br["t_bus"])
        union!(br["f_bus"], br["t_bus"])
    end

    anchors = Dict{Int,Vector{String}}()   # island root -> anchor list
    for (i, br) in vsc
        root = find(br["f_bus"])
        td = br["sadra_type_dc"]
        if td == 2
            push!(get!(anchors, root, String[]), "VSC branch $i: Vdc pin")
        elseif td == 3
            push!(get!(anchors, root, String[]), "VSC branch $i: droop")
        else
            get!(anchors, root, String[])
        end
    end

    for (root, a) in anchors
        nvdc = count(s -> occursin("Vdc pin", s), a)
        if isempty(a)
            error("SADRA PF: DC island (root bus $root) has no voltage anchor " *
                  "(no converter with type_dc==2 or 3). The PF is under-determined.")
        elseif nvdc > 1
            @warn "SADRA PF: DC island (root bus $root) has $nvdc Vdc pins: $(a). " *
                  "Over-determined unless setpoints are consistent."
        end
    end
    return nothing
end

# ------------------------------------------------------------------
# PF problem builder
# ------------------------------------------------------------------
function build_sadra_pf(pm::_PM.AbstractACPModel)

    # Controls are not optional in PF — they make the system square.
    if get(pm.setting, "sadra_controls", true) == false
        error("SADRA PF: setting[\"sadra_controls\"]=false is not valid for PF; " *
              "converter controls are the PF specification.")
    end

    sadra_build_ref!(pm)
    sadra_check_pf_determinacy(pm)

    # --- Variables: PM variables UNBOUNDED (as in PM's build_pf).
    # Bounded vars make PF infeasible whenever a vm fix needs qg beyond
    # gen limits (PV->PQ switching territory) or vm outside bus bounds.
    # The Vac pin no longer lives in bus bounds: it is an explicit
    # equality below. 
    _PM.variable_bus_voltage(pm, bounded=false)
    _PM.variable_gen_power(pm, bounded=false)
    _PM.variable_branch_power(pm, bounded=false)
    # SADRA converter variables: WIDE but finite bounds in PF.
    # PF must not enforce converter ratings (it reports violations), but
    # fully unbounded ma/phi degrade Ipopt conditioning (case5/case24 fail
    # to converge: phi is 2pi-periodic, ma scale-free). Wide finite bounds
    # keep the search region sane; iconv keeps only its structural >= 0.
    variable_sadra_converter(pm)
    for i in _PM.ids(pm, :sadra_vsc_branch)
        JuMP.set_lower_bound(_PM.var(pm, :ma)[i],  0.2)
        JuMP.set_upper_bound(_PM.var(pm, :ma)[i],  2.0)
        JuMP.set_lower_bound(_PM.var(pm, :phi)[i], -pi)
        JuMP.set_upper_bound(_PM.var(pm, :phi)[i],  pi)
        ic = _PM.var(pm, :iconv_ac)[i]
        JuMP.has_upper_bound(ic) && JuMP.delete_upper_bound(ic)
    end
    variable_xfmr_control(pm)

    # --- No objective: feasibility problem (square system) ---

    # --- Voltage model + AC slack ---
    _PM.constraint_model_voltage(pm)
    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
        _PM.constraint_voltage_magnitude_setpoint(pm, i)   # vm = bus["vm"]
    end

    # --- Power balance at all buses ---
    for i in _PM.ids(pm, :bus)
        _PM.constraint_power_balance(pm, i)
    end

    # Identify SADRA branches
    vsc_branch_ids = Set(keys(_PM.ref(pm, :sadra_vsc_branch)))
    dc_branch_ids  = Set(keys(_PM.ref(pm, :sadra_dc_branch)))
    xfmr_ctrl_ids  = Set(keys(_PM.ref(pm, :sadra_ctrl_xfmr)))
    sadra_ids      = union(vsc_branch_ids, dc_branch_ids, xfmr_ctrl_ids)

    # Buses whose vm is pinned by a converter Vac pin or a CTT:
    # the gen-V fix must be skipped there (converter setpoint wins).
    vm_pinned = Set{Int}()
    for i in vsc_branch_ids
        br = _PM.ref(pm, :branch, i)
        br["sadra_type_ac"] == 2 && push!(vm_pinned, br["sadra_ac_bus_i"])
        # DC-bus vm pinned by Vdc control: also exclude from gen-V logic
        br["sadra_type_dc"] == 2 && push!(vm_pinned, br["f_bus"])
    end
    for i in xfmr_ctrl_ids
        br = _PM.ref(pm, :branch, i)
        br["sadra_xfmr_ctrl"] == "ctt" && push!(vm_pinned, br["t_bus"])
    end

    # DC buses (host dummy gens only — no PV logic there)
    dc_buses = Set(br["f_bus"] for (_, br) in _PM.ref(pm, :sadra_vsc_branch))
    for (_, br) in _PM.ref(pm, :sadra_dc_branch)
        push!(dc_buses, br["f_bus"]); push!(dc_buses, br["t_bus"])
    end

    # --- PF specification: gen P fixes and PV-bus V fixes ---
    ref_buses = Set(_PM.ids(pm, :ref_buses))
    for (j, gen) in _PM.ref(pm, :gen)
        get(gen, "sadra_dummy_gen", false) && continue       # loss gens stay free
        b = gen["gen_bus"]
        b in ref_buses && continue                            # slack gens stay free
        _PM.constraint_gen_setpoint_active(pm, j)             # pg = gen["pg"]
        # A gen at a PQ bus gets no vm fix, so fix its qg too (PQ gen).
        if _PM.ref(pm, :bus, b)["bus_type"] == 1
            _PM.constraint_gen_setpoint_reactive(pm, j)       # qg = gen["qg"]
        end
    end
    for (b, bus) in _PM.ref(pm, :bus)
        b in ref_buses && continue
        b in dc_buses && continue
        b in vm_pinned && continue                            # converter/CTT wins
        if bus["bus_type"] == 2 && length(_PM.ref(pm, :bus_gens, b)) > 0
            _PM.constraint_voltage_magnitude_setpoint(pm, b)  # vm = bus["vm"]
        end
    end

    # --- Standard AC branches: Ohm's law only (no limits in PF) ---
    for i in _PM.ids(pm, :branch)
        i in sadra_ids && continue
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
    end

    # --- Controlled transformers: Ohm's law + mandatory control ---
    for i in xfmr_ctrl_ids
        constraint_xfmr_ohms(pm, i)
        kind = _PM.ref(pm, :branch, i)["sadra_xfmr_ctrl"]
        if kind == "pst"
            constraint_pst_pf(pm, i)
        elseif kind == "ctt"
            constraint_ctt_vt(pm, i)
        end
    end

    # --- VSC branches: physics + mandatory controls ---
    for i in vsc_branch_ids
        constraint_vsc_ohms_from(pm, i)
        constraint_vsc_ohms_to(pm, i)
        constraint_vsc_current(pm, i)
        constraint_vsc_losses(pm, i)
        constraint_vsc_qg_link(pm, i)

        branch  = _PM.ref(pm, :branch, i)
        type_dc = branch["sadra_type_dc"]
        type_ac = branch["sadra_type_ac"]

        if type_dc == 1
            constraint_vsc_p_setpoint(pm, i)
        elseif type_dc == 2
            constraint_vsc_vdc_setpoint(pm, i)
        elseif type_dc == 3
            constraint_vsc_droop(pm, i)
        elseif type_dc == 0
            JuMP.fix(_PM.var(pm, :phi, i), 0.0; force=true)
        end

        if type_ac == 1
            constraint_vsc_q_setpoint(pm, i)
        elseif type_ac == 2
            # Bounds are off in PF, so pin Vac explicitly to sadra_Vtar.
            vm_ac = _PM.var(pm, :vm, branch["sadra_ac_bus_i"])
            JuMP.@constraint(pm.model, vm_ac == branch["sadra_Vtar"])
        else
            @warn "SADRA PF: VSC branch $i has type_ac=$type_ac (no AC-side " *
                  "control). Its AC degree of freedom is unpinned; the PF " *
                  "may be under-determined."
        end
    end

    # --- DC branches: resistive Ohm's law (no limits in PF) ---
    for i in dc_branch_ids
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
    end
end

# ====================================================================
# Solve entry points
# ====================================================================

"""
    solve_sadra_pf(file::String, optimizer; setting=Dict())

Parse a PowerModelsACDC-format `.m` case, apply the SADRA transform,
and solve the AC/DC power flow (feasibility problem, no objective).
"""
function solve_sadra_pf(file::String, optimizer; setting=Dict{String,Any}())
    data = _PM.parse_file(file)
    return solve_sadra_pf(data, optimizer; setting=setting)
end

"""
    solve_sadra_fubm_pf(file::String, optimizer; setting=Dict())

As above for a MATPOWER-FUBM `.m` case (e.g. the 1354-bus PEGASE case).
"""
function solve_sadra_fubm_pf(file::String, optimizer; setting=Dict{String,Any}())
    data = parse_fubm_acdc(file)
    return solve_sadra_pf(data, optimizer; setting=setting)
end

"""
    solve_sadra_pf(data::Dict, optimizer; setting=Dict())

As above, but takes an already-parsed PowerModels data dictionary.
"""
function solve_sadra_pf(data::Dict{String,Any}, optimizer; setting=Dict{String,Any}())
    sadra_transform!(data)
    return _PM.solve_model(
        data,
        _PM.ACPPowerModel,
        optimizer,
        build_sadra_pf;
        setting = setting
    )
end
