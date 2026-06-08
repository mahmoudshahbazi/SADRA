# xfmr_control.jl
# Controlled transformers (PST and CTT), per SADRA paper Fig. 2(e)/2(f) and
# Table IV. These are ORDINARY AC branches (no DC bus, no dummy generator,
# no loss model, no current variable) — distinct from VSCs — so they get
# their own variables and constraints rather than reusing the VSC machinery.
#
#   PST (phase shifter): controls from-side active power P_f to a setpoint.
#                        Variable = phase shift phi_x; tap magnitude fixed = 1.
#   CTT (controlled tap): controls to-side voltage magnitude V_t to a setpoint.
#                        Variable = tap magnitude ma_x; phase shift fixed = 0.
#
# Both use the standard Matpower variable-tap branch equations (paper eq 7-8),
# identical in form to the VSC Ohm's law but with their own variables.

import PowerModels as _PM
using JuMP

# ------------------------------------------------------------------
# Ref set: controlled transformers, flagged sadra_xfmr_ctrl on the branch.
# Called from sadra_build_ref! (added there).
# ------------------------------------------------------------------
function sadra_build_xfmr_ref!(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)
    ref = _PM.ref(pm, nw)
    ref[:sadra_ctrl_xfmr] = Dict(
        i => branch
        for (i, branch) in ref[:branch]
        if haskey(branch, "sadra_xfmr_ctrl")
    )
    return nothing
end

# ------------------------------------------------------------------
# Variables: ma_x (tap) and phi_x (shift) for controlled transformers.
# Bounds taken from the branch data where present, else sensible defaults.
# For a PST, ma_x is fixed to 1 (only phi varies); for a CTT, phi_x is fixed
# to 0 (only ma varies). We declare both for all and pin the inactive one.
# ------------------------------------------------------------------
function variable_xfmr_control(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)
    ids = _PM.ids(pm, nw, :sadra_ctrl_xfmr)

    ma_x = _PM.var(pm, nw)[:ma_x] = JuMP.@variable(pm.model,
        [i in ids], base_name="$(nw)_ma_x", start=1.0)
    phi_x = _PM.var(pm, nw)[:phi_x] = JuMP.@variable(pm.model,
        [i in ids], base_name="$(nw)_phi_x", start=0.0)

    for (i, branch) in _PM.ref(pm, nw, :sadra_ctrl_xfmr)
        kind = branch["sadra_xfmr_ctrl"]
        if kind == "pst"
            # phase shifter: tap fixed at 1, phase free within +/- pi/2
            JuMP.fix(ma_x[i], 1.0; force=true)
            JuMP.set_lower_bound(phi_x[i], -pi/2)
            JuMP.set_upper_bound(phi_x[i],  pi/2)
        elseif kind == "ctt"
            # controlled tap: phase fixed at 0, tap free within bounds
            JuMP.fix(phi_x[i], 0.0; force=true)
            JuMP.set_lower_bound(ma_x[i], get(branch, "sadra_ma_min", 0.9))
            JuMP.set_upper_bound(ma_x[i], get(branch, "sadra_ma_max", 1.1))
        end
    end
    return nothing
end

# ------------------------------------------------------------------
# Ohm's law for a controlled transformer (variable tap/shift form).
# Same equations as the VSC Ohm's law, but reads ma_x/phi_x and there is
# no from-side shunt subtlety (uses branch b_fr/b_to as parsed).
# ------------------------------------------------------------------
function constraint_xfmr_ohms(pm::_PM.AbstractACPModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]; t_bus = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus); t_idx = (i, t_bus, f_bus)
    r = branch["br_r"]; x = branch["br_x"]
    denom = r^2 + x^2
    g = denom > 0 ? r/denom : 0.0
    b = denom > 0 ? -x/denom : 0.0
    b_fr = get(branch, "b_fr", 0.0); b_to = get(branch, "b_to", 0.0)

    p_fr = _PM.var(pm, nw, :p)[f_idx]; q_fr = _PM.var(pm, nw, :q)[f_idx]
    p_to = _PM.var(pm, nw, :p)[t_idx]; q_to = _PM.var(pm, nw, :q)[t_idx]
    vm_f = _PM.var(pm, nw, :vm, f_bus); vm_t = _PM.var(pm, nw, :vm, t_bus)
    va_f = _PM.var(pm, nw, :va, f_bus); va_t = _PM.var(pm, nw, :va, t_bus)
    ma   = _PM.var(pm, nw, :ma_x, i);   phi  = _PM.var(pm, nw, :phi_x, i)

    JuMP.@NLconstraint(pm.model, p_fr == g/ma^2*vm_f^2
        - (1/ma)*vm_f*vm_t*(g*cos(va_f-va_t-phi) + b*sin(va_f-va_t-phi)))
    JuMP.@NLconstraint(pm.model, q_fr == -(b+b_fr)/ma^2*vm_f^2
        - (1/ma)*vm_f*vm_t*(g*sin(va_f-va_t-phi) - b*cos(va_f-va_t-phi)))
    JuMP.@NLconstraint(pm.model, p_to == g*vm_t^2
        - (1/ma)*vm_f*vm_t*(g*cos(va_t-va_f+phi) + b*sin(va_t-va_f+phi)))
    JuMP.@NLconstraint(pm.model, q_to == -(b+b_to)*vm_t^2
        - (1/ma)*vm_f*vm_t*(g*sin(va_t-va_f+phi) - b*cos(va_t-va_f+phi)))
end

# ------------------------------------------------------------------
# PST control: from-side active power fixed to setpoint (paper Table IV).
# ------------------------------------------------------------------
function constraint_pst_pf(pm::_PM.AbstractACPModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]; t_bus = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus)
    p_fr   = _PM.var(pm, nw, :p)[f_idx]
    JuMP.@constraint(pm.model, p_fr == branch["sadra_pst_pf"])
end

# ------------------------------------------------------------------
# CTT control: to-side voltage magnitude fixed to setpoint (paper Table IV).
# ------------------------------------------------------------------
function constraint_ctt_vt(pm::_PM.AbstractACPModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    t_bus  = branch["t_bus"]
    vm_t   = _PM.var(pm, nw, :vm, t_bus)
    JuMP.@constraint(pm.model, vm_t == branch["sadra_ctt_vt"])
end
