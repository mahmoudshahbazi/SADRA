# constraint_template.jl
# Data-extraction layer. Each function pulls bus indices, impedances and
# parameters from the PowerModels ref, then calls the matching math
# function in constraints.jl with clean numeric arguments.
#
#   A. constraint_vsc_ohms_from   — VSC branch Ohm's law, from (DC) side
#   B. constraint_vsc_ohms_to     — VSC branch Ohm's law, to (AC) side
#   C. constraint_vsc_current     — AC-side current magnitude definition
#   D. constraint_vsc_losses      — converter loss equation (dummy gen)
#   E. constraint_vsc_zero_qf     — q_fr = 0 (keeps DC grid reactive-free)
#   F. constraint_vsc_p_setpoint  — optional AC active power control
#      constraint_vsc_q_setpoint  — optional AC reactive power control

import PowerModels as _PM

# Helper: series conductance/susceptance from branch r, x
@inline function _series_gb(branch)
    r = branch["br_r"]; x = branch["br_x"]
    denom = r^2 + x^2
    g = (denom > 0) ? r/denom : 0.0
    b = (denom > 0) ? -x/denom : 0.0
    return g, b
end

# A. VSC Ohm's law from-side
function constraint_vsc_ohms_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus)
    t_idx  = (i, t_bus, f_bus)
    g, b   = _series_gb(branch)
    constraint_vsc_ohms_from(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, branch["b_fr"])
end

# B. VSC Ohm's law to-side
function constraint_vsc_ohms_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus)
    t_idx  = (i, t_bus, f_bus)
    g, b   = _series_gb(branch)
    constraint_vsc_ohms_to(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, branch["b_to"])
end

# C. VSC current definition
function constraint_vsc_current(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    t_bus  = branch["t_bus"]
    f_bus  = branch["f_bus"]
    t_idx  = (i, t_bus, f_bus)
    conv_i = branch["sadra_conv_i"]
    gen_i  = _PM.ref(pm, nw, :sadra_conv)[conv_i][:gen_i]
    constraint_vsc_current(pm, nw, i, t_bus, t_idx, gen_i)
end

# D. VSC loss equation
function constraint_vsc_losses(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    conv_i = branch["sadra_conv_i"]
    gen_i  = _PM.ref(pm, nw, :sadra_conv)[conv_i][:gen_i]
    constraint_vsc_losses(pm, nw, i, gen_i,
        branch["sadra_LossA"], branch["sadra_LossB"], branch["sadra_LossC"])
end

# E. VSC from-side reactive = 0
function constraint_vsc_zero_qf(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus)
    constraint_vsc_zero_qf(pm, nw, i, f_idx)
end

# F. Optional control setpoints
# P control: p_fr == Pset + Pg (AIMMS PfShiftControl relation; +Pg = loss term).
function constraint_vsc_p_setpoint(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus)
    conv_i = branch["sadra_conv_i"]
    gen_i  = _PM.ref(pm, nw, :sadra_conv)[conv_i][:gen_i]
    constraint_vsc_p_setpoint(pm, nw, i, f_idx, gen_i, branch["sadra_Pset"])
end

# Q control (eq 23): reactive power on the to (AC) side fixed to Qset.
function constraint_vsc_q_setpoint(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    t_idx  = (i, t_bus, f_bus)
    constraint_vsc_q_setpoint(pm, nw, i, t_idx, branch["sadra_Qset"])
end

# DC voltage pin (type_dc=2): vm at the DC bus fixed to Vdcset.
function constraint_vsc_vdc_setpoint(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    dc_bus = branch["f_bus"]
    constraint_vsc_vdc_setpoint(pm, nw, dc_bus, branch["sadra_Vdcset"])
end

# Droop (type_dc=3, eq 24): P_f = Pdcset - kd*(vm_dc - Vdcref) + Pg.
function constraint_vsc_droop(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus)
    conv_i = branch["sadra_conv_i"]
    gen_i  = _PM.ref(pm, nw, :sadra_conv)[conv_i][:gen_i]
    vdc_ref = get(branch, "sadra_droop_vref", branch["sadra_Vdcset"])
    constraint_vsc_droop(pm, nw, i, f_idx, f_bus, gen_i,
        branch["sadra_droop_kd"], branch["sadra_Pdcset"], vdc_ref)
end

function constraint_vsc_qg_link(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (i, f_bus, t_bus)
    conv_i = branch["sadra_conv_i"]
    gen_i  = _PM.ref(pm, nw, :sadra_conv)[conv_i][:gen_i]
    constraint_vsc_qg_link(pm, nw, i, f_idx, gen_i)
end