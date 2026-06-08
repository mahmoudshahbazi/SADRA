# constraints.jl
# NLP constraint math for the ACP (polar AC) formulation of SADRA.
#
# These functions receive only JuMP variables and numeric parameters —
# no ref/data dict lookups (those happen in constraint_template.jl).
#
# The VSC branch is a standard Matpower branch (eq 7-8 of the paper) with
# variable tap ratio (ma) and phase shift (phi). Its series impedance is
# rtf+j*xtf from the converter data. The dummy generator at the DC bus
# models converter losses (eq 4) and provides the reactive balance that
# keeps the DC grid reactive-power-free (q_fr = 0, eq 3).

import PowerModels as _PM
using JuMP

# ====================================================================
# A. VSC Ohm's law — from-side (DC bus)
# ====================================================================
# Matpower branch equations (eq 7-8) with variable tap ma and phase phi.
# Series admittance g+jb from rtf+j*xtf. b_fr is the from-side shunt
# susceptance (zero for VSC; filter is on the AC/to side).

function constraint_vsc_ohms_from(pm::_PM.AbstractACPModel, n::Int, i::Int,
        f_bus, t_bus, f_idx, t_idx, g, b, b_fr)
    p_fr = _PM.var(pm, n, :p)[f_idx]
    q_fr = _PM.var(pm, n, :q)[f_idx]
    vm_f = _PM.var(pm, n, :vm, f_bus)
    vm_t = _PM.var(pm, n, :vm, t_bus)
    va_f = _PM.var(pm, n, :va, f_bus)
    va_t = _PM.var(pm, n, :va, t_bus)
    ma   = _PM.var(pm, n, :ma, i)
    phi  = _PM.var(pm, n, :phi, i)

    JuMP.@NLconstraint(pm.model,
        p_fr == g / ma^2 * vm_f^2
              - (1/ma) * vm_f * vm_t * (g * cos(va_f - va_t - phi)
                                       + b * sin(va_f - va_t - phi))
    )
    JuMP.@NLconstraint(pm.model,
        q_fr == -(b + b_fr) / ma^2 * vm_f^2
              - (1/ma) * vm_f * vm_t * (g * sin(va_f - va_t - phi)
                                       - b * cos(va_f - va_t - phi))
    )
end

# ====================================================================
# B. VSC Ohm's law — to-side (AC bus)
# ====================================================================
function constraint_vsc_ohms_to(pm::_PM.AbstractACPModel, n::Int, i::Int,
        f_bus, t_bus, f_idx, t_idx, g, b, b_to)
    p_to = _PM.var(pm, n, :p)[t_idx]
    q_to = _PM.var(pm, n, :q)[t_idx]
    vm_f = _PM.var(pm, n, :vm, f_bus)
    vm_t = _PM.var(pm, n, :vm, t_bus)
    va_f = _PM.var(pm, n, :va, f_bus)
    va_t = _PM.var(pm, n, :va, t_bus)
    ma   = _PM.var(pm, n, :ma, i)
    phi  = _PM.var(pm, n, :phi, i)

    JuMP.@NLconstraint(pm.model,
        p_to == g * vm_t^2
              - (1/ma) * vm_f * vm_t * (g * cos(va_t - va_f + phi)
                                       + b * sin(va_t - va_f + phi))
    )
    JuMP.@NLconstraint(pm.model,
        q_to == -(b + b_to) * vm_t^2
              - (1/ma) * vm_f * vm_t * (g * sin(va_t - va_f + phi)
                                       - b * cos(va_t - va_f + phi))
    )
end

# ====================================================================
# C. VSC current definition
# ====================================================================
# Current magnitude at the AC (to) terminal, used in the loss equation.
# Matpower convention |it| = |St|/|vt|, rearranged to avoid division:
#   p_to^2 + q_to^2 = vm_t^2 * iconv_ac^2

function constraint_vsc_current(pm::_PM.AbstractACPModel, n::Int, i::Int,
        t_bus, t_idx, gen_i)
    p_to     = _PM.var(pm, n, :p)[t_idx]
    q_to     = _PM.var(pm, n, :q)[t_idx]
    vm_t     = _PM.var(pm, n, :vm, t_bus)
    iconv_ac = _PM.var(pm, n, :iconv_ac, i)

    JuMP.@NLconstraint(pm.model,
        p_to^2 + q_to^2 == vm_t^2 * iconv_ac^2
    )
end

# ====================================================================
# D. VSC loss equation — eq (4) of paper
# ====================================================================
# P_gdc = -gamma*iconv^2 - beta*iconv - alpha
# Dummy generator active power equals the negative of the converter losses.

function constraint_vsc_losses(pm::_PM.AbstractACPModel, n::Int, i::Int,
        gen_i, a, b, c)
    pg       = _PM.var(pm, n, :pg, gen_i)
    iconv_ac = _PM.var(pm, n, :iconv_ac, i)

    JuMP.@NLconstraint(pm.model,
        pg == -c * iconv_ac^2 - b * iconv_ac - a
    )
end

# ====================================================================
# E. VSC from-side reactive = 0 — eq (3) of paper
# ====================================================================
# Forcing q_fr = 0 at the DC (from) side keeps the DC grid reactive-free.
# The dummy generator's qg is left free and absorbs the reactive balance
# at the DC bus automatically via PowerModels' bus power balance.

function constraint_vsc_zero_qf(pm::_PM.AbstractACPModel, n::Int,
        i::Int, f_idx)
    q_fr = _PM.var(pm, n, :q)[f_idx]
    JuMP.@constraint(pm.model, q_fr == 0)
end

# ====================================================================
# F. VSC control constraints — eq (22-24) of paper (optional)
# ====================================================================
# P control: AIMMS pins the CONTROL setpoint, not the raw branch flow:
#   Pf*Sbase - PfShiftControl*Sbase - Pg*Sbase = 0  ->  p_fr = p_set + pg
# i.e. the branch from-side flow equals the setpoint plus the dummy-gen power
# (the converter loss). Pinning p_fr == p_set directly (without +pg) is
# infeasible: the true feasible p_fr differs from the setpoint by exactly Pg
# (verified against AIMMS: VSC4 setpoint -5.0 -> Pf -5.0022 with Pg -0.0022).
function constraint_vsc_p_setpoint(pm::_PM.AbstractACPModel, n::Int, i::Int,
        f_idx, gen_i, p_set)
    p_fr = _PM.var(pm, n, :p)[f_idx]
    pg   = _PM.var(pm, n, :pg, gen_i)
    JuMP.@constraint(pm.model, p_fr == p_set + pg)
end

# Q control (eq 23): reactive power on the to (AC) side fixed.
# p_to/q_to > 0 is power from the AC bus into the branch; PMACDC Q_g is
# positive into the AC grid, hence q_to == -Q_set.
function constraint_vsc_q_setpoint(pm::_PM.AbstractACPModel, n::Int, i::Int,
        t_idx, q_set)
    q_to = _PM.var(pm, n, :q)[t_idx]
    JuMP.@constraint(pm.model, q_to == -q_set)
end

# DC voltage control (type_dc=2): pin DC bus vm to Vdcset.
function constraint_vsc_vdc_setpoint(pm::_PM.AbstractACPModel, n::Int,
        dc_bus, vdc_set)
    vm_dc = _PM.var(pm, n, :vm, dc_bus)
    JuMP.@constraint(pm.model, vm_dc == vdc_set)
end

# Droop control (type_dc=3, eq 24): matches AIMMS PfVdcDroop:
#   Pf = Pref - kd*(vm_dc - Vdcref) + Pg
# The +Pg (dummy-gen loss power) is essential: p_fr is the DC-side branch
# flow, which differs from the net DC injection by the converter loss. The
# AIMMS reference includes it, so we must too or the droop pins the wrong
# quantity. Pg < 0 (a loss), so it shifts p_fr by the loss magnitude.
function constraint_vsc_droop(pm::_PM.AbstractACPModel, n::Int, i::Int,
        f_idx, dc_bus, gen_i, kd, p_ref, vdc_ref)
    p_fr  = _PM.var(pm, n, :p)[f_idx]
    vm_dc = _PM.var(pm, n, :vm, dc_bus)
    pg    = _PM.var(pm, n, :pg, gen_i)
    JuMP.@constraint(pm.model, p_fr == p_ref - kd * (vm_dc - vdc_ref) + pg)
end


function constraint_vsc_qg_link(pm::_PM.AbstractACPModel, n::Int, i::Int,
        f_idx, gen_i)
    q_fr = _PM.var(pm, n, :q)[f_idx]
    qg   = _PM.var(pm, n, :qg, gen_i)
    JuMP.@constraint(pm.model, qg == q_fr)
end