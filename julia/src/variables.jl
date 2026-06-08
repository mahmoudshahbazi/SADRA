# variables.jl
# Declares the optimisation variables that SADRA adds on top of the
# standard PowerModels AC variables (vm, va, pg, qg, branch p/q).
#
# New variable families (one per VSC converter branch):
#   ma       — VSC tap ratio
#   phi      — VSC phase shift angle (radians)
#   iconv_ac — AC-side current magnitude (>= 0), used in the loss equation

import PowerModels as _PM
using JuMP

# ------------------------------------------------------------------
# ma: VSC tap ratio. Bounds from sadra_ma_min / sadra_ma_max
# (set from conv Vmmin / Vmmax in data.jl).
# ------------------------------------------------------------------
function variable_ma(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true)
    ma = _PM.var(pm, nw)[:ma] = JuMP.@variable(
        pm.model,
        [i in _PM.ids(pm, nw, :sadra_vsc_branch)],
        base_name = "$(nw)_ma",
        start = 0.95   # warm start: VSC taps converge to ~0.93-0.98 in practice
    )
    if bounded
        for (i, branch) in _PM.ref(pm, nw, :sadra_vsc_branch)
            JuMP.set_lower_bound(ma[i], branch["sadra_ma_min"])
            JuMP.set_upper_bound(ma[i], branch["sadra_ma_max"])
        end
    end
    return ma
end

# ------------------------------------------------------------------
# phi: VSC phase shift angle (radians).
# ------------------------------------------------------------------
function variable_phi(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true)
    phi = _PM.var(pm, nw)[:phi] = JuMP.@variable(
        pm.model,
        [i in _PM.ids(pm, nw, :sadra_vsc_branch)],
        base_name = "$(nw)_phi",
        start = 0.0
    )
    if bounded
        for (i, branch) in _PM.ref(pm, nw, :sadra_vsc_branch)
            JuMP.set_lower_bound(phi[i], branch["sadra_phi_min"])
            JuMP.set_upper_bound(phi[i], branch["sadra_phi_max"])
        end
    end
    return phi
end

# ------------------------------------------------------------------
# iconv_ac: AC-side current magnitude (>= 0). Defined by constraint C:
#   p_to^2 + q_to^2 = vm_t^2 * iconv_ac^2
# In the AIMMS reference, the converter current (ABSIt) is declared
# Range:free — it has NO upper bound; it is purely the loss-equation input.
# An artificial cap here makes high-power setpoints (e.g. Pf=-500MW ~ 5pu
# current) infeasible. Real limits come from branch thermal rate + voltages.
# We keep only the physical lower bound of 0 (it is a magnitude).
# ------------------------------------------------------------------
function variable_iconv_ac(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true)
    iconv_ac = _PM.var(pm, nw)[:iconv_ac] = JuMP.@variable(
        pm.model,
        [i in _PM.ids(pm, nw, :sadra_vsc_branch)],
        base_name = "$(nw)_iconv_ac",
        lower_bound = 0.0,
        start = 1.0
    )
    # Warm start: a current variable starting at 0 (or 1) is far from the
    # solution for high-power converters (e.g. VSC4 needs i~4.7), and the
    # nonlinear defining constraint p^2+q^2=vm^2*i^2 then takes many iterations
    # to climb. Seed iconv_ac from the converter's setpoint power |Pset|/vm
    # where a setpoint exists, else leave at 1.0.
    for (i, branch) in _PM.ref(pm, nw, :sadra_vsc_branch)
        pset = get(branch, "sadra_Pset", 0.0)
        if pset != 0.0
            JuMP.set_start_value(iconv_ac[i], abs(pset))  # vm~1 => i~|P|
        end
    end
    # Cap current only for PMACDC input (v0.1 behaviour). For FUBM input the
    # converter current is free (AIMMS ABSIt is Range:free); an artificial cap
    # makes high-power setpoints infeasible. The FUBM ingest sets sadra_loss_pu.
    if bounded
        for (i, branch) in _PM.ref(pm, nw, :sadra_vsc_branch)
            if !get(branch, "sadra_loss_pu", false)
                JuMP.set_upper_bound(iconv_ac[i], branch["sadra_imax"])
            end
        end
    end
    return iconv_ac
end

# ------------------------------------------------------------------
# Declare all SADRA variables at once.
# ------------------------------------------------------------------
function variable_sadra_converter(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true)
    variable_ma(pm; nw=nw, bounded=bounded)
    variable_phi(pm; nw=nw, bounded=bounded)
    variable_iconv_ac(pm; nw=nw, bounded=bounded)
end
