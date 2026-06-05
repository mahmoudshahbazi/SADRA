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
        start = 1.0
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
# Upper bound is the converter rated current (sadra_imax).
# ------------------------------------------------------------------
function variable_iconv_ac(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true)
    iconv_ac = _PM.var(pm, nw)[:iconv_ac] = JuMP.@variable(
        pm.model,
        [i in _PM.ids(pm, nw, :sadra_vsc_branch)],
        base_name = "$(nw)_iconv_ac",
        lower_bound = 0.0,
        start = 0.0
    )
    if bounded
        for (i, branch) in _PM.ref(pm, nw, :sadra_vsc_branch)
            JuMP.set_upper_bound(iconv_ac[i], branch["sadra_imax"])
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
