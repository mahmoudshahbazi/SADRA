# base.jl
# Builds SADRA-specific lookup sets directly from the PowerModels ref,
# called at the start of build_sadra_opf after PM has populated ref[:branch] etc.
#
# We do NOT use ref_extensions because InfrastructureModels calls those
# before ref_add_core! has populated :branch, :gen etc.
# Instead, build_sadra_opf calls sadra_build_ref!(pm) as its first step.

import PowerModels as _PM

function sadra_build_ref!(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)

    ref = _PM.ref(pm, nw)

    # VSC branches: flagged sadra_vsc == true
    ref[:sadra_vsc_branch] = Dict(
        i => branch
        for (i, branch) in ref[:branch]
        if get(branch, "sadra_vsc", false)
    )

    # DC branches: flagged sadra_dc_branch == true
    ref[:sadra_dc_branch] = Dict(
        i => branch
        for (i, branch) in ref[:branch]
        if get(branch, "sadra_dc_branch", false)
    )

    # Dummy generators: flagged sadra_dummy_gen == true
    ref[:sadra_dummy_gen] = Dict(
        i => gen
        for (i, gen) in ref[:gen]
        if get(gen, "sadra_dummy_gen", false)
    )

    # DC buses: flagged sadra_dc_bus == true
    ref[:sadra_dc_bus] = Dict(
        i => bus
        for (i, bus) in ref[:bus]
        if get(bus, "sadra_dc_bus", false)
    )

    # Per-converter lookup: conv_i -> (branch_i, gen_i, dc_bus_i, ac_bus_i)
    sadra_conv = Dict{Int, Dict{Symbol,Int}}()
    for (br_i, branch) in ref[:sadra_vsc_branch]
        conv_i   = branch["sadra_conv_i"]
        dc_bus_i = branch["sadra_dc_bus_i"]
        ac_bus_i = branch["sadra_ac_bus_i"]
        gen_i    = 0
        for (g_i, gen) in ref[:sadra_dummy_gen]
            if gen["sadra_conv_i"] == conv_i
                gen_i = g_i
                break
            end
        end
        sadra_conv[conv_i] = Dict{Symbol,Int}(
            :branch_i => br_i,
            :gen_i    => gen_i,
            :dc_bus_i => dc_bus_i,
            :ac_bus_i => ac_bus_i,
        )
    end
    ref[:sadra_conv] = sadra_conv

    return nothing
end
