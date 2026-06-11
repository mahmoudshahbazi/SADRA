# data.jl
# Transforms a PowerModelsACDC-parsed data dict into a SADRA-augmented
# pure-AC data dict. No JuMP, no PowerModels internals — pure data work.
#
# After this function runs, `data` contains:
#   - All original AC buses, branches, generators (unchanged)
#   - New AC buses for each DC bus (indexed from ac_bus_offset+1)
#   - New AC branches for each DC line (r from data, x=1e-6, b=0)
#   - New AC branches for each VSC (series Z conditional on transformer flag)
#   - New dummy generators at each DC bus (cost=0, pg/qg free)
#   - data["sadra"] dict with index maps and metadata for later use
#
# Convention (from SADRA paper and confirmed from PMACDC output):
#   VSC branch: from-bus = DC bus, to-bus = AC bus
#   Loss current It is at the AC (to) side

function sadra_transform!(data::Dict{String,Any})

    baseMVA = data["baseMVA"]

    # ------------------------------------------------------------------
    # 1. Compute safe index offsets so new entries don't collide
    # ------------------------------------------------------------------
    max_ac_bus   = maximum(parse(Int, k) for k in keys(data["bus"]))
    max_branch   = maximum(parse(Int, k) for k in keys(data["branch"]))
    max_gen      = maximum(parse(Int, k) for k in keys(data["gen"]))

    ac_bus_offset  = max_ac_bus    # DC bus i -> AC bus index (ac_bus_offset + i)
    dc_branch_offset = max_branch  # DC branch i -> AC branch index (dc_branch_offset + i)
    vsc_branch_offset = max_branch + length(data["branchdc"])  # VSC i -> AC branch
    dummy_gen_offset  = max_gen    # Conv i -> dummy gen index (dummy_gen_offset + i)

    # ------------------------------------------------------------------
    # 2. Index maps (stored for use in constraints and variables)
    # ------------------------------------------------------------------
    # dc_bus_i (1-based DC bus number) -> AC bus index in data["bus"]
    dc_to_ac_bus = Dict{Int,Int}()
    for (_, busdc) in data["busdc"]
        dc_i = busdc["busdc_i"]
        dc_to_ac_bus[dc_i] = ac_bus_offset + dc_i
    end

    # conv index -> AC branch index for the VSC branch
    conv_to_branch = Dict{Int,Int}()
    for (_, conv) in data["convdc"]
        c_i = conv["index"]
        conv_to_branch[c_i] = vsc_branch_offset + c_i
    end

    # conv index -> dummy gen index
    conv_to_gen = Dict{Int,Int}()
    for (_, conv) in data["convdc"]
        c_i = conv["index"]
        conv_to_gen[c_i] = dummy_gen_offset + c_i
    end

    # ------------------------------------------------------------------
    # 3. Add DC buses as AC buses
    # ------------------------------------------------------------------
    max_load = isempty(data["load"]) ? 0 : maximum(parse(Int, k) for k in keys(data["load"]))

    for (_, busdc) in data["busdc"]
        dc_i  = busdc["busdc_i"]
        new_i = dc_to_ac_bus[dc_i]

        # DC buses are PQ-type in the OPF (voltage controlled by VSC tap ma)
        # bus_type=1 (PQ). The VSC dummy gen will handle P/Q injection.
        # Check if any converter at this DC bus has type_dc=2 (voltage control)
        # If so, pin vm to Vdcset via tight bounds
        vdc_min = busdc["Vdcmin"]
        vdc_max = busdc["Vdcmax"]
        #= for (_, conv) in data["convdc"]
            if conv["busdc_i"] == dc_i && get(conv, "type_dc", 1) == 2
                vdcset = get(conv, "Vdcset", busdc["Vdc"])
                vdc_min = vdcset
                vdc_max = vdcset
                break
            end
        end =#

        data["bus"]["$new_i"] = Dict{String,Any}(
            "index"      => new_i,
            "bus_i"      => new_i,
            "bus_type"   => 1,
            "vm"         => 1.01,    # warm start near typical DC operating voltage
            "va"         => -0.06,   # warm start near converged common DC angle
            "vmax"       => vdc_max,
            "vmin"       => vdc_min,
            "base_kv"    => busdc["basekVdc"],
            "zone"       => 1,
            "area"       => 1,
            "sadra_dc_bus" => true,
            "sadra_dc_i"   => dc_i,
        )

        # Add Pdc as a load at the DC bus (if nonzero)
        # Pdc > 0 = power consumed from DC bus, Pdc < 0 = power injected
        pdc_raw = get(busdc, "Pdc", 0.0)
        if pdc_raw != 0.0
            max_load += 1
            data["load"]["$max_load"] = Dict{String,Any}(
                "index"      => max_load,
                "load_bus"   => new_i,
                "pd"         => pdc_raw / baseMVA,
                "qd"         => 0.0,
                "status"     => 1,
                "sadra_dc_load" => true,
            )
        end
    end

    # ------------------------------------------------------------------
    # 4. Add DC lines as AC branches
    # ------------------------------------------------------------------
    # Q=0 is enforced by explicit constraints in constraints.jl, not here.
    # We set x to a small nonzero value to keep Jacobian well-conditioned.
    DC_X_STUB = 0.0   # DC lines are purely resistive; Q=0 enforced explicitly

    for (_, brdc) in data["branchdc"]
        br_i  = brdc["index"]
        new_i = dc_branch_offset + br_i

        f_bus_ac = dc_to_ac_bus[brdc["fbusdc"]]
        t_bus_ac = dc_to_ac_bus[brdc["tbusdc"]]

        data["branch"]["$new_i"] = Dict{String,Any}(
            "index"    => new_i,
            "f_bus"    => f_bus_ac,
            "t_bus" => t_bus_ac,
            # new — MatACDC/PMACDC bipolar convention: P = dcpol·(1/r)·vf·(vf−vt),
            # so the equivalent AC branch resistance is r/dcpol. FUBM-parsed data has
            # no "dcpol" key and its r is already in the equivalent-AC convention -> /1.
            "br_r" => brdc["r"] / get(data, "dcpol", 1),
            "br_x" => DC_X_STUB,
            "br_b" => 0.0,
            "g_fr" => 0.0,
            "g_to" => 0.0,
            "b_fr"     => 0.0,
            "b_to"     => 0.0,
            "tap"      => 1.0,
            "shift"    => 0.0,
            "br_status"=> 1,
            "angmin"   => -pi,   # angle bounds wide open — Q=0 via explicit constraint
            "angmax"   =>  pi,
            "rate_a"   => brdc["rateA"],
            "rate_b"   => get(brdc, "rateB", brdc["rateA"]),
            "rate_c"   => get(brdc, "rateC", brdc["rateA"]),
            "transformer" => false,
            # Custom flags
            "sadra_dc_branch" => true,
            "sadra_dc_br_i"   => br_i,
        )
    end

    # ------------------------------------------------------------------
    # 5. Add VSC converters as AC branches (variable tap/shift)
    # ------------------------------------------------------------------
    # The branch represents the VSC series path:
    #   from-bus = DC bus (AC index), to-bus = AC bus
    #   Series impedance: rtf+j*xtf if transformer=1, else 0+j*DC_X_STUB
    # The tap ratio `ma` and phase shift `phi` are optimisation variables
    # declared in variables.jl. Here we store the data; the branch equations
    # are written in constraints.jl (not using PM's standard branch model).

    for (_, conv) in data["convdc"]
        c_i   = conv["index"]
        new_i = vsc_branch_offset + c_i

        f_bus_ac = dc_to_ac_bus[conv["busdc_i"]]   # DC side (from)
        t_bus_ac = conv["busac_i"]                   # AC side (to)

        # Always use rtf/xtf from data regardless of transformer flag.
        # The flag in PMACDC controls whether to model the transformer as a
        # separate internal π-circuit — it does NOT mean zero impedance.
        # The physical impedance values are always present in the data.
        # Fallback to DC_X_STUB only if both are truly zero (degenerate data).
        r_series = get(conv, "rtf", 0.0)
        x_series = get(conv, "xtf", 0.0)
        if r_series == 0.0 && x_series == 0.0
            x_series = DC_X_STUB
        end

        # Filter susceptance (shunt on AC side), honour filter flag
        b_filter = (conv["filter"] == 1) ? conv["bf"] : 0.0

        # Imax from data (PMACDC recalculates this — we use the stored value)
        # rateA = Imax * Vmmax (apparent power limit at max voltage)
        # We store Imax directly; the thermal limit constraint uses it
        imax   = get(conv, "Imax", 1.1)
        rate_a = get(conv, "Pacmax", 100.0) / baseMVA

        data["branch"]["$new_i"] = Dict{String,Any}(
            "index"    => new_i,
            "f_bus"    => f_bus_ac,   # DC bus
            "t_bus"    => t_bus_ac,   # AC bus
            "br_r"     => r_series,
            "br_x"     => x_series,
            "br_b"     => 0.0,        # no shunt on series element
            "g_fr"     => 0.0,
            "g_to"     => 0.0,
            "b_fr"     => 0.0,
            "b_to"     => b_filter,   # filter on AC (to) side
            "tap"      => 1.0,        # nominal; actual tap is variable ma
            "shift"    => 0.0,        # nominal; actual shift is variable phi
            "br_status"=> 1,
            "angmin"   => -pi,
            "angmax"   =>  pi,
            "rate_a"   => rate_a,
            "rate_b"   => rate_a,
            "rate_c"   => rate_a,
            "transformer" => true,    # tell PM this is a transformer branch
            # SADRA-specific metadata
            "sadra_vsc"    => true,
            "sadra_pset_side" => "grid",   # PMACDC convention: P_g defined at the PCC (grid side)
            "sadra_conv_i" => c_i,
            "sadra_imax"   => imax,
            "sadra_ma_min" => conv["Vmmin"],
            "sadra_ma_max" => conv["Vmmax"],
            "sadra_phi_min"=> -pi/2,
            "sadra_phi_max"=>  pi/2,
            # Control mode (from SADRA paper Table I / eq 22-24).
            # type_dc: 1=P_f setpoint, 2=Vdc pin, 3=droop
            # type_ac: 1=Q setpoint,   2=Vac pin
            "sadra_type_dc"  => conv["type_dc"],
            "sadra_type_ac"  => conv["type_ac"],
            "sadra_Pset"     => get(conv, "P_g", 0.0) / baseMVA,
            "sadra_Qset"     => get(conv, "Q_g", 0.0) / baseMVA,
            "sadra_Vdcset"   => get(conv, "Vdcset", 1.0),
            "sadra_Vtar"     => get(conv, "Vtar",   1.0),
            # Droop (type_dc=3): P_f = Pref - kd*(vm_dc - Vdcset). eq 24.
            # PMACDC stores droop slope as `droop`; Pdcset/Vdcset are the refs.
            "sadra_droop_kd" => get(conv, "droop", 0.0),
            "sadra_Pdcset"   => get(conv, "Pdcset", get(conv, "P_g", 0.0)) / baseMVA,
            "sadra_droop_vref" => get(conv, "droop_vref", get(conv, "Vdcset", 1.0)),
            # Loss coefficients. Two input conventions:
            #  * PMACDC: LossA [MW], LossB [MW/kA], LossC [MW/kA²] -> kA-normalise.
            #  * FUBM:   ALPHA1/2/3 already p.u. on p.u. current -> no normalisation.
            #    The FUBM ingest sets conv["sadra_loss_pu"]=true.
            "sadra_LossA"  => conv["LossA"] / baseMVA,
            "sadra_LossB"  => get(conv, "sadra_loss_pu", false) ?
                                  conv["LossB"] :
                                  conv["LossB"] * (baseMVA / (sqrt(3) * conv["basekVac"])) / baseMVA,
            "sadra_LossC"  => get(conv, "sadra_loss_pu", false) ?
                                  conv["LossCinv"] :
                                  conv["LossCinv"] * (baseMVA / (sqrt(3) * conv["basekVac"]))^2 / baseMVA,
            "sadra_dc_bus_i" => f_bus_ac,
            "sadra_ac_bus_i" => t_bus_ac,
            "sadra_loss_pu"  => get(conv, "sadra_loss_pu", false),
        )
    end

    # ------------------------------------------------------------------
    # 6. Add dummy generators at each DC bus
    # ------------------------------------------------------------------
    # Each dummy gen:
    #   - Sits at the DC bus
    #   - pg is free (bounded by ±Pacmax) — will be constrained to VSC losses
    #   - qg is free (bounded by ±Qacmax) — reactive slack for DC bus balance
    #   - cost = 0 (losses are physical, not a cost decision)

    for (_, conv) in data["convdc"]
        c_i   = conv["index"]
        new_i = dummy_gen_offset + c_i
        dc_bus_ac = dc_to_ac_bus[conv["busdc_i"]]

        # Losses are always consumed (never generated): pmax=0
        # ploss_max at rated current, using raw p.u. coefficients (matching PMACDC)
        if get(conv, "sadra_loss_pu", false)
            # FUBM path: converter current is uncapped (AIMMS-faithful), so
            # losses can be large. Use AIMMS dummy-gen bounds.
            pmax =  10.0; pmin = -10.0; qmax = 2.0; qmin = -2.0
        else
            # PMACDC path (v0.1): current is bounded by Imax, size pmin to the
            # worst-case loss at rated current, as before.
            imax_gen = get(conv, "Imax", 1.1)
            I_base   = baseMVA / (sqrt(3) * conv["basekVac"])
            lossA_pu = conv["LossA"] / baseMVA
            lossB_pu = conv["LossB"] * I_base / baseMVA
            lossC_pu = conv["LossCinv"] * I_base^2 / baseMVA
            ploss_max = lossA_pu + lossB_pu * imax_gen + lossC_pu * imax_gen^2
            pmax = 0.0
            pmin = -ploss_max
            qmax = get(conv, "Qacmax", 100.0) / baseMVA
            qmin = get(conv, "Qacmin", -100.0) / baseMVA
        end

        data["gen"]["$new_i"] = Dict{String,Any}(
            "index"      => new_i,
            "gen_bus"    => dc_bus_ac,
            "pg"         => 0.0,
            "qg"         => 0.0,
            "pmax"       => pmax,
            "pmin"       => pmin,
            "qmax"       => qmax,
            "qmin"       => qmin,
            "vg"         => 1.0,
            "mbase"      => baseMVA,
            "gen_status" => 1,
            "cost"       => [0.0, 0.0],
            "ncost"      => 2,
            "model"      => 2,
            "startup"    => 0.0,
            "shutdown"   => 0.0,
            # SADRA metadata
            "sadra_dummy_gen" => true,
            "sadra_conv_i"    => c_i,
            "sadra_vsc_branch_i" => conv_to_branch[c_i],  # link to VSC branch
        )
    end

    # ------------------------------------------------------------------
    # 7. Apply AC voltage setpoints for type_ac=2 converters
    #    Pin AC bus vm via tight bounds (eq: Vac control)
    # ------------------------------------------------------------------
    for (_, conv) in data["convdc"]
        if get(conv, "type_ac", 1) == 2
            ac_bus = conv["busac_i"]
            vtar = get(conv, "Vtar", 1.0)
            data["bus"]["$ac_bus"]["vmax"] = vtar
            data["bus"]["$ac_bus"]["vmin"] = vtar
        end
    end

    # ------------------------------------------------------------------
    # 8. Store the index maps in data["sadra"] for use by other files
    # ------------------------------------------------------------------
    data["sadra"] = Dict{String,Any}(
        "dc_to_ac_bus"    => dc_to_ac_bus,    # dc bus number -> AC bus index
        "conv_to_branch"  => conv_to_branch,   # conv index -> VSC branch index
        "conv_to_gen"     => conv_to_gen,      # conv index -> dummy gen index
        "ac_bus_offset"   => ac_bus_offset,
        "dc_branch_offset"=> dc_branch_offset,
        "vsc_branch_offset"=> vsc_branch_offset,
        "dummy_gen_offset" => dummy_gen_offset,
        "DC_X_STUB"       => DC_X_STUB,
    )

    return data
end


# ------------------------------------------------------------------
# Validation helper — call after sadra_transform! to sanity-check
# ------------------------------------------------------------------
function sadra_check(data::Dict{String,Any})
    s = data["sadra"]
    println("=== SADRA data check ===")
    println("DC buses added: ", length(data["busdc"]))
    println("DC branches added: ", length(data["branchdc"]))
    println("VSC branches added: ", length(data["convdc"]))
    println("Dummy generators added: ", length(data["convdc"]))
    println("\nDC bus -> AC bus mapping:")
    for (dc, ac) in sort(collect(s["dc_to_ac_bus"]))
        bus = data["bus"]["$ac"]
        println("  DC bus $dc -> AC bus $ac  (vm bounds: [$(bus["vmin"]), $(bus["vmax"])])")
    end
    println("\nVSC branches:")
    for (c_i, br_i) in sort(collect(s["conv_to_branch"]))
        br = data["branch"]["$br_i"]
        gen_i = s["conv_to_gen"][c_i]
        println("  Conv $c_i -> branch $br_i  (f=$(br["f_bus"]) DC -> t=$(br["t_bus"]) AC)")
        println("    r=$(br["br_r"]), x=$(br["br_x"]), ma∈[$(br["sadra_ma_min"]),$(br["sadra_ma_max"])]")
        println("    LossA=$(br["sadra_LossA"]), LossB=$(br["sadra_LossB"]), LossC=$(br["sadra_LossC"])")
        println("    dummy gen index: $gen_i at bus $(data["gen"]["$gen_i"]["gen_bus"])")
    end
    println("\nDC branches (AC representation):")
    for (_, brdc) in data["branchdc"]
        br_i = s["dc_branch_offset"] + brdc["index"]
        br   = data["branch"]["$br_i"]
        println("  DC branch $(brdc["index"]) -> AC branch $br_i  ($(br["f_bus"]) -> $(br["t_bus"]), r=$(br["br_r"]))")
    end
    println("========================")
end
