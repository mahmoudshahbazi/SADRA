## verify_control.jl
## Run the SADRA control-mode OPF on the FUBM 1354-bus PEGASE case and
## check the results against the paper's Table IV (control setpoints) and
## Tables V/VI (solved VSC voltages/powers).
##
## Usage (from the sadra package dir, env activated):
##   include("test/verify_control.jl")
##   verify_control("test/data/fubm_case1354pegase_2MTDC_ctrls_pf_qt_dp.m")

using SADRA, PowerModels, Ipopt, JuMP, Printf

# Paper Table IV setpoints, keyed by AC bus. (free => no AC constraint.)
const TARGETS = Dict(
    2072 => (vsc="VSC1", dc="theta=0",   ac_vac=1.08,  pf=nothing),
    8195 => (vsc="VSC2", dc="droop",     ac_vac=1.075, pf=nothing),
    6246 => (vsc="VSC3", dc="Vdc=1.01",  ac_vac=1.06,  pf=nothing),
    7282 => (vsc="VSC4", dc="Pf=-5.0pu", ac_vac=1.07,  pf=-5.0),
    352  => (vsc="VSC5", dc="Pf=-4.5pu", ac_vac=nothing, pf=-4.5),
)

function verify_control(casefile::String; baseMVA=100.0,
        optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6))
    println("="^70)
    println("SADRA control-mode verification — 1354-bus PEGASE")
    println("="^70)

    res = solve_sadra_fubm(casefile, optimizer)
    st  = res["termination_status"]
    obj = res["objective"]
    @printf("\nTermination: %s\nObjective:   %.2f\n\n", st, obj)

    sol = res["solution"]

    # Re-parse to recover index maps (which AC/DC bus each VSC touches).
    data = parse_fubm_acdc(casefile)
    SADRA.sadra_transform!(data)
    s = data["sadra"]

    println("VSC results vs paper Table IV/V/VI:")
    println("-"^70)
    @printf("%-5s %-12s %10s %10s %10s %10s\n",
            "VSC", "control", "vm_ac", "tgt_vac", "p_fr(pu)", "tgt_pf")
    for (c_i, br_i) in sort(collect(s["conv_to_branch"]))
        br = data["branch"]["$br_i"]
        ac_bus = br["sadra_ac_bus_i"]
        dc_bus = br["sadra_dc_bus_i"]
        tgt = get(TARGETS, ac_bus, nothing)

        vm_ac = haskey(sol["bus"], "$ac_bus") ? sol["bus"]["$ac_bus"]["vm"] : NaN
        vm_dc = haskey(sol["bus"], "$dc_bus") ? sol["bus"]["$dc_bus"]["vm"] : NaN
        # branch from-side active power
        p_fr = NaN
        if haskey(sol["branch"], "$br_i")
            p_fr = get(sol["branch"]["$br_i"], "pf", NaN)
        end

        name = tgt===nothing ? "?" : tgt.vsc
        ctrl = tgt===nothing ? "?" : tgt.dc
        tvac = (tgt===nothing || tgt.ac_vac===nothing) ? NaN : tgt.ac_vac
        tpf  = (tgt===nothing || tgt.pf===nothing)     ? NaN : tgt.pf
        @printf("%-5s %-12s %10.4f %10s %10.4f %10s\n",
                name, ctrl, vm_ac,
                isnan(tvac) ? "free" : @sprintf("%.4f",tvac),
                p_fr,
                isnan(tpf) ? "-" : @sprintf("%.4f",tpf))
    end

    println("\nDC bus voltages (should share one phase angle):")
    for (dc_i, ac_idx) in sort(collect(s["dc_to_ac_bus"]))
        if haskey(sol["bus"], "$ac_idx")
            b = sol["bus"]["$ac_idx"]
            @printf("  DC bus %d (idx %d): vm=%.4f  va=%.6f rad\n",
                    dc_i, ac_idx, b["vm"], b["va"])
        end
    end

    println("\nControlled transformers (Table IV: PST Pf=-1.73, CTT Vt=1.1):")
    for (i, br) in data["branch"]
        haskey(br, "sadra_xfmr_ctrl") || continue
        kind = br["sadra_xfmr_ctrl"]
        bi = parse(Int, i)
        if kind == "pst"
            p_fr = haskey(sol["branch"], i) ? get(sol["branch"][i], "pf", NaN) : NaN
            @printf("  PST branch %d: p_fr=%.4f  (target %.4f)\n", bi, p_fr, br["sadra_pst_pf"])
        else
            tb = br["t_bus"]
            vm = haskey(sol["bus"], "$tb") ? sol["bus"]["$tb"]["vm"] : NaN
            @printf("  CTT branch %d: vm_to=%.4f  (target %.4f)\n", bi, vm, br["sadra_ctt_vt"])
        end
    end

    # Reactive-free DC check: every DC line should carry q ~ 0.
    println("\nDC line reactive power (should be ~0):")
    maxq = 0.0
    for (_, brdc) in data["branchdc"]
        bi = s["dc_branch_offset"] + brdc["index"]
        if haskey(sol["branch"], "$bi")
            qf = get(sol["branch"]["$bi"], "qf", 0.0)
            qt = get(sol["branch"]["$bi"], "qt", 0.0)
            maxq = max(maxq, abs(qf), abs(qt))
        end
    end
    @printf("  max |q| on any DC line: %.2e pu\n", maxq)

    println("="^70)
    return res
end
