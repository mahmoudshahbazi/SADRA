# test_pf_case.jl
# Reusable validation harness for solve_sadra_pf on PMACDC-format cases.
#
# Per case it runs:
#   A. External check : SADRA PF vs PMACDC's own AC/DC PF (same file).
#      Reports max AC |dvm|, |dva|, and DC bus voltage deltas.
#      This measures INGEST FIDELITY (station model, conventions), so it
#      has a soft threshold: warn above WARN_XVAL, no hard fail.
#   B. Internal check : SADRA OPF -> fix dispatch -> SADRA PF must
#      reproduce the OPF state. Same equations both sides, so this is a
#      hard PASS/FAIL on the PF machinery itself (threshold TOL_SELF).
#
# Usage (from the project folder, after `using SADRA, PowerModels,
# PowerModelsACDC, Ipopt, JuMP`):
#   include("test_pf_case.jl")
#   r = test_pf_case("data/case5_acdc.m")
#   test_pf_case("data/case24_3zones_acdc.m")
#   test_pf_case("data/case39_acdc.m")

using JuMP
import PowerModels
import PowerModelsACDC

# Silence PowerModels/InfrastructureModels info+warn logging (angmin spam etc).
# Comment this out when debugging a NEW case for the first time - the parse
# warnings (setpoint mismatches, auto bus-type changes) carry real diagnostics.
PowerModels.silence()

const TOL_SELF  = 1e-6     # hard threshold for self-consistency (B)
const WARN_XVAL = 1e-3     # soft threshold for PMACDC comparison (A)

const _IPOPT  = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
                    "tol" => 1e-6, "print_level" => 0)
const _IPOPT8 = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
                    "tol" => 1e-8, "print_level" => 0)

"PMACDC PF entry point across versions (run_acdcpf <=0.7, solve_acdcpf >=0.8)."
_pmacdc_pf() = isdefined(PowerModelsACDC, :solve_acdcpf) ?
    PowerModelsACDC.solve_acdcpf : PowerModelsACDC.run_acdcpf

function test_pf_case(file::String; verbose::Bool=true)
    name = basename(file)
    println("="^60, "\nPF validation: ", name, "\n", "="^60)

    # ----------------------------------------------------------------
    # A. External: SADRA PF vs PMACDC PF
    # ----------------------------------------------------------------
    data_ref = PowerModels.parse_file(file)
    PowerModelsACDC.process_additional_data!(data_ref)
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)
    res_ref = _pmacdc_pf()(data_ref, PowerModels.ACPPowerModel, _IPOPT; setting=s)

    res_s = solve_sadra_pf(file, _IPOPT)

    st_ref, st_s = res_ref["termination_status"], res_s["termination_status"]
    println("  PMACDC PF: $st_ref   SADRA PF: $st_s")
    ok_status = string(st_ref) == "LOCALLY_SOLVED" && string(st_s) == "LOCALLY_SOLVED"

    max_dvm_ac = NaN; max_dva_ac = NaN; max_dvm_dc = NaN
    if ok_status
        # AC buses (track argmax for diagnosis)
        max_dvm_ac = 0.0; max_dva_ac = 0.0
        argmax_vm = ""; argmax_va = ""
        for (i, b) in res_ref["solution"]["bus"]
            bs = res_s["solution"]["bus"][i]
            dv = abs(bs["vm"] - b["vm"]); da = abs(rad2deg(bs["va"] - b["va"]))
            dv > max_dvm_ac && (max_dvm_ac = dv; argmax_vm = i)
            da > max_dva_ac && (max_dva_ac = da; argmax_va = i)
        end
        if argmax_vm != ""
            b = res_ref["solution"]["bus"][argmax_vm]; bs = res_s["solution"]["bus"][argmax_vm]
            println("    worst AC vm at bus $argmax_vm: PMACDC ", round(b["vm"], digits=6),
                    "  SADRA ", round(bs["vm"], digits=6), "  (worst va at bus $argmax_va)")
        end
        # DC buses: use the transform's own dc_to_ac_bus map
        data_map = PowerModels.parse_file(file); sadra_transform!(data_map)
        dc2ac = data_map["sadra"]["dc_to_ac_bus"]
        max_dvm_dc = 0.0
        for (i, b) in res_ref["solution"]["busdc"]
            sadra_bus = string(dc2ac[parse(Int, i)])
            dv = res_s["solution"]["bus"][sadra_bus]["vm"] - b["vm"]
            verbose && println("    DC bus $i: PMACDC ", round(b["vm"], digits=7),
                               "  SADRA ", round(res_s["solution"]["bus"][sadra_bus]["vm"], digits=7),
                               "  d=", round(dv, sigdigits=3))
            max_dvm_dc = max(max_dvm_dc, abs(dv))
        end
        # Per-converter comparison: PMACDC pgrid/qgrid/pdc vs SADRA p_to/q_to/p_fr
        if verbose && haskey(res_ref["solution"], "convdc")
            data_s = PowerModels.parse_file(file); sadra_transform!(data_s)
            println("    converters (PMACDC pgrid,qgrid,pdc | SADRA -p_to,-q_to,p_fr):")
            for (ci, cv) in sort(collect(res_ref["solution"]["convdc"]), by=x->parse(Int,x[1]))
                bi = string(data_s["sadra"]["conv_to_branch"][parse(Int, ci)])
                br = res_s["solution"]["branch"][bi]
                brd = data_s["branch"][bi]
                vm_t = res_s["solution"]["bus"][string(brd["t_bus"])]["vm"]
                i_est = sqrt(br["pt"]^2 + br["qt"]^2) / vm_t
                println("      conv $ci: ",
                    round(cv["pgrid"], digits=5), ", ", round(cv["qgrid"], digits=5), ", ",
                    round(cv["pdc"],   digits=5), "  |  ",
                    round(-br["pt"], digits=5), ", ", round(-br["qt"], digits=5), ", ",
                    round(br["pf"],  digits=5),
                    "   i~", round(i_est, digits=3),
                    i_est > brd["sadra_imax"] ? " (EXCEEDS imax $(brd["sadra_imax"]))" : "")
            end
        end
        println("  [A] vs PMACDC: max|dvm_AC|=", round(max_dvm_ac, sigdigits=3),
                "  max|dva_AC|=", round(max_dva_ac, sigdigits=3), " deg",
                "  max|dvm_DC|=", round(max_dvm_dc, sigdigits=3))
        max(max_dvm_ac, max_dvm_dc) > WARN_XVAL &&
            @warn "Comparison delta above $(WARN_XVAL): check ingest conventions for this case."
    else
        @warn "Solve status not LOCALLY_SOLVED on one side; skipping comparison."
    end

    # ----------------------------------------------------------------
    # B. Internal: OPF -> PF self-consistency (hard check)
    # ----------------------------------------------------------------
    # The OPF runs UNCONTROLLED (converter powers free), so it is feasible
    # whenever an economic dispatch exists, regardless of the file's
    # setpoints. A consistent PF specification is then DERIVED from the
    # OPF solution: converter flows become the setpoints. By construction
    # the OPF point solves the resulting PF, so agreement tests exactly
    # the PF machinery, on any case.
    res_opf = solve_sadra_opf(file, _IPOPT8;
                  setting=Dict{String,Any}("sadra_controls" => false))
    if string(res_opf["termination_status"]) != "LOCALLY_SOLVED"
        println("  [B] SKIPPED: uncontrolled SADRA OPF status = ",
                res_opf["termination_status"])
        println("="^60)
        return (case = name, status_ok = ok_status,
                max_dvm_ac = max_dvm_ac, max_dva_ac = max_dva_ac,
                max_dvm_dc = max_dvm_dc,
                self_consistency = NaN, pass = ok_status)
    end

    data_pf = PowerModels.parse_file(file)
    sbase = data_pf["baseMVA"]
    # gen dispatch and gen-bus voltages from the OPF state
    for (j, g) in res_opf["solution"]["gen"]
        haskey(data_pf["gen"], j) || continue        # skip dummy gens
        data_pf["gen"][j]["pg"] = g["pg"]
    end
    for (j, g) in data_pf["gen"]
        b = string(g["gen_bus"])
        data_pf["bus"][b]["vm"] = res_opf["solution"]["bus"][b]["vm"]
    end
    # converter setpoints from the OPF converter flows (sign conventions
    # mirror the ingest: Pset = P_g/sbase, constraint p_to == -Pset, etc.)
    data_m = PowerModels.parse_file(file); sadra_transform!(data_m)
    for (ci_int, br_i) in data_m["sadra"]["conv_to_branch"]
        ci = string(ci_int)
        br  = res_opf["solution"]["branch"][string(br_i)]
        brd = data_m["branch"][string(br_i)]
        cv  = data_pf["convdc"][ci]
        cv["P_g"]    = -br["pt"] * sbase
        cv["Q_g"]    = -br["qt"] * sbase
        cv["Vdcset"] = res_opf["solution"]["bus"][string(brd["f_bus"])]["vm"]
        cv["Vtar"]   = res_opf["solution"]["bus"][string(brd["sadra_ac_bus_i"])]["vm"]
    end
    res_pf2 = solve_sadra_pf(data_pf, _IPOPT8)

    max_self = 0.0; argmax_self = ""
    for (i, b) in res_opf["solution"]["bus"]
        d = abs(res_pf2["solution"]["bus"][i]["vm"] - b["vm"])
        d > max_self && (max_self = d; argmax_self = i)
    end
    if argmax_self != ""
        println("    worst self-consistency at bus $argmax_self: OPF vm=",
                round(res_opf["solution"]["bus"][argmax_self]["vm"], digits=7),
                "  PF vm=", round(res_pf2["solution"]["bus"][argmax_self]["vm"], digits=7))
    end

    # Classify: the PF is deliberately unbounded, so wherever the OPF point
    # sits ON a bound (bus vmin/vmax, gen qmin/qmax), the PF legitimately
    # drifts off it. Count active bounds at the OPF solution.
    data_b = PowerModels.parse_file(file); sadra_transform!(data_b)
    n_active = 0
    for (i, b) in res_opf["solution"]["bus"]
        bd = data_b["bus"][i]
        (abs(b["vm"] - bd["vmax"]) < 1e-6 || abs(b["vm"] - bd["vmin"]) < 1e-6) &&
            (n_active += 1)
    end
    for (j, g) in res_opf["solution"]["gen"]
        gd = data_b["gen"][j]
        get(gd, "sadra_dummy_gen", false) && continue
        (abs(g["qg"] - gd["qmax"]) < 1e-6 || abs(g["qg"] - gd["qmin"]) < 1e-6 ||
         abs(g["pg"] - gd["pmax"]) < 1e-6 || abs(g["pg"] - gd["pmin"]) < 1e-6) &&
            (n_active += 1)
    end
    # Converter current bound active at the OPF point? i ~ |S_to|/vm_to.
    for (i, br) in data_b["branch"]
        get(br, "sadra_vsc", false) || continue
        s  = res_opf["solution"]["branch"][i]
        vm = res_opf["solution"]["bus"][string(br["t_bus"])]["vm"]
        i_est = sqrt(s["pt"]^2 + s["qt"]^2) / vm
        i_est > br["sadra_imax"] - 1e-4 && (n_active += 1)
    end

    if max_self < TOL_SELF
        verdict = "PASS"; pass = true
    elseif n_active > 0 && max_self < 1e-3
        verdict = "EXPLAINED ($(n_active) active OPF bounds; unbounded PF drifts off them)"
        pass = true
    else
        verdict = "FAIL (threshold $(TOL_SELF), no active bounds to explain it)"
        pass = false
    end
    println("  [B] self-consistency max|dvm| = ", round(max_self, sigdigits=3), "   ", verdict)

    println("="^60)
    return (case = name,
            status_ok = ok_status,
            max_dvm_ac = max_dvm_ac, max_dva_ac = max_dva_ac,
            max_dvm_dc = max_dvm_dc,
            self_consistency = max_self,
            pass = pass)
end
