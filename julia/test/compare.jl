# compare.jl
# Compares SADRA results against PMACDC results loaded from pmacdc_results.json.
# Run run_pmacdc.jl first in a separate Julia session to generate the JSON file.

using SADRA, PowerModels, PowerModelsACDC, Ipopt, Printf, JSON
using JuMP: optimizer_with_attributes
import Logging
try
    using Memento
    setlevel!(getlogger(PowerModels), "error")
catch end
Logging.disable_logging(Logging.Warn)

ipopt = optimizer_with_attributes(Ipopt.Optimizer,
    "max_iter" => 500,
    "tol" => 1e-6,
    "acceptable_tol" => 1e-4,
    "acceptable_iter" => 5,
    "acceptable_dual_inf_tol" => 100.0,
    "acceptable_constr_viol_tol" => 1e-4,
    "nlp_scaling_method" => "gradient-based",
    "print_level" => 0)

# Load PMACDC results from JSON
pmacdc_json = joinpath(@__DIR__, "pmacdc_results.json")
if !isfile(pmacdc_json)
    error("pmacdc_results.json not found. Run test/run_pmacdc.jl first in a separate Julia session.")
end
pmacdc_results = JSON.parsefile(pmacdc_json)
println("Loaded PMACDC results for $(length(pmacdc_results)) cases.")

function is_acdc_case(data::Dict)
    return haskey(data, "busdc") && length(data["busdc"]) > 0 &&
           haskey(data, "convdc") && length(data["convdc"]) > 0 &&
           haskey(data, "branchdc") && length(data["branchdc"]) > 0
end

function run_sadra(case_file, optimizer)
    try
        data = PowerModels.parse_file(case_file)
        if !is_acdc_case(data)
            return nothing, "not AC/DC case"
        end
        result = solve_sadra_opf(case_file, optimizer)
        return result, nothing
    catch e
        return nothing, string(e)
    end
end

function print_row(label, sadra_val, pmacdc_val)
    @printf("  %-12s  %-22s  %-22s\n", label,
        sadra_val === nothing ? "—" : string(sadra_val),
        pmacdc_val === nothing ? "—" : string(pmacdc_val))
end

function compare_case(case_file, optimizer, pmacdc_results)
    name = basename(case_file)
    println("\n══════════════════════════════════════════════════════")
    println("  Case: $name")
    println("══════════════════════════════════════════════════════")

    r_sadra, err_sadra = run_sadra(case_file, optimizer)

    if err_sadra == "not AC/DC case"
        println("  Skipping — no DC components")
        return
    end

    pmacdc = get(pmacdc_results, name, nothing)

    @printf("  %-12s  %-22s  %-22s\n", "", "SADRA", "PMACDC")
    @printf("  %-12s  %-22s  %-22s\n", "─"^12, "─"^22, "─"^22)

    if err_sadra === nothing
        print_row("Status",
            string(r_sadra["termination_status"]),
            pmacdc === nothing ? "no result" : pmacdc["status"])
        print_row("Objective",
            @sprintf("%.4f", r_sadra["objective"]),
            pmacdc === nothing ? "—" : (pmacdc["status"] == "ERROR" ? "ERROR" : @sprintf("%.4f", pmacdc["objective"])))
        print_row("Solve time",
            @sprintf("%.3fs", r_sadra["solve_time"]),
            pmacdc === nothing ? "—" : (pmacdc["status"] == "ERROR" ? "ERROR" : @sprintf("%.3fs", pmacdc["solve_time"])))

        if pmacdc !== nothing && pmacdc["status"] != "ERROR" && haskey(pmacdc, "objective")
            gap = abs(r_sadra["objective"] - pmacdc["objective"]) /
                  max(abs(pmacdc["objective"]), 1e-10) * 100
            @printf("  %-12s  %.4f%%\n", "Obj gap", gap)
        end
    else
        println("  SADRA error:  $err_sadra")
        if pmacdc !== nothing
            println("  PMACDC:       $(pmacdc["status"])")
        end
    end
end

# Collect all cases to test
test_dir = joinpath(pkgdir(PowerModelsACDC), "test", "data")

cases = sort(filter(f -> endswith(f, ".m"), readdir(test_dir, join=true)))

# Add 3120-bus case from sadra test/data
sadra_case = joinpath(@__DIR__, "data", "case3120sp_acdc.m")
if isfile(sadra_case)
    insert!(cases, 1, sadra_case)
end

# Run comparison - save to file and print to screen
out_file = joinpath(@__DIR__, "comparison_results.txt")
open(out_file, "w") do io
    redirect_stdout(io) do
        for case in cases
            compare_case(case, ipopt, pmacdc_results)
        end
        println("\n\nDone.")
    end
end

# Also print to screen
for case in cases
    compare_case(case, ipopt, pmacdc_results)
end
println("\n\nDone.")
println("Results saved to $out_file")
