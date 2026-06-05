# run_pmacdc.jl
# Run this in a SEPARATE Julia session (NO Pkg.activate - use global environment)
#
# Instructions:
#   1. Open a NEW Julia terminal (not your SADRA one)
#   2. Run: include("C:/Users/hznh83/Downloads/sadra/test/run_pmacdc.jl")
#   3. Results saved to: C:/Users/hznh83/Downloads/sadra/test/pmacdc_results.json

using PowerModels, PowerModelsACDC, Ipopt, JSON
using JuMP: optimizer_with_attributes
import Logging
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

test_dir = joinpath(pkgdir(PowerModelsACDC), "test", "data")
cases = sort(filter(f -> endswith(f, ".m"), readdir(test_dir, join=true)))

# Add 3120-bus case
sadra_case = "C:/Users/hznh83/Downloads/sadra/test/data/case3120sp_acdc.m"
if isfile(sadra_case)
    insert!(cases, 1, sadra_case)
end

results = Dict{String,Any}()

for case in cases
    name = basename(case)
    data = PowerModels.parse_file(case)
    if !haskey(data, "busdc") || isempty(data["busdc"])
        continue
    end
    print("$name ... ")
    try
        try
            PowerModelsACDC.process_additional_data!(data)
        catch e
            if !occursin("Memento", string(e))
                rethrow(e)
            end
        end
        r = PowerModelsACDC.run_acdcopf(
            data, PowerModels.ACPPowerModel, ipopt;
            setting = Dict("output" => Dict("branch_flows" => true),
                           "conv_losses_mp" => true))
        results[name] = Dict(
            "status"     => string(r["termination_status"]),
            "objective"  => r["objective"],
            "solve_time" => r["solve_time"])
        println("$(r["termination_status"])  obj=$(round(r["objective"],digits=2))  t=$(round(r["solve_time"],digits=3))s")
    catch e
        results[name] = Dict("status" => "ERROR", "error" => string(e))
        println("ERROR: $e")
    end
end

out = "C:/Users/hznh83/Downloads/sadra/test/pmacdc_results.json"
open(out, "w") do f; JSON.print(f, results, 2); end
println("\nSaved to $out")
