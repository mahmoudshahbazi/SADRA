# test/test_3120.jl
# Run SADRA on the 3120-bus AC/DC case and compare against PMACDC reference.
# Expected objective: ~2,142,635 (matches PMACDC NLP and SADRA AIMMS results)

using SADRA
using Ipopt
using PowerModels

# Path to case file — adjust as needed
CASE_FILE = joinpath(@__DIR__, "data", "case3120sp_acdc.m")

println("=== SADRA OPF: case3120sp_acdc ===")
println("Parsing and transforming data...")

# Quick data check before solving
data = PowerModels.parse_file(CASE_FILE)
sadra_transform!(data)
sadra_check(data)

println("\nSolving...")
result = solve_sadra_opf(
    CASE_FILE,
    Ipopt.Optimizer;
    setting = Dict("output" => Dict("branch_flows" => true))
)

println("\n=== Results ===")
println("Status:    ", result["termination_status"])
println("Objective: ", result["objective"])
println("Solve time: ", result["solve_time"], " s")

# Reference value from paper Table III
ref_obj = 2_142_635.0
gap = abs(result["objective"] - ref_obj) / ref_obj * 100
println("\nReference (SADRA AIMMS): $ref_obj")
println("Gap: $(round(gap, digits=4))%")

# Print DC bus voltages
println("\n=== DC bus voltages ===")
s = data["sadra"]
for (dc_i, ac_i) in sort(collect(s["dc_to_ac_bus"]))
    vm = get(get(result["solution"]["bus"], "$ac_i", Dict()), "vm", "?")
    println("DC bus $dc_i (AC bus $ac_i): vm = $vm")
end

# Print converter results
println("\n=== Converter results ===")
for (conv_i, br_i) in sort(collect(s["conv_to_branch"]))
    gen_i = s["conv_to_gen"][conv_i]
    pg = get(get(result["solution"]["gen"], "$gen_i", Dict()), "pg", "?")
    qg = get(get(result["solution"]["gen"], "$gen_i", Dict()), "qg", "?")
    println("Conv $conv_i: pg_dummy=$pg, qg_dummy=$qg")
end
