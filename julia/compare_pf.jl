# compare_pf.jl
# Validate solve_sadra_pf against PowerModelsACDC's own AC/DC power flow
# on the same case file. Run from the sadra_control_pf project folder:
#   julia> include("compare_pf.jl")
#
# One-time setup in this environment:
#   julia> using Pkg; Pkg.activate("."); Pkg.add("PowerModelsACDC")

using Pkg; Pkg.activate(".")
using SADRA, PowerModels, PowerModelsACDC, Ipopt, JuMP

const _PMACDC = PowerModelsACDC

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
    "tol" => 1e-6, "print_level" => 0)

file = joinpath(@__DIR__, "test", "data", "case5_acdc.m")   # adjust path as needed . Other test case: case5_acdc_xfmr_only

# ------------------------------------------------------------------
# 1. PMACDC reference power flow
# ------------------------------------------------------------------
data_ref = PowerModels.parse_file(file)
_PMACDC.process_additional_data!(data_ref)
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)

# function name depends on PMACDC version: run_acdcpf (<=0.7) / solve_acdcpf
solve_ref = isdefined(_PMACDC, :solve_acdcpf) ? _PMACDC.solve_acdcpf :
                                                _PMACDC.run_acdcpf
res_ref = solve_ref(data_ref, ACPPowerModel, ipopt; setting = s)
println("PMACDC PF: ", res_ref["termination_status"])

# ------------------------------------------------------------------
# 2. SADRA power flow
# ------------------------------------------------------------------
res_s = solve_sadra_pf(file, ipopt)
println("SADRA  PF: ", res_s["termination_status"])

# ------------------------------------------------------------------
# 3. Compare AC bus voltages (magnitude and angle)
# ------------------------------------------------------------------
println("\nAC bus comparison (vm / va_deg):")
maxdv = 0.0; maxda = 0.0
ac_ids = sort(parse.(Int, collect(keys(res_ref["solution"]["bus"]))))
for i in ac_ids
    bref = res_ref["solution"]["bus"]["$i"]
    bs   = res_s["solution"]["bus"]["$i"]
    dv = bs["vm"] - bref["vm"]
    da = rad2deg(bs["va"] - bref["va"])
    global maxdv = max(maxdv, abs(dv)); global maxda = max(maxda, abs(da))
    println(rpad("bus $i", 9),
        "PMACDC ", rpad(round(bref["vm"], digits=6), 10),
        "SADRA ",  rpad(round(bs["vm"],  digits=6), 10),
        "dvm=", rpad(round(dv, sigdigits=3), 11),
        "dva=", round(da, sigdigits=3), " deg")
end
println("max |dvm| = $maxdv,  max |dva| = $maxda deg")

# ------------------------------------------------------------------
# 4. DC bus voltages
# PMACDC reports them under "busdc". SADRA's DC buses became AC buses
# with new indices; print both sides and match using the transform's
# mapping (sadra_check on the transformed data prints dc bus -> new id).
# ------------------------------------------------------------------
println("\nPMACDC DC bus voltages:")
for i in sort(parse.(Int, collect(keys(res_ref["solution"]["busdc"]))))
    println("  busdc $i  vm = ", res_ref["solution"]["busdc"]["$i"]["vm"])
end
println("\nSADRA DC bus voltages (DC 1,2,3 -> AC 6,7,8):")
for b in ["6","7","8"]
    println("  bus $b  vm = ", res_s["solution"]["bus"][b]["vm"])
end
println("\nSADRA transformed-case DC bus mapping (run once to see indices):")
data_s = PowerModels.parse_file(file)
sadra_transform!(data_s)
sadra_check(data_s)
