# fubm_ingest.jl
# Read a MATPOWER-FUBM .m case (DC data encoded inline in extended branch
# columns) and produce a data dict in the SAME shape that a PowerModelsACDC
# parse would give, so that sadra_transform!(data) runs downstream unchanged.
#
# This is needed because the SADRA paper's 1354-bus PEGASE case is distributed
# in FUBM format, not PMACDC format. FUBM stores DC buses as ordinary AC buses
# with index >= a high offset (10001..), DC lines as branches with x=0 between
# those buses, and VSCs as branches with CONV_A != 0.
#
# IMPORTANT differences from PMACDC input, handled here:
#   * Loss coefficients (ALPHA1/2/3) are in p.u. applied to p.u. current — NOT
#     MW/kA. We mark them sadra_loss_pu=true so sadra_transform! skips the kA
#     normalisation.
#   * Control assignment and all setpoints are read FROM THE FILE columns
#     (this .m has been updated to the AIMMS reference values — see the header
#     comment in the .m). Control TYPE comes from CONV_A; setpoints from
#     VT_SET / PF / KDP; loss coeffs from ALPHA1/2/3; DC bus bounds from the bus
#     matrix. The only value with no FUBM column is the droop voltage reference
#     Vref, which stays as a constant below.
#
# FUBM branch columns (1-indexed, per the case file header):
#   1 fbus 2 tbus 3 r 4 x 5 b 6 rateA 7 rateB 8 rateC 9 ratio/ma 10 angle
#   11 status 12 angmin 13 angmax 14 PF 15 QF 16 PT 17 QT 18 MU_SF 19 MU_ST
#   20 MU_ANGMIN 21 MU_ANGMAX 22 VF_SET 23 VT_SET 24 MA_MAX 25 MA_MIN
#   26 CONV_A 27 BEQ 28 K2 29 BEQ_MIN 30 BEQ_MAX 31 SH_MIN 32 SH_MAX
#   33 GSW 34 ALPHA1 35 ALPHA2 36 ALPHA3 37 KDP

# DC bus index threshold: FUBM buses at/above this are DC buses.
const FUBM_DC_BUS_MIN = 10000

# CONV_A (col 26) encodes the converter control type (FUBM convention):
#   1 = type I   : controls Pf (if a PF setpoint is given) else theta=0
#   2 = type II  : controls DC voltage (Vdc)
#   4 = type III : droop controlled
# The AC side controls Vac to VT_SET (col 23) when VT_SET != 0, else free.
# The type-I split (theta=0 vs Pf) is decided by whether PF (col 14) is nonzero:
# a type-I converter with a PF setpoint controls that power; one with PF=0 and
# only a voltage target holds its phase shift at zero.
#
# Droop voltage reference reads from VF_SET (col 22), the from-side (DC) voltage
# setpoint — semantically the droop's DC voltage reference. No constant needed.
#
# Controlled transformers are detected structurally (no hardcoded list): a
# branch that is NOT a VSC (CONV_A=0), IS a transformer (tap!=1 or shift!=0,
# which excludes DC lines), AND carries a control setpoint (PF!=0 => PST, or
# VT_SET!=0 => CTT). See the branch-classification loop below.


"""
    parse_fubm_matpower(path) -> Dict

Lightweight parser for the mpc.bus / mpc.gen / mpc.branch / mpc.gencost
matrices of a FUBM-format .m file. Handles both tab- and comma-delimited
rows and strips inline `%` comments. Returns raw numeric matrices plus baseMVA.
"""
function parse_fubm_matpower(path::AbstractString)
    blocks = Dict("bus"=>Vector{Vector{Float64}}(),
                  "gen"=>Vector{Vector{Float64}}(),
                  "branch"=>Vector{Vector{Float64}}(),
                  "gencost"=>Vector{Vector{Float64}}())
    baseMVA = 100.0
    cur = ""
    for raw in eachline(path)
        line = strip(raw)
        if startswith(line, "mpc.baseMVA")
            m = match(r"=\s*([0-9.eE+-]+)", line)
            m !== nothing && (baseMVA = parse(Float64, m.captures[1]))
            continue
        end
        matched = false
        for b in ("gencost","bus","gen","branch")
            if startswith(line, "mpc.$b")
                cur = b
                matched = true
                break
            end
        end
        matched && @goto nextline
        if cur != ""
            if startswith(line, "]")
                cur = ""
                @goto nextline
            end
            (isempty(line) || startswith(line, "%")) && @goto nextline
            body = split(line, "%")[1]
            body = replace(body, ";"=>"")
            parts = filter(!isempty, split(body, r"[\s,]+"))
            isempty(parts) && @goto nextline
            vals = Float64[]
            ok = true
            for p in parts
                v = tryparse(Float64, p)
                v === nothing ? (ok=false; break) : push!(vals, v)
            end
            ok && push!(blocks[cur], vals)
        end
        @label nextline
    end
    return Dict("baseMVA"=>baseMVA,
                "bus"=>blocks["bus"], "gen"=>blocks["gen"],
                "branch"=>blocks["branch"], "gencost"=>blocks["gencost"])
end


"""
    parse_fubm_acdc(path) -> data::Dict{String,Any}

Parse a FUBM .m case and build a PMACDC-shaped dict (bus/gen/branch/load +
busdc/branchdc/convdc) that sadra_transform! can consume directly.
Control type comes from the CONV_A column; setpoints, loss coeffs and DC bus
bounds are read from their respective columns (the .m has been updated to the
AIMMS reference values).
"""
function parse_fubm_acdc(path::AbstractString)
    raw = parse_fubm_matpower(path)
    baseMVA = raw["baseMVA"]

    data = Dict{String,Any}(
        "baseMVA"   => baseMVA,
        "per_unit"  => true,
        "name"      => "fubm_case",
        "bus"       => Dict{String,Any}(),
        "gen"       => Dict{String,Any}(),
        "branch"    => Dict{String,Any}(),
        "load"      => Dict{String,Any}(),
        "shunt"     => Dict{String,Any}(),
        "storage"   => Dict{String,Any}(),
        "switch"    => Dict{String,Any}(),
        "dcline"    => Dict{String,Any}(),
        "busdc"     => Dict{String,Any}(),
        "branchdc"  => Dict{String,Any}(),
        "convdc"    => Dict{String,Any}(),
    )

    # ---- AC buses (skip DC buses; those go in busdc) -------------------
    load_i = 0
    for r in raw["bus"]
        bi = Int(r[1])
        if bi >= FUBM_DC_BUS_MIN
            continue  # handled as DC bus below
        end
        data["bus"]["$bi"] = Dict{String,Any}(
            "index"=>bi, "bus_i"=>bi, "bus_type"=>Int(r[2]),
            "vm"=>r[8], "va"=>deg2rad(r[9]),
            "base_kv"=>r[10], "vmax"=>r[12], "vmin"=>r[13],
            "area"=>Int(r[7]), "zone"=>Int(r[11]),
        )
        # loads
        if r[3] != 0.0 || r[4] != 0.0
            load_i += 1
            data["load"]["$load_i"] = Dict{String,Any}(
                "index"=>load_i, "load_bus"=>bi,
                "pd"=>r[3]/baseMVA, "qd"=>r[4]/baseMVA, "status"=>1)
        end
        # shunts
        if r[5] != 0.0 || r[6] != 0.0
            si = length(data["shunt"]) + 1
            data["shunt"]["$si"] = Dict{String,Any}(
                "index"=>si, "shunt_bus"=>bi,
                "gs"=>r[5]/baseMVA, "bs"=>r[6]/baseMVA, "status"=>1)
        end
    end

    # ---- DC buses (FUBM bus index >= threshold) ------------------------
    # busdc_i is the LOCAL dc index (1..n). Map fubm bus number -> local index.
    fubm_to_dc = Dict{Int,Int}()
    dc_count = 0
    for r in raw["bus"]
        bi = Int(r[1])
        bi < FUBM_DC_BUS_MIN && continue
        dc_count += 1
        fubm_to_dc[bi] = dc_count
        data["busdc"]["$dc_count"] = Dict{String,Any}(
            "index"=>dc_count, "busdc_i"=>dc_count,
            "Vdc"=>r[8],
            "Vdcmax"=>r[12],   # bus matrix VmMax (now AIMMS-correct in the .m)
            "Vdcmin"=>r[13],   # bus matrix VmMin
            "basekVdc"=>r[10], "Pdc"=>0.0,
            "fubm_bus"=>bi)
    end

    # ---- Generators ----------------------------------------------------
    # gen cols: 1 bus 2 Pg 3 Qg 4 Qmax 5 Qmin 6 Vg 7 mBase 8 status
    #           9 Pmax 10 Pmin ...
    # NOTE: some gens (e.g. the slack) carry +/-Inf reactive limits in the .m.
    # JuMP cannot set a bound to +/-Inf, so we clamp to a large finite value.
    # PowerModels' own parser skips infinite bounds; clamping is behaviourally
    # equivalent for the optimizer and keeps the variable builder happy.
    BIG = 1e6 / baseMVA   # large p.u. bound
    clampbnd(x) = isinf(x) ? sign(x)*BIG*baseMVA : x
    for (gi, r) in enumerate(raw["gen"])
        gbus = Int(r[1])
        data["gen"]["$gi"] = Dict{String,Any}(
            "index"=>gi, "gen_bus"=>gbus,
            "pg"=>r[2]/baseMVA, "qg"=>r[3]/baseMVA,
            "qmax"=>clampbnd(r[4])/baseMVA, "qmin"=>clampbnd(r[5])/baseMVA,
            "vg"=>r[6], "mbase"=>r[7], "gen_status"=>Int(r[8]),
            "pmax"=>clampbnd(r[9])/baseMVA, "pmin"=>clampbnd(r[10])/baseMVA,
            "model"=>2, "ncost"=>2, "cost"=>[0.0,0.0],
            "startup"=>0.0, "shutdown"=>0.0)
    end
    # ---- Generator costs (gencost rows align with gen rows) ------------
    # gencost cols: 1 model 2 startup 3 shutdown 4 ncost 5.. coeffs
    # MATPOWER cost coefficients are defined for power in MW, but PowerModels
    # evaluates the objective against power in p.u. For a polynomial cost
    # sum_k c_k * P^(n-k), each coefficient on P^p must be multiplied by
    # baseMVA^p to keep the cost value invariant under P_MW = P_pu*baseMVA.
    # (_PM.parse_file does this automatically; our hand-rolled ingest must too.)
    for (gi, r) in enumerate(raw["gencost"])
        haskey(data["gen"], "$gi") || continue
        model = Int(r[1]); ncost = Int(r[4])
        coeffs = collect(r[5:5+ncost-1])
        # coeffs[1] is the highest-order term (P^(ncost-1)); scale by baseMVA^power
        scaled = similar(coeffs)
        for (k, c) in enumerate(coeffs)
            power = (ncost - 1) - (k - 1)   # exponent of P for this coefficient
            scaled[k] = c * baseMVA^power
        end
        data["gen"]["$gi"]["model"]   = model
        data["gen"]["$gi"]["ncost"]   = ncost
        data["gen"]["$gi"]["cost"]    = scaled
        data["gen"]["$gi"]["startup"] = r[2]
        data["gen"]["$gi"]["shutdown"]= r[3]
    end

    # ---- Branches: classify into AC line / DC line / VSC ---------------
    # column indices (1-based)
    C_CONV_A = 26; C_ALPHA1 = 34; C_ALPHA2 = 35; C_ALPHA3 = 36; C_KDP = 37
    C_VTSET = 23; C_VFSET = 22; C_PF = 14
    ac_bi = 0; dc_bi = 0; conv_i = 0
    for r in raw["branch"]
        fb = Int(r[1]); tb = Int(r[2])
        conv_a = length(r) >= C_CONV_A ? r[C_CONV_A] : 0.0
        is_dc_dc = (fb >= FUBM_DC_BUS_MIN) && (tb >= FUBM_DC_BUS_MIN)
        is_vsc   = conv_a != 0.0

        if is_vsc
            # VSC branch: from = DC bus (FUBM>=thr), to = AC bus
            conv_i += 1
            @assert fb >= FUBM_DC_BUS_MIN "VSC from-bus $fb not a DC bus"
            dcloc = fubm_to_dc[fb]

            vtset = r[C_VTSET]    # AC voltage target (0 => AC side free)
            pf    = r[C_PF]       # active power column (MW)
            vfset = r[C_VFSET]    # from-side (DC) voltage setpoint; droop Vref
            kdp   = length(r) >= C_KDP ? r[C_KDP] : 0.0

            # --- Control TYPE from CONV_A, with the documented type-I split ---
            #   type_dc codes used downstream: 0=theta0, 1=Pf, 2=Vdc, 3=droop
            if conv_a == 4            # type III: droop
                type_dc = 3
            elseif conv_a == 2        # type II: DC voltage
                type_dc = 2
            elseif conv_a == 1        # type I: Pf if a setpoint is given, else theta=0
                type_dc = (pf != 0.0) ? 1 : 0
            else
                error("Unrecognised CONV_A=$conv_a on VSC $fb->$tb")
            end
            # AC side: Vac control if VT_SET nonzero, else free.
            type_ac = (vtset != 0.0) ? 2 : 0

            data["convdc"]["$conv_i"] = Dict{String,Any}(
                "index"=>conv_i, "busdc_i"=>dcloc, "busac_i"=>tb,
                "rtf"=>r[3], "xtf"=>r[4], "bf"=>r[5], "filter"=> r[5]!=0 ? 1 : 0,
                "Vmmax"=>r[24], "Vmmin"=>r[25],
                "Imax"=>1.1,
                "Pacmax"=> r[6]!=0 ? r[6] : 1e4, "Pacmin"=>-(r[6]!=0 ? r[6] : 1e4),
                "Qacmax"=>1e4, "Qacmin"=>-1e4,
                "basekVac"=> get(data["bus"]["$tb"], "base_kv", 220.0),
                # loss coefficients from ALPHA1/2/3 columns (AIMMS values, p.u.
                # on p.u. current). transform reads sadra_loss_pu to skip kA norm.
                "LossA"=> r[C_ALPHA1],
                "LossB"=> r[C_ALPHA2],
                "LossCinv"=> r[C_ALPHA3],
                "sadra_loss_pu"=> true,
                # control fields, all from the file columns
                "type_dc"=>type_dc, "type_ac"=>type_ac,
                "P_g"=> pf,                 # MW; transform divides by baseMVA -> pu
                "Q_g"=> 0.0,
                "Vtar"=> vtset != 0.0 ? vtset : 1.0,
                "Vdcset"=> data["busdc"]["$dcloc"]["Vdc"],  # set below for type II
                # droop: slope from KDP col, Pref from PF col; Vref is the constant
                "droop"=> type_dc==3 ? kdp : 0.0,
                "Pdcset"=> type_dc==3 ? pf : 0.0,          # MW; /baseMVA in transform
                "droop_vref"=> vfset != 0.0 ? vfset : 1.0,
                "fubm_conv_a"=> conv_a,
            )
            # For type II (Vdc control), the DC voltage setpoint is the DC bus's
            # target. AIMMS pins it via the DC bus bound (VmMin=VmMax in the .m),
            # which we already read; use that bound as the explicit Vdcset too.
            if type_dc == 2
                data["convdc"]["$conv_i"]["Vdcset"] = data["busdc"]["$dcloc"]["Vdcmin"]
            end
        elseif is_dc_dc
            # DC line: between two DC buses, x=0 expected
            dc_bi += 1
            data["branchdc"]["$dc_bi"] = Dict{String,Any}(
                "index"=>dc_bi,
                "fbusdc"=>fubm_to_dc[fb], "tbusdc"=>fubm_to_dc[tb],
                "r"=>r[3], "rateA"=>r[6]!=0 ? r[6]/baseMVA : 1e4/baseMVA,
                "status"=>Int(r[11]))
        else
            # ordinary AC branch (may be a controlled PST/CTT — flagged below)
            ac_bi += 1
            br = Dict{String,Any}(
                "index"=>ac_bi, "f_bus"=>fb, "t_bus"=>tb,
                "br_r"=>r[3], "br_x"=>r[4], "br_b"=>r[5],
                "g_fr"=>0.0, "g_to"=>0.0, "b_fr"=>r[5]/2, "b_to"=>r[5]/2,
                "tap"=> r[9]==0 ? 1.0 : r[9], "shift"=>deg2rad(r[10]),
                "br_status"=>Int(r[11]),
                "angmin"=>deg2rad(r[12]), "angmax"=>deg2rad(r[13]),
                "rate_a"=> r[6]/baseMVA, "rate_b"=> (r[7]!=0 ? r[7] : r[6])/baseMVA,
                "rate_c"=> (r[8]!=0 ? r[8] : r[6])/baseMVA,
                "transformer"=> (r[9]!=0 && r[9]!=1) || r[10]!=0)
            # Controlled transformer? Detected structurally (no hardcoded list):
            #   - it is a transformer: tap != 1 or shift != 0 (DC lines, with
            #     tap=1 & shift=0, are already handled above and excluded here);
            #   - and it carries a control setpoint:
            #       PF (col 14) != 0   -> PST, controls from-side active power;
            #       VT_SET (col 23)!=0 -> CTT, controls to-side voltage.
            pf_x = r[C_PF]; vt_x = r[C_VTSET]
            is_xfmr = (r[9] != 0 && r[9] != 1) || r[10] != 0
            if is_xfmr && pf_x != 0.0
                br["sadra_xfmr_ctrl"] = "pst"
                br["sadra_pst_pf"] = pf_x / baseMVA   # MW -> pu (from-side P)
            elseif is_xfmr && vt_x != 0.0
                br["sadra_xfmr_ctrl"] = "ctt"
                br["sadra_ctt_vt"] = vt_x             # pu voltage target
            end
            data["branch"]["$ac_bi"] = br
        end
    end

    return data
end
