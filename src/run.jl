
"point"
function run!(ebm::EBM, cn::Constants, snowdepth::Vector{Float64}, SWE::Vector{Float64}, Tsurf::Vector{Float64}, in::Input)

    D = zeros(ebm.Nsmax)
    S = zeros(ebm.Nsmax)
    W = zeros(ebm.Nsmax)

    for i in 1:length(in.year)

        year = in.year[i]
        month = in.month[i]
        day = in.day[i]
        hour = in.hour[i]
        SW = in.SW[i]
        LW = in.LW[i]
        Sf = in.Sf[i]
        Rf = in.Rf[i]
        Ta = in.Ta[i]
        RH = in.RH[i]
        Ua = in.Ua[i]
        Ps = in.Ps[i]

        Qs = qsat(true, Ps, Ta, cn)
        Qa = (RH / 100) * Qs

        Ua = max(Ua, 0.1)

        surf_props(ebm, cn, Sf)

        for i in 1:6
            surf_exch(ebm, cn, Ta, Ua)
            surf_ebal(ebm, cn, Ta, Qa, Ua, Ps, SW, LW)
        end

        snowdepth[i], SWE[i] = snow(ebm, cn, Sf, Rf, Ta, D, S, W)

        soil(ebm)

        Tsurf[i] = ebm.Tsurf[row, col]

    end

end

"spatial"
function run!(ebm::EBM, cn::Constants, snowdepth::Matrix{Float64}, SWE::Matrix{Float64}, Tsurf::Matrix{Float64},
    SW::Matrix{Float64}, LW::Matrix{Float64}, Sf::Matrix{Float64}, Rf::Matrix{Float64}, Ta::Matrix{Float64}, RH::Matrix{Float64}, Ua::Matrix{Float64}, Ps::Matrix{Float64}, nrow::Int, ncol::Int)

    D = zeros(Float64, nrow, ncol, ebm.Nsmax)
    S = zeros(Float64, nrow, ncol, ebm.Nsmax)
    W = zeros(Float64, nrow, ncol, ebm.Nsmax)

    Ua = max.(Ua, 0.1)

    Qs = qsat(true, Ps, Ta, cn)
    Qa = (RH ./ 100) .* Qs

    #println("surf_props")
    #@time 
    surf_props(ebm, cn, Sf, nrow, ncol)

    #println("surf_exch and _ebal 6 times")
    kdx = 1
    while kdx < 7
        #@time 
        surf_exch(ebm, cn, Ta, Ua, nrow, ncol)
        #println("surf_exch done")
        #@time 
        surf_ebal(ebm, cn, Ta, Qa, Ua, Ps, SW, LW, nrow, ncol)
        #println("surf_ebal done")
        kdx = kdx + 1
    end

    #println("snow")
    #@time
    (snowdepth, SWE) = snow(ebm, cn, Sf, Rf, Ta, D, S, W, nrow, ncol)

    #println("soil")
    #@time
    soil(ebm, nrow, ncol)

    Tsurf = ebm.Tsurf

end
