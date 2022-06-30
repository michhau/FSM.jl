function surf_ebal(ebm::EBM, cn::Constants, Ta, Qa, Ua, Ps, SW, LW, irow::Int = 1, icol::Int = 1)
    
    Qs = qsat(false, Ps, ebm.Tsurf, cn)
 
    psi = ebm.gs ./ (ebm.gs .+ ebm.CH .* Ua)
    psi[Qs .< Qa .|| ebm.Sice[:,:,1] .> 0] .= 1.0
    
    rho = Ps ./ (cn.Rgas .* Ta)
    rKH = rho .* ebm.CH .* Ua
    #rho = nothing

    Lh = fill(cn.Ls, irow, icol)
    Lh[ebm.Tsurf .> cn.Tm] .= cn.Lc

    # Surface energy balance without melt
    D = Lh .* Qs ./ (cn.Rwat .* ebm.Tsurf .^2)
    Esurf = psi .* rKH .* (Qs .- Qa)
    ebm.Gsurf = 2 .* ebm.ksurf .* (ebm.Tsurf .- ebm.Ts1) ./ ebm.Dz1

    ebm.Hsurf = cn.cp .* rKH .* (ebm.Tsurf .- Ta)
    ebm.Lesrf = Lh .* Esurf
    ebm.Melt .= 0.0
    ebm.Rnet = (1 .- ebm.alb) .* SW .+ LW .- cn.sb .* ebm.Tsurf .^4
    dTs = (ebm.Rnet .- ebm.Hsurf .- ebm.Lesrf .- ebm.Gsurf) ./ ((cn.cp .+ Lh .* psi .* D) .* rKH .+ 2 .* ebm.ksurf ./ ebm.Dz1 .+ 4.0 .* cn.sb .* ebm.Tsurf .^3)
    dE = psi .* rKH .* D .* dTs
    dG = 2 .* ebm.ksurf .* dTs ./ ebm.Dz1
    dH = cn.cp .* rKH .* dTs
    dR = -cn.sb .* ebm.Tsurf .^3 .* dTs
    
    # Surface melting
    #define gridpoints with surface melting
    srfmelt = ebm.Tsurf .+ dTs .> cn.Tm .&& ebm.Sice[:,:,1] .> 0

    ebm.Melt[srfmelt] = sum(ebm.Sice[srfmelt,:], dims=2) ./ ebm.dt
    dTs[srfmelt] .-= (cn.Lf .* ebm.Melt[srfmelt]) ./ ((cn.cp .+ cn.Ls .* psi[srfmelt] .* D[srfmelt]) .* rKH[srfmelt] .+ 2 .* ebm.ksurf[srfmelt] ./ ebm.Dz1[srfmelt] .+ 4 .* cn.sb .* ebm.Tsurf[srfmelt] .^3)
    dE[srfmelt] ./= psi[srfmelt]
    #dG[srfmelt] = 2 .* ebm.ksurf[srfmelt] .* dTs[srfmelt] ./ ebm.Dz1[srfmelt]
    #dH[srfmelt] = cn.cp .* rKH .* dTs

    cond_tmp = ebm.Tsurf .+ dTs .< cn.Tm
    cond2 = srfmelt .& cond_tmp

    Qs[cond2] .= qsat(false, Ps[cond2], cn.Tm, cn)
    Esurf[cond2] .= rKH[cond2] .* (Qs[cond2] - Qa[cond2])
    ebm.Gsurf[cond2] .*= (cn.Tm .- ebm.Ts1[cond2])./(ebm.Tsurf[cond2] .- ebm.Ts1[cond2])
    ebm.Hsurf[cond2] .*= (cn.Tm .- Ta[cond2])./(ebm.Tsurf[cond2] .- Ta[cond2])
    ebm.Lesrf[cond2] .= cn.Ls .* Esurf[cond2]
    ebm.Rnet[cond2] .+= cn.sb .* (ebm.Tsurf[cond2] .^4 .- cn.Tm .^4)
    ebm.Melt[cond2] .= (ebm.Rnet[cond2] .- ebm.Hsurf[cond2] .- ebm.Lesrf[cond2] .- ebm.Gsurf[cond2]) ./ cn.Lf
    ebm.Melt[cond2] .= max.(ebm.Melt[cond2], 0.0)
    dE[cond2] .= 0.0
    dG[cond2] .= 0.0
    dH[cond2] .= 0.0
    dR[cond2] .= 0.0
    dTs[cond2] .= cn.Tm .- ebm.Tsurf[cond2]

    # Update surface temperature and fluxes
    ebm.Tsurf +=  dTs
    Esurf += dE
    ebm.Gsurf += dG
    ebm.Hsurf += dH
    ebm.Rnet += dR
    ebm.Esnow .= 0.0
    #Esoil = zeros(Float64, irow, icol)

    cond3 = (ebm.Sice[:,:,1] .> 0 .|| ebm.Tsurf .< Tm)

    ebm.Esnow[cond3] .= Esurf[cond3]
    ebm.Lesrf[cond3] .= cn.Ls .* Esurf[cond3]
    #Esoil[.!(cond3)] .= Esurf[cond3]
    ebm.Lesrf[.!(cond3)] .= cn.Lc .* Esurf[.!(cond3)]

    return nothing

end
