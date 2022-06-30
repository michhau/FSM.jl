function surf_props(ebm::EBM, cn::Constants, Sf, irow::Int = 1, icol::Int = 1)

    ebm.albs[ebm.am .== 0] = ebm.asmn .+ (ebm.asmx - ebm.asmn) .* (ebm.Tsurf[ebm.am .== 0] .- cn.Tm) ./ ebm.Talb
    cond_prog_alb = ebm.am .== 1  
    tau = similar(ebm.albs)
    fill!(tau, 3600.0 * ebm.tcld)
    tau[ebm.Tsurf .>= cn.Tm] .= 3600 * ebm.tmlt
    tau = tau[cond_prog_alb]
    rt = 1 ./ tau .+ Sf[cond_prog_alb] ./ ebm.Salb
    alim = (ebm.asmn ./ tau .+ Sf[cond_prog_alb] .* ebm.asmx ./ ebm.Salb) ./ rt
    ebm.albs[cond_prog_alb] = alim .+ (ebm.albs[cond_prog_alb] .- alim) .* exp.(-1 .* rt .* ebm.dt)

    # Snow albedo
    ebm.albs[ebm.albs .< min(ebm.asmx, ebm.asmn)] .= min(ebm.asmx, ebm.asmn)
    ebm.albs[ebm.albs .> max(ebm.asmx, ebm.asmn)] .= max(ebm.asmx, ebm.asmn)

    # Density of fresh snow
    ebm.rfs[ebm.dm .== 0] .= ebm.rho0
    ebm.rfs[ebm.dm .== 1] .= ebm.rhof

    # Thermal conductivity of snow
    ebm.ksnow .= ebm.kfix
    #ebm.ksnow[ebm.cm .== 0] .= ebm.kfix
    cond_dens_func = ebm.cm .== 1
    loc_cond_dens_func = findall(x->x==1, cond_dens_func)
    
    
    for idx in loc_cond_dens_func
        for k in 1:ebm.Nsnow[idx]
            rhos = ebm.rfs[idx]
            if (ebm.dm[idx] == 1 && ebm.Ds[idx, k] > eps(Float64))
                rhos = (ebm.Sice[idx, k] + ebm.Sliq[idx, k]) / ebm.Ds[idx, k]
            end
            ebm.ksnow[idx, k] = cn.hcon_ice * (rhos / cn.rho_ice)^ebm.bthr
        end
    end

    # Partial snow cover
    snowdepth = dropdims(sum(ebm.Ds, dims = 3), dims=3)
    ebm.fsnow = tanh.(snowdepth ./ ebm.hfsn)
    ebm.alb = ebm.fsnow .* ebm.albs .+ (1 .- ebm.fsnow) .* ebm.alb0
    ebm.z0 = (ebm.z0sn.^ebm.fsnow) .* (ebm.z0sf.^(1 .- ebm.fsnow))
    
    # Soil
    dPsidT = -cn.rho_ice * cn.Lf / (cn.rho_wat * cn.g * cn.Tm)

    for k in 1:ebm.Nsoil
        ebm.csoil[:,:, k] .= ebm.hcap_soil * ebm.Dzsoil[k]
        ebm.ksoil[:,:, k] .= ebm.hcon_soil

        cond_soil = ebm.theta[:,:,k] .> eps(Float64)
        dthudT = zeros(Float64, irow, icol)
        sthu = similar(dthudT)
        sthu[cond_soil] = reshape(ebm.theta[cond_soil, k], irow, icol)
        sthf = zeros(Float64, irow, icol)
        Tc = reshape(ebm.Tsoil[cond_soil, k] .- cn.Tm, irow, icol)
        Tmax = reshape(cn.Tm .+ (ebm.sathh / dPsidT) .* (ebm.Vsat ./ ebm.theta[cond_soil, k]).^ebm.b, irow, icol)

        cond_soil_2 = cond_soil .&& reshape(ebm.Tsoil[cond_soil,k], irow, icol) .< Tmax
        dthudT[cond_soil_2] = (-dPsidT * ebm.Vsat / (ebm.b * ebm.sathh)) .* (dPsidT .* Tc[cond_soil_2] ./ ebm.sathh).^(-1 / ebm.b - 1)
        sthu[cond_soil_2] = ebm.Vsat .* (dPsidT .* Tc[cond_soil_2] ./ ebm.sathh).^(-1 ./ ebm.b)
        sthu[cond_soil_2] = min.(sthu[cond_soil_2], ebm.theta[cond_soil_2, k])
        sthf[cond_soil_2] = (ebm.theta[cond_soil_2, k] .- sthu[cond_soil_2]) .* cn.rho_wat / cn.rho_ice

        Mf = cn.rho_ice * ebm.Dzsoil[k] .* sthf
        Mu = cn.rho_wat * ebm.Dzsoil[k] .* sthu
        ebm.csoil[cond_soil, k] = vec(ebm.hcap_soil * ebm.Dzsoil[k] .+ cn.hcap_ice .* Mf .+ cn.hcap_wat .* Mu .+ cn.rho_wat .* ebm.Dzsoil[k] .* ((cn.hcap_wat - cn.hcap_ice) .* Tc .+ cn.Lf) .* reshape(dthudT[cond_soil],irow, icol))
        Smf = cn.rho_ice .* sthf / (cn.rho_wat * ebm.Vsat)
        Smu = sthu ./ ebm.Vsat
        
        thice = similar(Smf)
        fill!(thice, 0.0)
        thice[Smf .> 0] = ebm.Vsat .+ Smf[Smf .> 0] ./ (Smu[Smf .> 0] .+ Smf[Smf .> 0])

        thwat = similar(Smu)
        fill!(thwat, 0.0)
        thwat[Smu .> 0.0] = ebm.Vsat .* Smu[Smu .> 0.0] ./ (Smu[Smu .> 0.0] .+ Smf[Smu .> 0.0])
        hcon_sat = ebm.hcon_soil * (cn.hcon_wat.^thwat) .* (cn.hcon_ice.^thice) ./ (cn.hcon_air.^ebm.Vsat)
        ebm.ksoil[cond_soil, k] = vec((hcon_sat .- ebm.hcon_soil) .* (Smf .+ Smu) .+ ebm.hcon_soil)
        if (k == 1)
            ebm.gs[cond_soil] = ebm.gsat .* max.((Smu .* ebm.Vsat ./ ebm.Vcrit).^2, 1)
        end
    end


    # Surface layer
    ebm.Dz1 = max.(ebm.Dzsoil[1], ebm.Ds[:,:, 1])
    ebm.Ts1 = ebm.Tsoil[:,:,1] .+ (ebm.Tsnow[:,:,1] .- ebm.Tsoil[:,:,1]) .* ebm.Ds[:,:,1] ./ ebm.Dzsoil[1]
    ebm.ksurf = ebm.Dzsoil[1] ./ (2 .* ebm.Ds[:,:,1] ./ ebm.ksnow[:,:,1] .+ (ebm.Dzsoil[1] .- 2 .* ebm.Ds[:,:,1]) ./ ebm.ksoil[:,:,1])

    ebm.ksurf[ebm.Ds[:,:,1] .> 0.5 * ebm.Dzsoil[1]] = ebm.ksnow[ebm.Ds[:,:,1] .> 0.5 * ebm.Dzsoil[1],1]

    ebm.Ts1[ebm.Ds[:,:,1] .> ebm.Dzsoil[1]] = ebm.Tsnow[ebm.Ds[:,:,1] .> ebm.Dzsoil[1],1]

    return nothing

end
