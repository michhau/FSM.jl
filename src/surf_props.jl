function surf_props(ebm::EBM, cn::Constants, Sf, row::Int = 1, col::Int = 1)

    # Snow albedo

    if ebm.am[row, col] == 0  # Diagnosed snow albedo
        ebm.albs[row, col] = ebm.asmn + (ebm.asmx - ebm.asmn) * (ebm.Tsurf[row, col] - cn.Tm) / ebm.Talb
    elseif ebm.am[row, col] == 1  # Prognostic snow albedo
        tau = 3600.0 * ebm.tcld
        if (ebm.Tsurf[row, col] >= cn.Tm)
            tau = 3600.0 * ebm.tmlt
        end
        rt = 1.0 / tau + Sf / ebm.Salb
        alim = (ebm.asmn / tau + Sf * ebm.asmx / ebm.Salb) / rt
        ebm.albs[row, col] = alim + (ebm.albs[row, col] - alim) * exp(-rt * ebm.dt)
    end

    if (ebm.albs[row, col] < min(ebm.asmx, ebm.asmn))
        ebm.albs[row, col] = min(ebm.asmx, ebm.asmn)
    end

    if (ebm.albs[row, col] > max(ebm.asmx, ebm.asmn))
        ebm.albs[row, col] = max(ebm.asmx, ebm.asmn)
    end

    # Density of fresh snow

    if (ebm.dm[row, col] == 0)
        ebm.rfs[row, col] = ebm.rho0
    elseif (ebm.dm[row, col] == 1)
        ebm.rfs[row, col] = ebm.rhof
    end

    # Thermal conductivity of snow

    ebm.ksnow[row, col, :] .= ebm.kfix
    if (ebm.cm[row, col] == 0) # Fixed
        ebm.ksnow[row, col, :] .= ebm.kfix
    elseif (ebm.cm[row, col] == 1)  # Density function
        for k in 1:ebm.Nsnow[row, col]
            rhos = ebm.rfs[row, col]
            if (ebm.dm[row, col] == 1 && ebm.Ds[row, col, k] > eps(Float64))
                rhos = (ebm.Sice[row, col, k] + ebm.Sliq[row, col, k]) / ebm.Ds[row, col, k]
            end
            ebm.ksnow[row, col, k] = cn.hcon_ice * (rhos / cn.rho_ice)^ebm.bthr
        end
    end

    # Partial snow cover

    snowdepth = sum(ebm.Ds[row, col, :])
    ebm.fsnow[row, col] = tanh(snowdepth / ebm.hfsn)
    ebm.alb[row, col] = ebm.fsnow[row, col] * ebm.albs[row, col] + (1.0 - ebm.fsnow[row, col]) * ebm.alb0
    ebm.z0[row, col] = (ebm.z0sn^ebm.fsnow[row, col]) * (ebm.z0sf^(1.0 - ebm.fsnow[row, col]))

    # Soil

    dPsidT = -cn.rho_ice * cn.Lf / (cn.rho_wat * cn.g * cn.Tm)

    for k in 1:ebm.Nsoil
        ebm.csoil[row, col, k] = ebm.hcap_soil * ebm.Dzsoil[k]
        ebm.ksoil[row, col, k] = ebm.hcon_soil
        if (ebm.theta[row, col, k] > eps(Float64))
            dthudT = 0.0
            sthu = ebm.theta[row, col, k]
            sthf = 0.0
            Tc = ebm.Tsoil[row, col, k] - cn.Tm
            Tmax = cn.Tm + (ebm.sathh / dPsidT) * (ebm.Vsat / ebm.theta[row, col, k])^ebm.b
            if (ebm.Tsoil[row, col, k] < Tmax)
                dthudT = (-dPsidT * ebm.Vsat / (ebm.b * ebm.sathh)) * (dPsidT * Tc / ebm.sathh)^(-1 / ebm.b - 1)
                sthu = ebm.Vsat * (dPsidT * Tc / ebm.sathh)^(-1.0 / ebm.b)
                sthu = min(sthu, ebm.theta[row, col, k])
                sthf = (ebm.theta[row, col, k] - sthu) * cn.rho_wat / cn.rho_ice
            end
            Mf = cn.rho_ice * ebm.Dzsoil[k] * sthf
            Mu = cn.rho_wat * ebm.Dzsoil[k] * sthu
            ebm.csoil[row, col, k] = ebm.hcap_soil * ebm.Dzsoil[k] + cn.hcap_ice * Mf + cn.hcap_wat * Mu + cn.rho_wat * ebm.Dzsoil[k] * ((cn.hcap_wat - cn.hcap_ice) * Tc + cn.Lf) * dthudT
            Smf = cn.rho_ice * sthf / (cn.rho_wat * ebm.Vsat)
            Smu = sthu / ebm.Vsat
            thice = 0.0
            if (Smf > 0.0)
                thice = ebm.Vsat * Smf / (Smu + Smf)
            end
            thwat = 0.0
            if (Smu > 0.0)
                thwat = ebm.Vsat * Smu / (Smu + Smf)
            end
            hcon_sat = ebm.hcon_soil * (cn.hcon_wat^thwat) * (cn.hcon_ice^thice) / (cn.hcon_air^ebm.Vsat)
            ebm.ksoil[row, col, k] = (hcon_sat - ebm.hcon_soil) * (Smf + Smu) + ebm.hcon_soil
            if (k == 1)
                ebm.gs[row, col] = ebm.gsat * max((Smu * ebm.Vsat / ebm.Vcrit)^2.0, 1.0)
            end
        end
    end

    # Surface layer

    ebm.Dz1[row, col] = max(ebm.Dzsoil[1], ebm.Ds[row, col, 1])

    ebm.Ts1[row, col] = ebm.Tsoil[row, col, 1] + (ebm.Tsnow[row, col, 1] - ebm.Tsoil[row, col, 1]) * ebm.Ds[row, col, 1] / ebm.Dzsoil[1]

    ebm.ksurf[row, col] = ebm.Dzsoil[1] / (2.0 * ebm.Ds[row, col, 1] / ebm.ksnow[row, col, 1] + (ebm.Dzsoil[1] - 2 * ebm.Ds[row, col, 1]) / ebm.ksoil[row, col, 1])

    if (ebm.Ds[row, col, 1] > 0.5 * ebm.Dzsoil[1])
        ebm.ksurf[row, col] = ebm.ksnow[row, col, 1]
    end

    if (ebm.Ds[row, col, 1] > ebm.Dzsoil[1])
        ebm.Ts1[row, col] = ebm.Tsnow[row, col, 1]
    end

    return nothing

end
