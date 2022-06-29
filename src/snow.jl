function snow(ebm::EBM, cn::Constants, Sf, Rf, Ta, D, S, W, row::Int = 1, col::Int = 1)

    ebm.Gsoil[row, col] = ebm.Gsurf[row, col]
    Roff = Rf * ebm.dt

    csnow = zeros(Float64, ebm.Nsmax)
    E = zeros(Float64, ebm.Nsmax)
    U = zeros(Float64, ebm.Nsmax)
    Gs = zeros(Float64, ebm.Nsmax)
    dTs = zeros(Float64, ebm.Nsmax)
    a = zeros(Float64, ebm.Nsmax)
    b = zeros(Float64, ebm.Nsmax)
    c = zeros(Float64, ebm.Nsmax)
    rhs = zeros(Float64, ebm.Nsmax)
    snow_gamma = zeros(Float64, ebm.Nsmax)

    if (ebm.Nsnow[row, col] > 0)   # Existing snowpack
        # Heat capacity
        for k = 1:ebm.Nsnow[row, col]
            csnow[k] = ebm.Sice[row, col, k] * cn.hcap_ice + ebm.Sliq[row, col, k] * cn.hcap_wat
        end

        # Heat conduction
        if (ebm.Nsnow[row, col] == 1)
            Gs[1] = 2.0 / (ebm.Ds[row, col, 1] / ebm.ksnow[row, col, 1] + ebm.Dzsoil[1] / ebm.ksoil[row, col, 1])
            dTs[1] = (ebm.Gsurf[row, col] + Gs[1] * (ebm.Tsoil[row, col, 1] - ebm.Tsnow[row, col, 1])) * ebm.dt / (csnow[1] + Gs[1] * ebm.dt)
        else
            for k = 1:(ebm.Nsnow[row, col] - 1)
                Gs[k] = 2.0 / (ebm.Ds[row, col, k] / ebm.ksnow[row, col, k] + ebm.Ds[row, col, k + 1] / ebm.ksnow[row, col, k + 1])
            end
            a[1] = 0.0
            b[1] = csnow[1] + Gs[1] * ebm.dt
            c[1] = - Gs[1] * ebm.dt
            rhs[1] = (ebm.Gsurf[row, col] - Gs[1] * (ebm.Tsnow[row, col, 1] - ebm.Tsnow[row, col, 2])) * ebm.dt
            for k = 2:(ebm.Nsnow[row, col] - 1)
                a[k] = c[k - 1]
                b[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * ebm.dt
                c[k] = - Gs[k] * ebm.dt
                rhs[k] = Gs[k - 1] * (ebm.Tsnow[row, col, k - 1] - ebm.Tsnow[row, col, k]) * ebm.dt + Gs[k] * (ebm.Tsnow[row, col, k + 1] - ebm.Tsnow[row, col, k]) * ebm.dt
            end
            k = ebm.Nsnow[row, col]

            Gs[k] = 2.0 / (ebm.Ds[row, col, k] / ebm.ksnow[row, col, k] + ebm.Dzsoil[1] / ebm.ksoil[row, col, 1])
            a[k] = c[k - 1]
            b[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * ebm.dt
            c[k] = 0.0
            rhs[k] = Gs[k - 1] * (ebm.Tsnow[row, col, k - 1] - ebm.Tsnow[row, col, k]) * ebm.dt + Gs[k] * (ebm.Tsoil[row, col, 1] - ebm.Tsnow[row, col, k]) * ebm.dt
            tridiag(ebm.Nsnow[row, col], snow_gamma, a, b, c, rhs, dTs)
        end


        for k = 1:ebm.Nsnow[row, col]
            ebm.Tsnow[row, col, k] = ebm.Tsnow[row, col, k] + dTs[k]
        end
        ebm.Gsoil[row, col] = Gs[ebm.Nsnow[row, col]] * (ebm.Tsnow[row, col, ebm.Nsnow[row, col]] - ebm.Tsoil[row, col, 1])



        # Convert melting ice to liquid water
        dSice = ebm.Melt[row, col] * ebm.dt
        for k = 1:ebm.Nsnow[row, col]
            coldcont = csnow[k] * (cn.Tm - ebm.Tsnow[row, col, k])
            if (coldcont < 0.0)
                # global dSice  ### HACK
                dSice = dSice - coldcont / cn.Lf
                ebm.Tsnow[row, col, k] = cn.Tm
            end
            if (dSice > 0.0)
                if (dSice > ebm.Sice[row, col, k])       # Layer melts completely
                    dSice = dSice - ebm.Sice[row, col, k]
                    ebm.Ds[row, col, k] = 0.0
                    ebm.Sliq[row, col, k] = ebm.Sliq[row, col, k] + ebm.Sice[row, col, k]
                    ebm.Sice[row, col, k] = 0.0
                else                       # Layer melts partially
                    ebm.Ds[row, col, k] = (1.0 - dSice / ebm.Sice[row, col, k]) * ebm.Ds[row, col, k]
                    ebm.Sice[row, col, k] = ebm.Sice[row, col, k] - dSice
                    ebm.Sliq[row, col, k] = ebm.Sliq[row, col, k] + dSice
                    dSice = 0.0                # Melt exhausted
                end
            end
        end


        # Remove snow by sublimation
        dSice = max(ebm.Esnow[row, col], 0.0) * ebm.dt
        if (dSice > 0.0)
            for k = 1:ebm.Nsnow[row, col]
                if (dSice > ebm.Sice[row, col, k])       # Layer sublimates completely
                    # global dSice #### hack
                    dSice = dSice - ebm.Sice[row, col, k]
                    ebm.Ds[row, col, k] = 0.0
                    ebm.Sice[row, col, k] = 0.0
                else                       # Layer sublimates partially
                    ebm.Ds[row, col, k] = (1.0 - dSice / ebm.Sice[row, col, k]) * ebm.Ds[row, col, k]
                    ebm.Sice[row, col, k] = ebm.Sice[row, col, k] - dSice
                    dSice = 0.0                # Sublimation exhausted
                end
            end
        end

        # Snow hydraulics
        if ebm.hm[row, col] == 0  #  Free-draining snow
            for k = 1:ebm.Nsnow[row, col]
                # global Roff ### hack
                Roff = Roff + ebm.Sliq[row, col, k]
                ebm.Sliq[row, col, k] = 0.0
            end
        elseif ebm.hm[row, col] == 1  #  Bucket storage
            for k = 1:ebm.Nsnow[row, col]
                # global Roff ### hack
                phi = 0.0
                if (ebm.Ds[row, col, k] > eps(Float64))
                    phi = 1.0 - ebm.Sice[row, col, k] / (cn.rho_ice * ebm.Ds[row, col, k])
                end
                SliqMax = cn.rho_wat * ebm.Ds[row, col, k] * phi * ebm.Wirr
                ebm.Sliq[row, col, k] = ebm.Sliq[row, col, k] + Roff
                Roff = 0.0
                if (ebm.Sliq[row, col, k] > SliqMax)       # Liquid capacity exceeded
                    Roff = ebm.Sliq[row, col, k] - SliqMax   # so drainage to next layer
                    ebm.Sliq[row, col, k] = SliqMax
                end
                coldcont = csnow[k] * (cn.Tm - ebm.Tsnow[row, col, k])
                if (coldcont > 0.0)            # Liquid can freeze
                    # global dSice #### hack
                    dSice = min(ebm.Sliq[row, col, k], coldcont / cn.Lf)
                    ebm.Sliq[row, col, k] = ebm.Sliq[row, col, k] - dSice
                    ebm.Sice[row, col, k] = ebm.Sice[row, col, k] + dSice
                    ebm.Tsnow[row, col, k] = ebm.Tsnow[row, col, k] + cn.Lf * dSice / csnow[k]
                end
            end
        end


        # Snow compaction
        if ebm.dm[row, col] == 0  # Fixed snow density
            for k = 1:ebm.Nsnow[row, col]
                ebm.Ds[row, col, k] = (ebm.Sice[row, col, k] + ebm.Sliq[row, col, k]) / ebm.rho0
            end
        elseif ebm.dm[row, col] == 1  # Prognostic snow density
            tau = 3600.0 * ebm.trho
            for k = 1:ebm.Nsnow[row, col]
                if (ebm.Ds[row, col, k] > eps(Float64))
                    rhos = (ebm.Sice[row, col, k] + ebm.Sliq[row, col, k]) / ebm.Ds[row, col, k]
                    if (ebm.Tsnow[row, col, k] >= cn.Tm)
                        if (rhos < ebm.rmlt)
                            rhos = ebm.rmlt + (rhos - ebm.rmlt) * exp(-ebm.dt / tau)
                        end
                    else
                        if (rhos < ebm.rcld)
                            rhos = ebm.rcld + (rhos - ebm.rcld) * exp(-ebm.dt / tau)
                        end
                    end
                    ebm.Ds[row, col, k] = (ebm.Sice[row, col, k] + ebm.Sliq[row, col, k]) / rhos
                end
            end
        end
    end  # Existing snowpack

    # Add snowfall and frost to layer 1
    dSice = Sf * ebm.dt - min(ebm.Esnow[row, col], 0.0) * ebm.dt
    ebm.Ds[row, col, 1] = ebm.Ds[row, col, 1] + dSice / ebm.rfs[row, col]
    ebm.Sice[row, col, 1] = ebm.Sice[row, col, 1] + dSice

    # New snowpack
    if (ebm.Nsnow[row, col] == 0 && ebm.Sice[row, col, 1] > 0.0)
        ebm.Nsnow[row, col] = 1
        ebm.Tsnow[row, col, 1] = min(Ta, cn.Tm)
    end

    # Calculate snow depth and SWE
    snowdepth = 0.0
    SWE = 0.0
    for k = 1:ebm.Nsnow[row, col]
       # global SWE
       # global snowdepth
        snowdepth = snowdepth + ebm.Ds[row, col, k]
        SWE = SWE + ebm.Sice[row, col, k] + ebm.Sliq[row, col, k]
    end

    # Store state of old layers
    D .= ebm.Ds[row, col, :]
    S .= ebm.Sice[row, col, :]
    W .= ebm.Sliq[row, col, :]

    for k = 1:ebm.Nsnow[row, col]
        csnow[k] = ebm.Sice[row, col, k] * cn.hcap_ice + ebm.Sliq[row, col, k] * cn.hcap_wat
        E[k] = csnow[k] * (ebm.Tsnow[row, col, k] - cn.Tm)
    end
    Nold = ebm.Nsnow[row, col]

    # Initialise new layers
    ebm.Ds[row, col, :] .= 0.0
    ebm.Sice[row, col, :] .= 0.0
    ebm.Sliq[row, col, :] .= 0.0
    ebm.Tsnow[row, col, :] .= cn.Tm
    U[:] .= 0.0
    ebm.Nsnow[row, col] = 0.0

    if (SWE > 0.0)       # Existing or new snowpack

        # Re-assign and count snow layers
        dnew = snowdepth
        ebm.Ds[row, col, 1] = dnew
        k = 1
        if (ebm.Ds[row, col, 1] > ebm.Dzsnow[1])
            for k = 1:ebm.Nsmax
                # global dnew   ##########  hack

                ebm.Ds[row, col, k] = ebm.Dzsnow[k]
                dnew = dnew - ebm.Dzsnow[k]
                if (dnew <= ebm.Dzsnow[k] || k == ebm.Nsmax)
                    ebm.Ds[row, col, k] = ebm.Ds[row, col, k] + dnew
                    break   ##### hack: was previously exit
                end
            end
        end
        # Nsnow = k   hackhackhack
        # ebm.Nsnow[row, col] = sum(ebm.Ds .> 0)

        ebm.Nsnow[row, col] = 0
        for Ds in ebm.Ds[row, col, :]
            if Ds > 0.0
                ebm.Nsnow[row, col] += 1
            end
        end


        # Fill new layers from the top downwards
        knew = 1
        dnew = ebm.Ds[row, col, 1]
        for kold = 1:Nold
            while true
                if (D[kold] < dnew)

                    # global dnew   ###### hack
                    # global knew   ###### hack


                    # Transfer all snow from old layer and move to next old layer
                    ebm.Sice[row, col, knew] = ebm.Sice[row, col, knew] + S[kold]
                    ebm.Sliq[row, col, knew] = ebm.Sliq[row, col, knew] + W[kold]
                    U[knew] = U[knew] + E[kold]
                    dnew = dnew - D[kold]
                    break   ##### hack   exit
                else
                    # Transfer some snow from old layer and move to next new layer
                    wt = dnew / D[kold]
                    ebm.Sice[row, col, knew] = ebm.Sice[row, col, knew] + wt * S[kold]
                    ebm.Sliq[row, col, knew] = ebm.Sliq[row, col, knew] + wt * W[kold]
                    U[knew] = U[knew] + wt * E[kold]
                    D[kold] = (1 - wt) * D[kold]
                    E[kold] = (1 - wt) * E[kold]
                    S[kold] = (1 - wt) * S[kold]
                    W[kold] = (1 - wt) * W[kold]
                    knew = knew + 1
                    if (knew > ebm.Nsnow[row, col])
                        break #### hack  exit
                    end
                    dnew = ebm.Ds[row, col, knew]
                end

            end
        end

        # Diagnose snow layer temperatures
        for k = 1:ebm.Nsnow[row, col]
            csnow[k] = ebm.Sice[row, col, k] * cn.hcap_ice + ebm.Sliq[row, col, k] * cn.hcap_wat
            if (csnow[k] > eps(Float64))
                ebm.Tsnow[row, col, k] = cn.Tm + U[k] / csnow[k]
            end
        end

    end

    return snowdepth, SWE
end
