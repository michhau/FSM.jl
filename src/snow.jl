function snow(ebm::EBM, cn::Constants, Sf, Rf, Ta, D, S, W)

    ebm.Gsoil = ebm.Gsurf
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

    if (ebm.Nsnow > 0)   # Existing snowpack
        # Heat capacity
        for k = 1:ebm.Nsnow
            csnow[k] = ebm.Sice[k] * cn.hcap_ice + ebm.Sliq[k] * cn.hcap_wat
        end

        # Heat conduction
        if (ebm.Nsnow == 1)
            Gs[1] = 2.0 / (ebm.Ds[1] / ebm.ksnow[1] + ebm.Dzsoil[1] / ebm.ksoil[1])
            dTs[1] = (ebm.Gsurf + Gs[1] * (ebm.Tsoil[1] - ebm.Tsnow[1])) * ebm.dt / (csnow[1] + Gs[1] * ebm.dt)
        else
            for k = 1:(ebm.Nsnow - 1)
                Gs[k] = 2.0 / (ebm.Ds[k] / ebm.ksnow[k] + ebm.Ds[k + 1] / ebm.ksnow[k + 1])
            end
            a[1] = 0.0
            b[1] = csnow[1] + Gs[1] * ebm.dt
            c[1] = - Gs[1] * ebm.dt
            rhs[1] = (ebm.Gsurf - Gs[1] * (ebm.Tsnow[1] - ebm.Tsnow[2])) * ebm.dt
            for k = 2:(ebm.Nsnow - 1)
                a[k] = c[k - 1]
                b[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * ebm.dt
                c[k] = - Gs[k] * ebm.dt
                rhs[k] = Gs[k - 1] * (ebm.Tsnow[k - 1] - ebm.Tsnow[k]) * ebm.dt + Gs[k] * (ebm.Tsnow[k + 1] - ebm.Tsnow[k]) * ebm.dt
            end
            k = ebm.Nsnow

            Gs[k] = 2.0 / (ebm.Ds[k] / ebm.ksnow[k] + ebm.Dzsoil[1] / ebm.ksoil[1])
            a[k] = c[k - 1]
            b[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * ebm.dt
            c[k] = 0.0
            rhs[k] = Gs[k - 1] * (ebm.Tsnow[k - 1] - ebm.Tsnow[k]) * ebm.dt + Gs[k] * (ebm.Tsoil[1] - ebm.Tsnow[k]) * ebm.dt
            tridiag(ebm.Nsnow, snow_gamma, a, b, c, rhs, dTs)
        end


        for k = 1:ebm.Nsnow
            ebm.Tsnow[k] = ebm.Tsnow[k] + dTs[k]
        end
        ebm.Gsoil = Gs[ebm.Nsnow] * (ebm.Tsnow[ebm.Nsnow] - ebm.Tsoil[1])



        # Convert melting ice to liquid water
        dSice = ebm.Melt * ebm.dt
        for k = 1:ebm.Nsnow
            coldcont = csnow[k] * (cn.Tm - ebm.Tsnow[k])
            if (coldcont < 0.0)
                # global dSice  ### HACK
                dSice = dSice - coldcont / cn.Lf
                ebm.Tsnow[k] = cn.Tm
            end
            if (dSice > 0.0)
                if (dSice > ebm.Sice[k])       # Layer melts completely
                    dSice = dSice - ebm.Sice[k]
                    ebm.Ds[k] = 0.0
                    ebm.Sliq[k] = ebm.Sliq[k] + ebm.Sice[k]
                    ebm.Sice[k] = 0.0
                else                       # Layer melts partially
                    ebm.Ds[k] = (1.0 - dSice / ebm.Sice[k]) * ebm.Ds[k]
                    ebm.Sice[k] = ebm.Sice[k] - dSice
                    ebm.Sliq[k] = ebm.Sliq[k] + dSice
                    dSice = 0.0                # Melt exhausted
                end
            end
        end


        # Remove snow by sublimation
        dSice = max(ebm.Esnow, 0.0) * ebm.dt
        if (dSice > 0.0)
            for k = 1:ebm.Nsnow
                if (dSice > ebm.Sice[k])       # Layer sublimates completely
                    # global dSice #### hack
                    dSice = dSice - ebm.Sice[k]
                    ebm.Ds[k] = 0.0
                    ebm.Sice[k] = 0.0
                else                       # Layer sublimates partially
                    ebm.Ds[k] = (1.0 - dSice / ebm.Sice[k]) * ebm.Ds[k]
                    ebm.Sice[k] = ebm.Sice[k] - dSice
                    dSice = 0.0                # Sublimation exhausted
                end
            end
        end

        # Snow hydraulics
        if ebm.hm == 0  #  Free-draining snow
            for k = 1:ebm.Nsnow
                # global Roff ### hack
                Roff = Roff + ebm.Sliq[k]
                ebm.Sliq[k] = 0.0
            end
        elseif ebm.hm == 1  #  Bucket storage
            for k = 1:ebm.Nsnow
                # global Roff ### hack
                phi = 0.0
                if (ebm.Ds[k] > eps(Float64))
                    phi = 1.0 - ebm.Sice[k] / (cn.rho_ice * ebm.Ds[k])
                end
                SliqMax = cn.rho_wat * ebm.Ds[k] * phi * ebm.Wirr
                ebm.Sliq[k] = ebm.Sliq[k] + Roff
                Roff = 0.0
                if (ebm.Sliq[k] > SliqMax)       # Liquid capacity exceeded
                    Roff = ebm.Sliq[k] - SliqMax   # so drainage to next layer
                    ebm.Sliq[k] = SliqMax
                end
                coldcont = csnow[k] * (cn.Tm - ebm.Tsnow[k])
                if (coldcont > 0.0)            # Liquid can freeze
                    # global dSice #### hack
                    dSice = min(ebm.Sliq[k], coldcont / cn.Lf)
                    ebm.Sliq[k] = ebm.Sliq[k] - dSice
                    ebm.Sice[k] = ebm.Sice[k] + dSice
                    ebm.Tsnow[k] = ebm.Tsnow[k] + cn.Lf * dSice / csnow[k]
                end
            end
        end


        # Snow compaction
        if ebm.dm == 0  # Fixed snow density
            for k = 1:ebm.Nsnow
                ebm.Ds[k] = (ebm.Sice[k] + ebm.Sliq[k]) / ebm.rho0
            end
        elseif ebm.dm == 1  # Prognostic snow density
            tau = 3600.0 * ebm.trho
            for k = 1:ebm.Nsnow
                if (ebm.Ds[k] > eps(Float64))
                    rhos = (ebm.Sice[k] + ebm.Sliq[k]) / ebm.Ds[k]
                    if (ebm.Tsnow[k] >= cn.Tm)
                        if (rhos < ebm.rmlt)
                            rhos = ebm.rmlt + (rhos - ebm.rmlt) * exp(-ebm.dt / tau)
                        end
                    else
                        if (rhos < ebm.rcld)
                            rhos = ebm.rcld + (rhos - ebm.rcld) * exp(-ebm.dt / tau)
                        end
                    end
                    ebm.Ds[k] = (ebm.Sice[k] + ebm.Sliq[k]) / rhos
                end
            end
        end
    end  # Existing snowpack

    # Add snowfall and frost to layer 1
    dSice = Sf * ebm.dt - min(ebm.Esnow, 0.0) * ebm.dt
    ebm.Ds[1] = ebm.Ds[1] + dSice / ebm.rfs
    ebm.Sice[1] = ebm.Sice[1] + dSice

    # New snowpack
    if (ebm.Nsnow == 0 && ebm.Sice[1] > 0.0)
        ebm.Nsnow = 1
        ebm.Tsnow[1] = min(Ta, cn.Tm)
    end

    # Calculate snow depth and SWE
    snowdepth = 0.0
    SWE = 0.0
    for k = 1:ebm.Nsnow
       # global SWE
       # global snowdepth
        snowdepth = snowdepth + ebm.Ds[k]
        SWE = SWE + ebm.Sice[k] + ebm.Sliq[k]
    end

    # Store state of old layers
    D .= ebm.Ds
    S .= ebm.Sice
    W .= ebm.Sliq

    for k = 1:ebm.Nsnow
        csnow[k] = ebm.Sice[k] * cn.hcap_ice + ebm.Sliq[k] * cn.hcap_wat
        E[k] = csnow[k] * (ebm.Tsnow[k] - cn.Tm)
    end
    Nold = ebm.Nsnow

    # Initialise new layers
    ebm.Ds[:] .= 0.0
    ebm.Sice[:] .= 0.0
    ebm.Sliq[:] .= 0.0
    ebm.Tsnow[:] .= cn.Tm
    U[:] .= 0.0
    ebm.Nsnow = 0.0

    if (SWE > 0.0)       # Existing or new snowpack

        # Re-assign and count snow layers
        dnew = snowdepth
        ebm.Ds[1] = dnew
        k = 1
        if (ebm.Ds[1] > ebm.Dzsnow[1])
            for k = 1:ebm.Nsmax
                # global dnew   ##########  hack

                ebm.Ds[k] = ebm.Dzsnow[k]
                dnew = dnew - ebm.Dzsnow[k]
                if (dnew <= ebm.Dzsnow[k] || k == ebm.Nsmax)
                    ebm.Ds[k] = ebm.Ds[k] + dnew
                    break   ##### hack: was previously exit
                end
            end
        end
        # Nsnow = k   hackhackhack
        # ebm.Nsnow = sum(ebm.Ds .> 0)

        ebm.Nsnow = 0
        for Ds in ebm.Ds
            if Ds > 0.0
                ebm.Nsnow += 1
            end
        end


        # Fill new layers from the top downwards
        knew = 1
        dnew = ebm.Ds[1]
        for kold = 1:Nold
            while true
                if (D[kold] < dnew)

                    # global dnew   ###### hack
                    # global knew   ###### hack


                    # Transfer all snow from old layer and move to next old layer
                    ebm.Sice[knew] = ebm.Sice[knew] + S[kold]
                    ebm.Sliq[knew] = ebm.Sliq[knew] + W[kold]
                    U[knew] = U[knew] + E[kold]
                    dnew = dnew - D[kold]
                    break   ##### hack   exit
                else
                    # Transfer some snow from old layer and move to next new layer
                    wt = dnew / D[kold]
                    ebm.Sice[knew] = ebm.Sice[knew] + wt * S[kold]
                    ebm.Sliq[knew] = ebm.Sliq[knew] + wt * W[kold]
                    U[knew] = U[knew] + wt * E[kold]
                    D[kold] = (1 - wt) * D[kold]
                    E[kold] = (1 - wt) * E[kold]
                    S[kold] = (1 - wt) * S[kold]
                    W[kold] = (1 - wt) * W[kold]
                    knew = knew + 1
                    if (knew > ebm.Nsnow)
                        break #### hack  exit
                    end
                    dnew = ebm.Ds[knew]
                end

            end
        end

        # Diagnose snow layer temperatures
        for k = 1:ebm.Nsnow
            csnow[k] = ebm.Sice[k] * cn.hcap_ice + ebm.Sliq[k] * cn.hcap_wat
            if (csnow[k] > eps(Float64))
                ebm.Tsnow[k] = cn.Tm + U[k] / csnow[k]
            end
        end

    end

    return snowdepth, SWE
end
