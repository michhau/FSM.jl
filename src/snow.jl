function snow(ebm::EBM, cn::Constants, Sf, Rf, Ta, D, S, W, irow::Int=1, icol::Int=1)

    snowdepth = zeros(Float64, irow, icol)
    SWE = zeros(Float64, irow, icol)

    ebm.Gsoil = ebm.Gsurf
    Roff = Rf .* ebm.dt

    coldcont = zeros(Float64, irow, icol, ebm.Nsmax)

    ################################################
    #condition existing snowpack
    cond_ex_sp = ebm.Nsnow .> 0

    nsnow_array = repeat(cond_ex_sp, 1, 1, ebm.Nsmax)
    for k in 1:ebm.Nsmax
        nsnow_array[:, :, k] .= nsnow_array[:, :, k] .&& k .== ebm.Nsnow
    end

    #Heat capacity
    for i in 1:ebm.Nsmax
        kind = cond_ex_sp .&& ebm.Nsnow .<= i
        ebm.snow_csnow[kind, i] = ebm.Sice[kind, i] .* cn.hcap_ice .+ ebm.Sliq[kind, i] .* cn.hcap_wat
    end
    #Heat conduction
    cond_heat_cond = cond_ex_sp .&& ebm.Nsnow .== 1
    ebm.snow_Gs[cond_heat_cond, 1] = 2 ./ (ebm.Ds[cond_heat_cond, 1] ./ ebm.ksnow[cond_heat_cond, 1] .+ ebm.Dzsoil[1] ./ ebm.ksoil[cond_heat_cond, 1])
    ebm.snow_dTs[cond_heat_cond, 1] = (ebm.Gsurf[cond_heat_cond] .+ ebm.snow_Gs[cond_heat_cond, 1] .* (ebm.Tsoil[cond_heat_cond, 1] .- ebm.Tsnow[cond_heat_cond, 1])) .* ebm.dt ./ (ebm.snow_csnow[1] .+ ebm.snow_Gs[cond_heat_cond, 1] .* ebm.dt)
    #else
    cond_heat_else = cond_ex_sp .&& ebm.Nsnow .> 1
    for i in 1:ebm.Nsmax-1
        kind = cond_heat_else .&& i .<= ebm.Nsnow .- 1
        ebm.snow_Gs[kind, i] = 2 ./ (ebm.Ds[kind, i] ./ ebm.ksnow[kind, i] .+ ebm.Ds[kind, i+1] ./ ebm.ksnow[kind, i+1])
    end
    ebm.snow_a[cond_heat_else, 1] .= 0
    ebm.snow_b[cond_heat_else, 1] = ebm.snow_csnow[cond_heat_else, 1] .+ ebm.snow_Gs[cond_heat_else, 1] .* ebm.dt
    ebm.snow_c[cond_heat_else, 1] = -1 .* ebm.snow_Gs[cond_heat_else, 1] .* ebm.dt
    ebm.snow_rhs[cond_heat_else, 1] = (ebm.Gsurf[cond_heat_else] .- ebm.snow_Gs[cond_heat_else, 1] .* (ebm.Tsnow[cond_heat_else, 1] .- ebm.Tsnow[cond_heat_else, 2])) .* ebm.dt
    for i in 2:ebm.Nsmax.-1
        kind = cond_heat_else .&& i .<= ebm.Nsnow .- 1
        ebm.snow_a[kind, i] = ebm.snow_c[kind, i-1]
        ebm.snow_b[kind, i] = ebm.snow_csnow[kind, i] .+ (ebm.snow_Gs[kind, i-1] .+ ebm.snow_Gs[kind, i]) .* ebm.dt
        ebm.snow_c[kind, i] = -1 .* ebm.snow_Gs[kind, i] .* ebm.dt
        ebm.snow_rhs[kind, i] = ebm.snow_Gs[kind, i-1] .* (ebm.Tsnow[kind, i-1] .- ebm.Tsnow[kind, i]) .* ebm.dt .+ ebm.snow_Gs[kind, i] .* (ebm.Tsnow[kind, i+1] .- ebm.Tsnow[kind, i]) .* ebm.dt
    end
    for i in 2:ebm.Nsmax
        kind = cond_heat_else .&& i .== ebm.Nsnow
        ebm.snow_Gs[kind, i] = 2 ./ (ebm.Ds[kind, i] ./ ebm.ksnow[kind, i] .+ ebm.Dzsoil[1] ./ ebm.ksoil[kind, 1])
        ebm.snow_a[kind, i] = ebm.snow_c[kind, i-1]
        ebm.snow_b[kind, i] = ebm.snow_csnow[kind, i] .+ (ebm.snow_Gs[kind, i-1] .+ ebm.snow_Gs[kind, i]) .* ebm.dt
        ebm.snow_c[kind, i] .= 0
        ebm.snow_rhs[kind, i] = ebm.snow_Gs[kind, i-1] .* (ebm.Tsnow[kind, i-1] .- ebm.Tsnow[kind, i]) .* ebm.dt .+ ebm.snow_Gs[kind, i] .* (ebm.Tsoil[kind, 1] .- ebm.Tsnow[kind, i]) .* ebm.dt
        tridiag(ebm.Nsnow, ebm.snow_gamma, ebm.snow_a, ebm.snow_b, ebm.snow_c, ebm.snow_rhs, ebm.snow_dTs)
    end
    #end else

    ebm.Tsnow .+= ebm.snow_dTs
    ebm.Gsoil[cond_ex_sp] = ebm.snow_Gs[nsnow_array] .* (ebm.Tsnow[nsnow_array] .- ebm.Tsoil[cond_ex_sp, 1])

    # Convert melting ice to liquid water
    dSice = ebm.Melt .* ebm.dt
    for k in 1:ebm.Nsmax
        cond_melt_ice = cond_ex_sp .&& ebm.Nsnow .<= k
        coldcont[cond_melt_ice, k] = ebm.snow_csnow[cond_melt_ice, k] .* (cn.Tm .- ebm.Tsnow[cond_melt_ice, k])
        cond_coldconts1 = cond_melt_ice .&& coldcont[:,:, k] .< 0
        dSice[cond_coldconts1] .-= coldcont[cond_coldconts1, k] ./ cn.Lf
        ebm.Tsnow[cond_coldconts1, k] .= cn.Tm

        cond_dSiceg0 = cond_ex_sp .&& dSice .> 0
        cond_dSicegebmSice = cond_dSiceg0 .&& dSice .> ebm.Sice[:,:,k] #layer melts completely
        ebm.Ds[cond_dSicegebmSice, k] .= 0
        ebm.Sliq[cond_dSicegebmSice, k] .+= ebm.Sice[cond_dSicegebmSice, k]
        ebm.Sice[cond_dSicegebmSice, k] .= 0
        #else layer melts partially
        cond_dSicegebmSiceelse = cond_dSiceg0 .&& .!(dSice .> ebm.Sice[:,:,k])
        ebm.Ds[cond_dSicegebmSiceelse, k] = (1 .- dSice[cond_dSicegebmSiceelse] ./ ebm.Sice[cond_dSicegebmSiceelse, k]) .* ebm.Ds[cond_dSicegebmSiceelse, k]
        ebm.Sice[cond_dSicegebmSiceelse, k] = ebm.Sice[cond_dSicegebmSiceelse, k] .- dSice[cond_dSicegebmSiceelse]
        ebm.Sliq[cond_dSicegebmSiceelse, k] = ebm.Sliq[cond_dSicegebmSiceelse, k] .+ dSice[cond_dSicegebmSiceelse]
        dSice[cond_dSicegebmSiceelse] .= 0.0                # Melt exhausted


    end

    #Remove snow by sublimation
    dSice = max.(ebm.Esnow, 0) .* ebm.dt
    cond_dSiceg0 = cond_ex_sp .&& dSice .> 0
    for k in 1:ebm.Nsmax
        kind = cond_dSiceg0 .&& k .<= ebm.Nsnow
        cond_comp_subl = kind .&& dSice .> ebm.Sice[:,:,k] #layer sublimates completely
        dSice[cond_comp_subl] .-= ebm.Sice[cond_comp_subl, k]
        ebm.Ds[cond_comp_subl, k] .= 0
        ebm.Sice[cond_comp_subl, k] .= 0
        #else partially sublimating layer
        cond_part_subl = kind .&& .!(dSice .> ebm.Sice[:,:,k])
        ebm.Ds[cond_part_subl, k] = (1 .- dSice[cond_part_subl] ./ ebm.Sice[cond_part_subl, k]) .* ebm.Ds[cond_part_subl, k]
        ebm.Sice[cond_part_subl, k] = ebm.Sice[cond_part_subl, k] .- dSice[cond_part_subl]
        dSice[cond_part_subl] .= 0 #sublimation exhausted
    end

    #Snow hydraulics
    cond_free_draining = cond_ex_sp .&& ebm.hm .== 0 #Free-draining snow
    for k in 1:ebm.Nsmax
        kind = cond_free_draining .&& k .<= ebm.Nsnow
        Roff[kind] += ebm.Sliq[kind, k]
        ebm.Sliq[kind, k] .= 0
    end
    cond_bucket_storage = cond_ex_sp .&& ebm.hm .== 1 #Bucket storage
    phi = zeros(Float64, irow, icol)
    SliqMax = zeros(Float64, irow, icol)
    for k in 1:ebm.Nsmax
        kind = cond_bucket_storage .&& k .<= ebm.Nsnow
        cond_tmp2 = kind .&& ebm.Ds[:,:,k] .> eps(Float64)
        phi[cond_tmp2] = 1 .- ebm.Sice[cond_tmp2, k] ./ (cn.rho_ice .* ebm.Ds[cond_tmp2, k])
        SliqMax[kind] = cn.rho_wat .* ebm.Ds[kind, k] .* phi[kind] .* ebm.Wirr
        ebm.Sliq[kind, k] = ebm.Sliq[kind, k] .+ Roff[kind]
        Roff[kind] .= 0
        cond_liq_excee = kind .&& ebm.Sliq[:,:,k] .> SliqMax #Liquid capacity exceeded
        Roff[cond_liq_excee] = ebm.Sliq[cond_liq_excee, k] .- SliqMax[cond_liq_excee] # so drainage to next layer
        ebm.Sliq[cond_liq_excee, k] = SliqMax[cond_liq_excee]

        coldcont[kind, k] = ebm.snow_csnow[kind, k] .* (cn.Tm .- ebm.Tsnow[kind, k])
        cond_coldcontg0 = kind .&& coldcont[:,:, k] .>0 #Liquid can freeze
        dSice[cond_coldcontg0] = min.(ebm.Sliq[cond_coldcontg0, k], coldcont[cond_coldcontg0, k] ./ cn.Lf)
        ebm.Sliq[cond_coldcontg0, k] -= dSice[cond_coldcontg0]
        ebm.Sice[cond_coldcontg0, k] += dSice[cond_coldcontg0]
        ebm.Tsnow[cond_coldcontg0, k] += cn.Lf .* dSice[cond_coldcontg0] ./ ebm.snow_csnow[cond_coldcontg0, k]
    end

    #Snow compaction
    cond_fixed_sndens = cond_ex_sp .&& ebm.dm .== 0 #fixed snow density
    for k in 1:ebm.Nsmax
        kind = cond_fixed_sndens .&& k .<= ebm.Nsnow
        ebm.Ds[kind, k] = (ebm.Sice[kind, k] .+ ebm.Sliq[kind, k]) ./ ebm.rho0
    end
    cond_prog_sndens = cond_ex_sp .&& ebm.dm .== 1 #prognostic snow density
    tau = zeros(Float64, irow, icol)
    rhos = zeros(Float64, irow, icol)
    for k in 1:ebm.Nsmax
        kind = cond_prog_sndens .&& k .<= ebm.Nsnow
        cond_Dsgeps = kind .&& ebm.Ds[:,:,k] .> eps(Float64)
        rhos[cond_Dsgeps] = (ebm.Sice[cond_Dsgeps, k] .+ ebm.Sliq[cond_Dsgeps, k]) ./ ebm.Ds[cond_Dsgeps, k]
        cond_tmp = cond_Dsgeps .&& ebm.Tsnow[:,:, k] .>= cn.Tm
        cond_tmp2 = cond_tmp .&& rhos .< ebm.rmlt
        rhos[cond_tmp2] = ebm.rmlt .+ (rhos[cond_tmp2] .- ebm.rmlt) .* exp.(-ebm.dt ./ tau[cond_tmp2])
        cond_elsetmp = cond_Dsgeps .&& ebm.Tsnow[:,:, k] .< cn.Tm
        cond_elsetmp2 = cond_elsetmp .&& rhos .< ebm.rcld
        rhos[cond_elsetmp2] = ebm.rcld .+ (rhos[cond_elsetmp2] .- ebm.rcld) .* exp.(-ebm.dt ./ tau[cond_elsetmp2])
        ebm.Ds[cond_Dsgeps, k] = (ebm.Sice[cond_Dsgeps, k] .+ ebm.Sliq[cond_Dsgeps, k]) ./ rhos[cond_Dsgeps]
    end
    #end existing snowpack

    # Add snowfall and frost to layer 1
    dSice = Sf .* ebm.dt .- min.(ebm.Esnow, 0) .* ebm.dt
    ebm.Ds[:,:, 1] += dSice ./ ebm.rfs
    ebm.Sice[:, :, 1] += dSice

    # New snowpack
    cond_new_snow = ebm.Nsnow .== 0 .&& ebm.Sice[:,:, 1] .> 0
    ebm.Nsnow[cond_new_snow] .= 1
    ebm.Tsnow[cond_new_snow, 1] = min.(Ta[cond_new_snow], cn.Tm)

    # Calculate snow depth and SWEs
    for k = 1:ebm.Nsmax
        kind = k .<= ebm.Nsnow
        snowdepth[kind] += ebm.Ds[kind, k]
        SWE[kind] += ebm.Sice[kind, k] .+ ebm.Sliq[kind, k]
    end

    # Store state of old layers
    D .= ebm.Ds
    S .= ebm.Sice
    W .= ebm.Sliq

    for k = 1:ebm.Nsmax
        kind = k .<= ebm.Nsnow
        ebm.snow_csnow[kind, k] = ebm.Sice[kind, k] .* cn.hcap_ice .+ ebm.Sliq[kind, k] .* cn.hcap_wat
        ebm.snow_E[kind, k] = ebm.snow_csnow[kind, k] .* (ebm.Tsnow[kind, k] .- cn.Tm)
    end
    Nold = ebm.Nsnow
    
    # Initialise new layers
    ebm.Ds .= 0.0
    ebm.Sice .= 0.0
    ebm.Sliq .= 0.0
    ebm.Tsnow .= cn.Tm
    ebm.snow_U .= 0.0
    ebm.Nsnow .= 0.0


    cond_ex_sp2 = SWE .> 0 #Existing or new snowpack
    #Re-assign and count snow layers
    dnew = zeros(Float64, irow, icol)
    dnew[cond_ex_sp2] = snowdepth[cond_ex_sp2]
    ebm.Ds[cond_ex_sp2, 1] = dnew[cond_ex_sp2]

    cond_DsgDz = cond_ex_sp2 .&& ebm.Ds[:,:,1] .> ebm.Dzsnow[1]
    kind2 = cond_DsgDz
    for k in 1:ebm.Nsmax
        ebm.Ds[kind2, k] .= ebm.Dzsnow[k]
        dnew[kind2] .-= ebm.Dzsnow[k]
        cond_2 = kind2 .&& (dnew .<= ebm.Dzsnow[k] .|| k.== ebm.Nsmax)
        ebm.Ds[cond_2, k] += dnew[cond_2]
        #kind2[cond_2] .= 0 #break
        for ix in eachindex(kind2)
            if cond_2[ix]
                kind2[ix] = 0
            end
        end
    end

    ebm.Nsnow[cond_ex_sp2] .= 0
    for i in 1:ebm.Nsmax
        kind = cond_ex_sp2 .&& ebm.Ds[:,:,i] .> 0
        ebm.Nsnow[kind] .+= 1
    end

    for col in 1:icol
        for row in 1:irow
            if (SWE[row, col] > 0.0)       # Existing or new snowpack

                # Fill new layers from the top downwards
                knew = 1
                dnew[row, col] = ebm.Ds[row, col, 1]
                for kold = 1:Nold[row, col]
                    while true
                        if (D[row, col, kold] < dnew[row, col])

                            # Transfer all snow from old layer and move to next old layer
                            ebm.Sice[row, col, knew] = ebm.Sice[row, col, knew] + S[row, col, kold]
                            ebm.Sliq[row, col, knew] = ebm.Sliq[row, col, knew] + W[row, col, kold]
                            ebm.snow_U[row, col, knew] += ebm.snow_E[row, col, kold]
                            dnew[row, col] -= D[row, col, kold]
                            break   ##### hack   exit
                        else
                            # Transfer some snow from old layer and move to next new layer
                            wt = dnew[row, col] / D[row, col, kold]
                            ebm.Sice[row, col, knew] = ebm.Sice[row, col, knew] + wt * S[row, col, kold]
                            ebm.Sliq[row, col, knew] = ebm.Sliq[row, col, knew] + wt * W[row, col, kold]
                            ebm.snow_U[row, col, knew] += wt * ebm.snow_E[row, col, kold]
                            D[row, col, kold] = (1 - wt) * D[row, col, kold]
                            ebm.snow_E[kold] = (1 - wt) * ebm.snow_E[row, col, kold]
                            S[row, col, kold] = (1 - wt) * S[row, col, kold]
                            W[row, col, kold] = (1 - wt) * W[row, col, kold]
                            knew += 1
                            if (knew > ebm.Nsnow[row, col])
                                break #### hack  exit
                            end
                            dnew[row, col] = ebm.Ds[row, col, knew]
                        end

                    end
                end
            end
        end
    end
    #=
    #Fill new layers from the top downwards
    knew = ones(Int64, irow, icol)
    dnew[cond_ex_sp2] = ebm.Ds[cond_ex_sp2, 1]
    for kold in 1:maximum(Nold)
        kind = cond_ex_sp2 .&& kold .<= Nold
        while true
            cond_while = kind .&& D[:,:,kold] .< dnew
            ebm.Sice[cond_while, knew] = ebm.Sice[cond_while, knew] .+ S[cond_while, kold]
            ebm.Sliq[cond_while, knew] = ebm.Sliq[cond_while, knew] .+ W[cond_while, kold]
            ebm.snow_U[cond_while, knew] = ebm.snow_U[cond_while, knew] .+ ebm.snow_E[cond_while, kold]
            dnew[cond_while] -= D[cond_while, kold]
            for ix in eachindex(kind)
                if cond_while[ix]
                    kind[ix] = 0 #break
                end
            end

            # Transfer some snow from old layer and move to next new layer
            cond_elsewhile = kind .&& D[:,:,kold] .>= dnew
            wt = dnew ./ D[:,:,kold]
            ebm.Sice[cond_elsewhile, knew] = ebm.Sice[cond_elsewhile, knew] .+ wt[cond_elsewhile] .* S[cond_elsewhile, kold]
            ebm.Sliq[cond_elsewhile, knew] = ebm.Sliq[cond_elsewhile, knew] .+ wt[cond_elsewhile] .* W[cond_elsewhile, kold]
            ebm.snow_U[cond_elsewhile, knew] += wt[cond_elsewhile] .* ebm.snow_E[cond_elsewhile, kold]
            D[cond_elsewhile, kold] = (1 .- wt[cond_elsewhile]) .* D[cond_elsewhile, kold]
            ebm.snow_E[cond_elsewhile, kold] = (1 .- wt[cond_elsewhile]) .* ebm.snow_E[cond_elsewhile, kold]
            S[cond_elsewhile, kold] = (1 .- wt[cond_elsewhile]) .* S[cond_elsewhile, kold]
            W[cond_elsewhile, kold] = (1 .- wt[cond_elsewhile]) .* W[cond_elsewhile, kold]
            knew[cond_elsewhile] .+= 1
            kind[cond_elsewhile .&& knew .> ebm.Nsnow] .= 0 #break
            dnew[cond_elsewhile] = ebm.Ds[cond_elsewhile, knew]
        end
    end
    =#
    #Diagnose snow layer temperatures
    for k in 1:ebm.Nsmax
        kind = cond_ex_sp2 .&& k .<= ebm.Nsnow
        ebm.snow_csnow[kind, k] = ebm.Sice[kind, k] .* cn.hcap_ice .+ ebm.Sliq[kind, k] .+ cn.hcap_wat
        cond_eps = kind .&& ebm.snow_csnow[:,:,k] .> eps(Float64)
        ebm.Tsnow[cond_eps, k] = cn.Tm .+ ebm.snow_U[cond_eps, k] ./ ebm.snow_csnow[cond_eps, k]
    end
    return snowdepth, SWE
end