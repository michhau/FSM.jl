function soil(ebm::EBM, irow::Int = 1, icol::Int = 1)

    for k in 1:(ebm.Nsoil - 1)
        ebm.soil_Gs[:,:, k] = 2 ./ (ebm.Dzsoil[k] ./ ebm.ksoil[:,:, k] .+ ebm.Dzsoil[k .+ 1] ./ ebm.ksoil[:,:, k .+ 1])
    end

    ebm.soil_a[:,:,1] .= 0.0
    ebm.soil_b[:,:,1] = ebm.csoil[:,:,1] .+ ebm.soil_Gs[:,:,1] .* ebm.dt
    ebm.soil_c[:,:,1] = -1 .* ebm.soil_Gs[:,:,1] .* ebm.dt
    ebm.soil_rhs[:,:,1] = (ebm.Gsoil .- ebm.soil_Gs[:,:,1] .* (ebm.Tsoil[:,:, 1] .- ebm.Tsoil[:,:, 2])) .* ebm.dt

    k = 2:(ebm.Nsoil - 1)
    ebm.soil_a[:,:,k] = ebm.soil_c[:,:,k .- 1]
    ebm.soil_b[:,:,k] = ebm.csoil[:,:, k] .+ (ebm.soil_Gs[:,:,k .- 1] .+ ebm.soil_Gs[:,:,k]) .* ebm.dt
    ebm.soil_c[:,:,k] = -1 .* ebm.soil_Gs[:,:,k] .* ebm.dt
    ebm.soil_rhs[:,:,k] = ebm.soil_Gs[:,:,k .- 1] .* (ebm.Tsoil[:,:, k .- 1] .- ebm.Tsoil[:,:, k]) .* ebm.dt .+ ebm.soil_Gs[:,:,k] .* (ebm.Tsoil[:,:, k .+ 1] .- ebm.Tsoil[:,:, k]) .* ebm.dt
    

    k = ebm.Nsoil
    ebm.soil_Gs[:,:,k] = ebm.ksoil[:,:, k] ./ ebm.Dzsoil[k]
    ebm.soil_a[:,:,k] = ebm.soil_c[:,:,k - 1]
    ebm.soil_b[:,:,k] = ebm.csoil[:,:, k] .+ (ebm.soil_Gs[:,:,k - 1] .+ ebm.soil_Gs[:,:,k]) .* ebm.dt
    ebm.soil_c[:,:,k] .= 0.0
    ebm.soil_rhs[:,:,k] = ebm.soil_Gs[:,:,k - 1] .* (ebm.Tsoil[:,:, k - 1] .- ebm.Tsoil[:,:, k]) .* ebm.dt
    tridiag(ebm.Nsoil, ebm.soil_gamma, ebm.soil_a, ebm.soil_b, ebm.soil_c, ebm.soil_rhs, ebm.soil_dTs)

    k = 1:ebm.Nsoil
    ebm.Tsoil[:,:, k] = ebm.Tsoil[:,:, k] + ebm.soil_dTs[:,:,k]

    return  nothing

end