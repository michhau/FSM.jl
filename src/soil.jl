function soil(ebm::EBM, row::Int = 1, col::Int = 1)

    a = zeros(Float64, ebm.Nsoil)
    b = zeros(Float64, ebm.Nsoil)
    c = zeros(Float64, ebm.Nsoil)
    Gs = zeros(Float64, ebm.Nsoil)
    rhs = zeros(Float64, ebm.Nsoil)
    dTs = zeros(Float64, ebm.Nsoil)
    soil_gamma = zeros(Float64, ebm.Nsoil)

    for k = 1:(ebm.Nsoil - 1)
        Gs[k] = 2.0 / (ebm.Dzsoil[k] / ebm.ksoil[row, col, k] + ebm.Dzsoil[k + 1] / ebm.ksoil[row, col, k + 1])
    end
    a[1] = 0.0
    b[1] = ebm.csoil[row, col, 1] + Gs[1] * ebm.dt
    c[1] = - Gs[1] * ebm.dt
    rhs[1] = (ebm.Gsoil[row, col] - Gs[1] * (ebm.Tsoil[row, col, 1] - ebm.Tsoil[row, col, 2])) * ebm.dt
    for k = 2:(ebm.Nsoil - 1)
        a[k] = c[k - 1]
        b[k] = ebm.csoil[row, col, k] + (Gs[k - 1] + Gs[k]) * ebm.dt
        c[k] = - Gs[k] * ebm.dt
        rhs[k] = Gs[k - 1] * (ebm.Tsoil[row, col, k - 1] - ebm.Tsoil[row, col, k]) * ebm.dt + Gs[k] * (ebm.Tsoil[row, col, k + 1] - ebm.Tsoil[row, col, k]) * ebm.dt
    end
    k = ebm.Nsoil
    Gs[k] = ebm.ksoil[row, col, k] / ebm.Dzsoil[k]
    a[k] = c[k - 1]
    b[k] = ebm.csoil[row, col, k] + (Gs[k - 1] + Gs[k]) * ebm.dt
    c[k] = 0.0
    rhs[k] = Gs[k - 1] * (ebm.Tsoil[row, col, k - 1] - ebm.Tsoil[row, col, k]) * ebm.dt
    tridiag(ebm.Nsoil, soil_gamma, a, b, c, rhs, dTs)
    for k = 1:ebm.Nsoil
        ebm.Tsoil[row, col, k] = ebm.Tsoil[row, col, k] + dTs[k]
    end

    return  nothing

end
