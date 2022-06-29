function surf_ebal(ebm::EBM, cn::Constants, Ta, Qa, Ua, Ps, SW, LW, row::Int = 1, col::Int = 1)

    Qs = qsat(false, Ps, ebm.Tsurf[row, col], cn)

    psi = ebm.gs[row, col] / (ebm.gs[row, col] + ebm.CH[row, col] * Ua)
    if (Qs < Qa || ebm.Sice[row, col, 1] > 0.0)
        psi = 1.0
    end
    Lh = cn.Ls
    if (ebm.Tsurf[row, col] > cn.Tm)
        Lh = cn.Lc
    end
    rho = Ps / (cn.Rgas * Ta)
    rKH = rho * ebm.CH[row, col] * Ua

    # Surface energy balance without melt
    D = Lh * Qs / (cn.Rwat * ebm.Tsurf[row, col]^2.0)
    Esurf = psi * rKH * (Qs - Qa)
    ebm.Gsurf[row, col] = 2.0 * ebm.ksurf[row, col] * (ebm.Tsurf[row, col] - ebm.Ts1[row, col]) / ebm.Dz1[row, col]

    ebm.Hsurf[row, col] = cn.cp * rKH * (ebm.Tsurf[row, col] - Ta)
    ebm.Lesrf[row, col] = Lh * Esurf
    ebm.Melt[row, col] = 0.0
    ebm.Rnet[row, col] = (1.0 - ebm.alb[row, col]) * SW + LW - cn.sb * ebm.Tsurf[row, col]^4.0
    dTs = (ebm.Rnet[row, col] - ebm.Hsurf[row, col] - ebm.Lesrf[row, col] - ebm.Gsurf[row, col]) / ((cn.cp + Lh * psi * D) * rKH + 2 * ebm.ksurf[row, col] / ebm.Dz1[row, col] + 4.0 * cn.sb * ebm.Tsurf[row, col]^3.0)
    dE = psi * rKH * D * dTs
    dG = 2.0 * ebm.ksurf[row, col] * dTs / ebm.Dz1[row, col]
    dH = cn.cp * rKH * dTs
    dR = -cn.sb * ebm.Tsurf[row, col]^3.0 * dTs

    # Surface melting
    if (ebm.Tsurf[row, col] + dTs > cn.Tm && ebm.Sice[row, col, 1] > 0.0)
        ebm.Melt[row, col] = sum(ebm.Sice[row, col, :]) / ebm.dt
        dTs = (ebm.Rnet[row, col] - ebm.Hsurf[row, col] - ebm.Lesrf[row, col] - ebm.Gsurf[row, col] - cn.Lf * ebm.Melt[row, col]) / ((cn.cp + cn.Ls * psi * D) * rKH + 2.0 * ebm.ksurf[row, col] / ebm.Dz1[row, col] + 4.0 * cn.sb * ebm.Tsurf[row, col]^3.0)
        dE = rKH * D * dTs
        dG = 2.0 * ebm.ksurf[row, col] * dTs / ebm.Dz1[row, col]
        dH = cn.cp * rKH * dTs
        if (ebm.Tsurf[row, col] + dTs < cn.Tm)
            Qs = qsat(false, Ps, cn.Tm, cn)
            Esurf = rKH * (Qs - Qa)
            ebm.Gsurf[row, col] = 2.0 * ebm.ksurf[row, col] * (cn.Tm - ebm.Ts1[row, col]) / ebm.Dz1[row, col]
            ebm.Hsurf[row, col] = cn.cp * rKH * (cn.Tm - Ta)
            ebm.Lesrf[row, col] = cn.Ls * Esurf
            ebm.Rnet[row, col] = (1.0 - ebm.alb[row, col]) * SW + LW - cn.sb * cn.Tm^4.0
            ebm.Melt[row, col] = (ebm.Rnet[row, col] - ebm.Hsurf[row, col] - ebm.Lesrf[row, col] - ebm.Gsurf[row, col]) / cn.Lf
            ebm.Melt[row, col] = max(ebm.Melt[row, col], 0.0)
            dE = 0.0
            dG = 0.0
            dH = 0.0
            dR = 0.0
            dTs = cn.Tm - ebm.Tsurf[row, col]
        end
    end

    # Update surface temperature and fluxes
    ebm.Tsurf[row, col] = ebm.Tsurf[row, col] + dTs
    Esurf = Esurf + dE
    ebm.Gsurf[row, col] = ebm.Gsurf[row, col] + dG
    ebm.Hsurf[row, col] = ebm.Hsurf[row, col] + dH
    ebm.Rnet[row, col]  = ebm.Rnet[row, col] + dR
    ebm.Esnow[row, col] = 0.0
    Esoil = 0.0
    if (ebm.Sice[row, col, 1] > 0.0 || ebm.Tsurf[row, col] < cn.Tm)
        ebm.Esnow[row, col] = Esurf
        ebm.Lesrf[row, col] = cn.Ls * Esurf
    else
        Esoil = Esurf
        ebm.Lesrf[row, col] = cn.Lc * Esurf
    end

    return nothing

end
