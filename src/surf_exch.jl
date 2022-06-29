function surf_exch(ebm::EBM, cn::Constants, Ta, Ua, row::Int = 1, col::Int = 1)

    zTs = ebm.zT
    if (ebm.zvar)
        zTs = ebm.zT - sum(ebm.Ds[row, col, :])
        zTs = max(zTs, 1.0)
    end

    # Neutral exchange coefficients
    z0h = 0.1 * ebm.z0[row, col]
    CD = (cn.vkman / log(ebm.zU / ebm.z0[row, col]))^2.0
    ebm.CH[row, col] = cn.vkman^2.0 / (log(ebm.zU / ebm.z0[row, col]) * log(zTs / z0h))

    # Stability correction (Louis et al. 1982, quoted by Beljaars 1992)
    if (ebm.em == 1)
        RiB = cn.g * (Ta - ebm.Tsurf[row, col]) * ebm.zU^2.0 / (zTs * Ta * Ua^2.0)
        if (RiB > 0.0)
            fh = 1.0 / (1.0 + 3.0 * ebm.bstb * RiB * sqrt(1.0 + ebm.bstb * RiB))
        else
            fh = 1.0 - 3.0 * ebm.bstb * RiB / (1.0 + 3.0 * ebm.bstb^2.0 * CD * sqrt(-RiB * ebm.zU / ebm.z0[row, col]))
        end
        ebm.CH[row, col] = fh * ebm.CH[row, col]
    end

    return nothing

end
