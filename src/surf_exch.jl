function surf_exch(ebm::EBM, cn::Constants, Ta, Ua, irow::Int = 1, icol::Int = 1)

    zTs = fill(ebm.zT, irow, icol)
    if (ebm.zvar)
        zTs = ebm.zT .- dropdims(sum(ebm.Ds, dims=3), dims=3)
        zTs = max.(zTs, 1)
    end

    # Neutral exchange coefficients
    z0h = 0.1 .* ebm.z0
    CD = (cn.vkman ./ log.(ebm.zU ./ ebm.z0)).^2
    ebm.CH = cn.vkman^2 ./ (log.(ebm.zU ./ ebm.z0) .* log.(zTs ./ z0h))
   
    # Stability correction (Louis et al. 1982, quoted by Beljaars 1992)
    stabcorr_cond = ebm.em .== 1
    RiB = cn.g .* (Ta[stabcorr_cond] .- ebm.Tsurf[stabcorr_cond]) .* ebm.zU^2 ./ (zTs[stabcorr_cond] .* Ta[stabcorr_cond] .* Ua[stabcorr_cond].^2)
    fh = similar(RiB)
    cond_RiB_g0 = RiB .> 0
    fh[cond_RiB_g0] = 1 ./ (1 .+ 3 .* ebm.bstb .* RiB[cond_RiB_g0] .* sqrt.(1 .+ ebm.bstb .* RiB[cond_RiB_g0]))
    fh[.!(cond_RiB_g0)] = 1 .- 3 .* ebm.bstb .* RiB[.!(cond_RiB_g0)] ./ (1 .+ 3 .* ebm.bstb.^2 .* CD[.!(cond_RiB_g0)] .* sqrt.(-1 .* RiB[.!(cond_RiB_g0)] .* ebm.zU ./ ebm.z0[.!(cond_RiB_g0)]))
    ebm.CH = reshape(fh, irow, icol) .* ebm.CH
    
    return nothing

end
