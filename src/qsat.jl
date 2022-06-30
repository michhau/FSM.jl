function qsat(water, P, T, cn::Constants)
    Tc = T .- cn.Tm
    if water
        es = cn.e0 .* exp.(17.5043 .* Tc ./ (241.3 .+ Tc))
        Qs = cn.eps .* es ./ P
        return Qs
    end
    if length(Tc) > 1
        es = similar(Tc)
        TCgr0 = Tc .> 0
        es[TCgr0] .= cn.e0 .* exp.(17.5043 .* Tc[TCgr0] ./ (241.3 .+ Tc[TCgr0]))
        es[.!(TCgr0)] .= cn.e0 .* exp.(22.4422 * Tc[.!(TCgr0)] ./ (272.186 .+ Tc[.!(TCgr0)]))
        Qs = cn.eps .* es ./ P
    else
        es = 0.0
        if Tc > 0
            es = cn.e0 * exp.(17.5043 * Tc / (241.3 + Tc))
        else
            es = cn.e0 * exp.(22.4422 * Tc / (272.186 + Tc))
        end
        Qs = cn.eps .* es ./ P
    end
    
    return Qs

end
