function tridiag(Nvec, gamma, a, b, c, r, x)

    gamma .= 0.0

    beta = b[:,:,1]
    x[:,:,1] = r[:,:,1] ./ beta

    if isa(Nvec, Array)
        for n = 2:maximum(Nvec)
            kind = Nvec .<= n
            gamma[kind,n] = c[kind,n - 1] ./ beta
            beta = b[kind, n] .- a[kind, n] .* gamma[kind, n]
            x[kind, n] = (r[kind, n] .- a[kind, n] .* x[kind, n - 1]) ./ beta
        end

        for n = (maximum(Nvec) - 1):-1:1
            kind = Nvec[:,:] .<= n
            x[kind,n] = x[kind,n] .- gamma[kind,n + 1] .* x[kind,n + 1]
        end
    else
        for n = 2:Nvec
            gamma[:,:,n] = c[:,:,n - 1] ./ beta
            beta = b[:,:, n] .- a[:,:, n] .* gamma[:,:, n]
            x[:,:, n] = (r[:,:, n] .- a[:,:, n] .* x[:,:, n - 1]) ./ beta
        end

        for n = (Nvec -1):-1:1
            x[:,:,n] = x[:,:,n] .- gamma[:,:,n + 1] .* x[:,:,n + 1]
        end
    end

    return nothing

end
