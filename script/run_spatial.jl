using FSM
using Plots
using ProgressMeter
using MAT
using Dates#, BenchmarkTools

if gethostname() == "LINUX24"
    meteosource = "/home/haugened/Documents/data/FSM_input/spatial"
elseif gethostname() == "Michi-T450s"
    meteosource = "/home/michi/Documents/slf/data_input_FSM/spatial/"
else 
    meteosource = "K:/DATA_COSMO/OUTPUT_GRID_OSHD_0250/PROCESSED_ANALYSIS/COSMO_1EFA"
end

# Helper functions

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

function compute_psolid(ptot, ta, thres_prec=274.19, m_prec=0.1500)

    p_corr = 1.0
    tp = @. (ta - thres_prec) / m_prec
    p_multi = @. p_corr / (1 + exp(tp))

    return @. p_multi * ptot

end

function compute_pliquid(ptot, ta, thres_prec=274.19, m_prec=0.1500)

    p_corr = 1.0

    Tp = @. (ta - thres_prec) / m_prec
    p_multi = @. p_corr * exp(Tp) / (1 + exp(Tp))
    return @. p_multi * ptot

end

function read_meteo(t, source)

    folder = joinpath(source, Dates.format(t, "yyyy.mm"))
    filename = searchdir(folder, "COSMODATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")

    meteo = matread(joinpath(folder, filename[1]))

    Ta = meteo["tais"]["data"]
    RH = meteo["rhus"]["data"]
    Ua = meteo["wnsc"]["data"]
    SW = meteo["sdrd"]["data"] .+ meteo["sdfd"]["data"]
    LW = meteo["lwrc"]["data"]
    Sf = compute_psolid(meteo["prcs"]["data"], meteo["tais"]["data"]) ./ 3600
    Rf = compute_pliquid(meteo["prcs"]["data"], meteo["tais"]["data"]) ./ 3600
    Ps = meteo["pail"]["data"]

    return Ta, RH, Ua, SW, LW, Sf, Rf, Ps

end

# Setting up model

times = DateTime(2021,11,01,01):Hour(1):DateTime(2021,11,02,01)
#DateTime(2020, 9, 1, 1):Hour(1):DateTime(2021, 7, 1, 1)

nrows = 1088
ncols = 1460

function setupspatialrun(nrows::Int, ncols::Int)
    @info("Setting up spatial run")
    tsoil_tmp = zeros(Float64, nrows, ncols, 4)
    for icol in 1:ncols
        for jrow in 1:nrows
            tsoil_tmp[jrow, icol, :] .= [282.98, 284.17, 284.70, 284.70]
        end
    end

    ebm_mat = EBM{Float64}(
            am=ones(Int64, nrows, ncols),
            cm=ones(Int64, nrows, ncols),
            dm=ones(Int64, nrows, ncols),
            em=ones(Int64, nrows, ncols),
            hm=ones(Int64, nrows, ncols),
            zT=10.5,
            zvar=false,
            Tsoil=tsoil_tmp
        )

    cn = Constants{Float64}()

    snowdepth = ones(nrows, ncols)
    SWE = ones(nrows, ncols)
    Tsurf = ones(nrows, ncols)

    return ebm_mat, cn, snowdepth, SWE, Tsurf
end

@time (ebm_mat, cn, snowdepth, SWE, Tsurf) = setupspatialrun(nrows, ncols)

# Run the model

@showprogress "Running model..." for t in times

    Ta, RH, Ua, SW, LW, Sf, Rf, Ps = read_meteo(t, meteosource)

    run!(ebm_mat, cn, snowdepth, SWE, Tsurf, SW, LW, Sf, Rf, Ta, RH, Ua, Ps)

  #=  if hour(t) == 0
        file = matopen("D:/FSMJL/grid/" * Dates.format(t, "yyyymmdd") * "_hs.mat", "w")
        write(file, "snowdepth", snowdepth)
        close(file)
    end
    =#
end
