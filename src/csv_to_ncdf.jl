###########################################
###      BENCHMARK/PROFILE FSM.JL       ###
###########################################
#=
run runtests.jl before to see, if code
works properly.

needed packages ('import Pkg; Pkg.add($PKGNAME)'):
- BenchmarkTools
=#

using BenchmarkTools, CSV, DataFrames, FSM
"Read files to perform input. Content taken from runtests.jl"
function readinput()
    data_force = CSV.File(joinpath(pwd(), "data", "met_CdP_0506.csv")) |> DataFrame

    input = Input{Float64}(
        data_force.year,
        data_force.month,
        data_force.day,
        data_force.hour,
        data_force.SW,
        data_force.LW,
        data_force.Sf,
        data_force.Rf,
        data_force.Ta,
        data_force.RH,
        data_force.Ua,
        data_force.Ps,)

    path_ref = joinpath(pwd(), "fortran", "output/") #../fortran/output/"
    files_ref = readdir(path_ref)
    return input, files_ref
end

"Perform the FSM run as in runtests.jl"
function performFSM(input, files_ref)
    for file_ref in files_ref

        # Initilize model

        ebm = EBM{Float64}(
            am=parse(Int, file_ref[14]),
            cm=parse(Int, file_ref[15]),
            dm=parse(Int, file_ref[16]),
            em=parse(Int, file_ref[17]),
            hm=parse(Int, file_ref[18]),
            zT=1.5,
            zvar=false,
            Tsoil=[282.98, 284.17, 284.70, 284.70]
        )

        cn = Constants{Float64}()

        snowdepth = similar(input.Ta)
        SWE = similar(input.Ta)
        Tsurf = similar(input.Ta)

        # Run model

        run!(ebm, cn, snowdepth, SWE, Tsurf, input)
    end
end

(input, files_ref) = readinput()
performFSM(input, files_ref)

br = @benchmark readinput()
bp = @benchmark performFSM(input, files_ref)
