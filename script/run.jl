using FSM
using CSV
using DataFrames

data_force = CSV.File("../data/met_CdP_0506.csv") |> DataFrame

input = Input{Float32}(
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
    data_force.Ps,
    )

ebm = EBM{Float32}()

cn = Constants{Float32}()

snowdepth = similar(input.Ta)
SWE = similar(input.Ta)
Tsurf = similar(input.Ta)

run!(ebm, cn, snowdepth, SWE, Tsurf, input)
