struct Input{T}
    year::Vector{Int}
    month::Vector{Int}
    day::Vector{Int}
    hour::Vector{Int}
    SW::Vector{T}
    LW::Vector{T}
    Sf::Vector{T}
    Rf::Vector{T}
    Ta::Vector{T}
    RH::Vector{T}
    Ua::Vector{T}
    Ps::Vector{T}
end


Base.@kwdef struct Constants{T}

    cp::T = 1005          # Specific heat capacity of dry air (J/K/kg)
    eps::T = 0.622       # Ratio of molecular weights of water and dry air    BEFORE eps
    e0::T = 611.213       # Saturation vapour pressure at Tm (Pa)
    g::T = 9.81           # Acceleration due to gravity (m/s^2)
    hcap_ice::T = 2100.   # Specific heat capacity of ice (J/K/kg)
    hcap_wat::T = 4180.   # Specific heat capacity of water (J/K/kg)
    hcon_air::T = 0.025   # Thermal conductivity of air (W/m/K)
    hcon_clay::T = 1.16   # Thermal conductivity of clay (W/m/K)
    hcon_ice::T = 2.24    # Thermal conducivity of ice (W/m/K)
    hcon_sand::T = 1.57   # Thermal conductivity of sand (W/m/K)
    hcon_wat::T = 0.56    # Thermal conductivity of water (W/m/K)
    Lc::T = 2.501e6       # Latent heat of condensation (J/kg)
    Lf::T = 0.334e6       # Latent heat of fusion (J/kg)
    Ls::T = Lc + Lf       # Latent heat of sublimation (J/kg)
    R::T = 8.3145         # Universal gas constant (J/K/mol)
    Rgas::T = 287         # Gas constant for dry air (J/K/kg)
    Rwat::T = 462         # Gas constant for water vapour (J/K/kg)
    rho_ice::T = 917.     # Density of ice (kg/m^3)
    rho_wat::T = 1000.    # Density of water (kg/m^3)
    sb::T = 5.67e-8       # Stefan-Boltzmann constant (W/m^2/K^4)
    Tm::T = 273.15        # Melting point (K)
    vkman::T = 0.4        # Von Karman constant

end


Base.@kwdef mutable struct EBM{Float64}

    nrow = 1088
    ncol = 1460
    
    ### Settings ##############################################################

    # Maximum number of snow layers
    Nsmax::Int = 3

    # Number of soil layers
    Nsoil::Int = 4

    # Time step (s)
    dt::Float64 = 3600

    # Minimum snow layer thicknesses (m)
    Dzsnow::Vector{Float64} = [0.1, 0.2, 0.4]

    # Soil layer thicknesses (m)
    Dzsoil::Vector{Float64} = [0.1, 0.2, 0.4, 0.8]

    # Temperature measurement height (m)
    zT::Float64 = 2

    # Wind measurement height (m)
    zU::Float64 = 10

    # Subtract snow depth from measurement height
    zvar::Bool = true

    ### Snow parameters #######################################################

    # Maximum albedo for fresh snow
    asmx::Float64 = 0.8

    # Minimum albedo for melting snow
    asmn::Float64 = 0.5

    # Stability slope parameter
    bstb::Float64 = 5

    # Snow thermal conductivity exponent
    bthr::Float64 = 2

    # Snow cover fraction depth scale (m)
    hfsn::Float64 = 0.1

    # Fixed thermal conductivity of snow (W/m/K)
    kfix::Float64 = 0.24

    # Fixed snow density (kg/m^3)
    rho0::Float64 = 300

    # Fresh snow density (kg/m^3)
    rhof::Float64 = 100

    # Maximum density for cold snow (kg/m^3)
    rcld::Float64 = 300

    # Maximum density for melting snow (kg/m^3)
    rmlt::Float64 = 500

    # Snowfall to refresh albedo (kg/m^2)
    Salb::Float64 = 10

    # Albedo decay temperature threshold (C)
    Talb::Float64 = -2

    # Cold snow albedo decay timescale (h)
    tcld::Float64 = 1000

    # Melting snow albedo decay timescale (h)
    tmlt::Float64 = 100

    # Snow compaction time scale (h)
    trho::Float64 = 200

    # Irreducible liquid water content of snow
    Wirr::Float64 = 0.03

    # Snow roughness length (m)
    z0sn::Float64 = 0.01

    ### Surface parameters ####################################################

    # Snow-free ground albedo
    alb0::Float64 = 0.2

    # Surface conductance for saturated soil (m/s)
    gsat::Float64 = 0.01

    # Snow-free roughness length (m)
    z0sf::Float64 = 0.1

    ### Soil properties #######################################################

    # Soil clay fraction
    fcly::Float64 = 0.3

    # Soil sand fraction
    fsnd::Float64 = 0.6

    # Clapp-Hornberger exponent
    b::Float64 = 7.630000000000001

    # Volumetric heat capacity of dry soil (J/K/m^3)
    hcap_soil::Float64 = 2.2993333333333335e6

    # Saturated soil water pressure (m)
    sathh::Float64 = 0.10789467222298288

    # Volumetric soil moisture concentration at saturation
    Vsat::Float64 = 0.4087

    # Volumetric soil moisture concentration at critical point
    Vcrit::Float64 = 0.2603859120641063

    # Thermal conductivity of dry soil (W/m/K)
    hcon_soil::Float64 = 0.2740041303112452

    ### State variables #######################################################

    # Snow albedo
    albs::Array{Float64,2} = fill(0.8, nrow, ncol)

    # Snow layer thicknesses (m)
    Ds = zeros(Float64, nrow, ncol, Nsmax)

    # Number of snow layers
    Nsnow = zeros(Int64, nrow, ncol)

    # Ice content of snow layers (kg/m^2)
    Sice = zeros(Float64, nrow, ncol, Nsmax)

    # Liquid content of snow layers (kg/m^2)
    Sliq = zeros(Float64, nrow, ncol, Nsmax)

    # Volumetric moisture content of soil layers
    theta::Array{Float64, 3} = fill(0.5*Vsat, nrow, ncol, Nsoil)

    # Snow layer temperatures (K)
    Tsnow = fill(273.15, nrow, ncol, Nsmax)

    # Soil layer temperatures (K)
    Tsoil::Array{Float64, 3} = fill(285.0, nrow, ncol, Nsoil)

    # Surface skin temperature (K) (needs to be filled after construction)
    Tsurf = Tsoil[:,:,1]
    #for icol in 1:ncol
    #    for jrow in 1:nrow
    #        Tsurf[jrow, icol] = Tsoil[jrow, icol, 1]
    #    end
    #end

    ### Configurations ########################################################

    # Snow albedo model
    am = zeros(Int64, nrow, ncol)

    # Snow conductivity model
    cm = zeros(Int64, nrow, ncol)

    # Snow density model
    dm = zeros(Int64, nrow, ncol)

    # Surface exchange model
    em = zeros(Int64, nrow, ncol)

    # Snow hydraulics model
    hm = zeros(Int64, nrow, ncol)

    ### Internal variables ####################################################

    # Thermal conductivity of snow (W/m/K)
    ksnow = zeros(Float64, nrow, ncol, Nsmax)

    # Thermal conductivity of soil (W/m/K)
    ksoil = zeros(Float64, nrow, ncol, Nsoil)

    # Areal heat capacity of soil (J/K/m^2)
    csoil= zeros(Float64, nrow, ncol, Nsoil)

    # Surface moisture conductance (m/s)
    gs = zeros(Float64, nrow, ncol)

    # Effective albedo
    alb = zeros(Float64, nrow, ncol)

    # Fresh snow density (kg/m^3)
    rfs = zeros(Float64, nrow, ncol)

    # Snow cover fraction
    fsnow = zeros(Float64, nrow, ncol)

    # Transfer coefficient for heat and moisture
    CH = zeros(Float64, nrow, ncol)

    # Surface roughness length (m)
    z0 = zeros(Float64, nrow, ncol)

    # Surface thermal conductivity (W/m/K)
    ksurf = zeros(Float64, nrow, ncol)

    # Surface layer temperature (K)
    Ts1 = zeros(Float64, nrow, ncol)

    # Surface layer thickness (m)
    Dz1 = zeros(Float64, nrow, ncol)

    # Snow sublimation rate (kg/m^2/s)
    Esnow = zeros(Float64, nrow, ncol)

    # Heat flux into surface (W/m^2)
    Gsurf = zeros(Float64, nrow, ncol)

    # Sensible heat flux (W/m^2)
    Hsurf = zeros(Float64, nrow, ncol)

    # Latent heat flux (W/m^2)
    Lesrf = zeros(Float64, nrow, ncol)

    # Surface melt rate (kg/m^2/s)
    Melt = zeros(Float64, nrow, ncol)

    # Net radiation (W/m^2)
    Rnet = zeros(Float64, nrow, ncol)

    # Heat flux into soil (W/m^2)
    Gsoil = zeros(Float64, nrow, ncol)
    #=
    # Temporary soil variables
    soil_a::Vector{T} = fill(0, Nsoil)
    soil_b::Vector{T} = fill(0, Nsoil)
    soil_c::Vector{T} = fill(0, Nsoil)
    soil_Gs::Vector{T} = fill(0, Nsoil)
    soil_rhs::Vector{T} = fill(0, Nsoil)
    soil_dTs::Vector{T} = fill(0, Nsoil)
    soil_gamma::Vector{T} = fill(0, Nsoil)

    # Temporary snow variables
    snow_csnow::Vector{T} = fill(0, Nsmax)
    snow_E::Vector{T} = fill(0, Nsmax)
    snow_U::Vector{T} = fill(0, Nsmax)
    snow_Gs::Vector{T} = fill(0, Nsmax)
    snow_dTs::Vector{T} = fill(0, Nsmax)
    snow_a::Vector{T} = fill(0, Nsmax)
    snow_b::Vector{T} = fill(0, Nsmax)
    snow_c::Vector{T} = fill(0, Nsmax)
    snow_rhs::Vector{T} = fill(0, Nsmax)
    snow_gamma::Vector{T} = fill(0, Nsmax)
    =#
end
