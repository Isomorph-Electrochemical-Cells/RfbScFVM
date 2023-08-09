@kwdef struct CurrentCollectorParameters{T<:AbstractFloat}
    σₑ::T # electric conductivity
    hₜ::T # heat transfer coefficient
    λₜ::T # thermal conductivity
    cpᵥ::T # volumetric thermal capacity
end

@kwdef struct ReactionParameters{T<:AbstractFloat, NSPEC}
    name::String = "" # name of reaction
    Δs::T = 0.0 # entropy change
    α::T = 0.5 # symmetry coefficient
    ki::T = 0.0 # kinetic number
    ∂Δϕ₀_∂temp::T = 0.0 # temperature dependence of ki
    temp_ref::T = 1.0 # reference temperature of ki and Δϕ₀
    ν_el::IndexType # number of exchanged electrons (absolute value)
    # stoichiometric coefficients of species in the electrolyte participating in the reaction
    ν_coeff::SVector{NSPEC, IndexType} = SVector{NSPEC, IndexType}(zeros(IndexType, NSPEC))
    idx_ox::IndexType # index of oxidized species in ν_coeff
    idx_red::IndexType # index of reduced species in ν_coeff
    Δϕ₀::T # (formal) standard reduction potential
end

@kwdef struct ElectrodeParameters{T<:AbstractFloat, NSPEC}
    λₜ::T # thermal conductivity
    cpᵥ::T # thermal volumetric capacity
    εₗ::T # electrolyte volume fraction (porosity)
    kₕ::T # hydraulic permeability
    σₑ::T # electric conductivity
    aᵥ::T # specific surface area
    μ::T # dynamic viscosity
    Cᵥ::T # volumetric capacitance
    Deff_const::T # effect of pore scale geometry on effective diffusivity
    Deff_linear::T # linear dependency off mechanical dispersion on velocity
    Deff_quadratic::T # quadratic dependency of mechanical dispersion on velocity
    sh_factor::T # pre-factor of the mass-transfer model
    sh_exponent::T # exponent of the mass-transfer model
    reactions::ReactionParameters{T, NSPEC}
end

@kwdef struct SeparatorParameters{T<:AbstractFloat}
    λₜ::T # thermal conductivity
    cpᵥ::T # thermal volumetric capacity
    μ::T # dynamic viscosity
    kh::T # hydraulic permeability
    kϕ::T # electrokinetic permeability
    σₑ::T # electric conductivity
    ∂σₑ∂temp::T # derivative of conductivity with temperature
    temp_ref::T # reference temperature of model parameters
    cf::T # molar density of fixed ionic groups
    zf::T # charge value of fixed ionic groups
end

@kwdef struct BoundaryConditions{T<:AbstractFloat, NSPEC}
    # species_neg::SVector{NSPEC, T}
    # species_pos::SVector{NSPEC, T}
    species_neg::KeyedArray{T, 1, NamedDimsArray{(:row,), T, 1, SVector{NSPEC, T}}, Base.RefValue{Vector{String}}} # inlet concentrations in negative electrode
    species_pos::KeyedArray{T, 1, NamedDimsArray{(:row,), T, 1, SVector{NSPEC, T}}, Base.RefValue{Vector{String}}} # inlet concentrations in positive electrode
    temp_amb::T # ambient temperature
    p_in_neg::T # electrolyte pressure at inlet in negative electrode
    p_in_pos::T # electrolyte pressure at inlet in positive electrode
    v_out_neg::T #TODO: add support for outlet pressure
    v_out_pos::T #TODO: add support for outlet pressure
    ϕₛ_neg::Base.RefValue{T} # applied electrostatic potential at negative electrode
    ϕₛ_pos::Base.RefValue{T} # applied electrostatic potential at positive electrode
end

@kwdef struct PolarizationParameters{T<:AbstractFloat}
    voltage_start::T
    voltage_stop::T
    voltage_step::T
end

@kwdef struct StudyParameters{T<:AbstractFloat}
    polarization::PolarizationParameters{T}
    non_isothermal::Bool
    migration::Bool
    output_folder::String
    output_file_name::String
    generate_figures::Bool
end

@kwdef struct CharacteristicScales{T<:AbstractFloat}
    TEMP0::Unitful.Temperature{T} = T(273.15+25.0)u"K" # reference temperature [K]
    L0::Unitful.Length{T} = T(1e-2)*u"m" # default reference macroscopic length [m]
    LP0::Unitful.Length{T} = T(1e-5)*u"m" # reference pore-scale length [m]
    LDL0::Unitful.Length{T} = T(1e-9)*u"m" # reference double layer length [m]
    C0::Unitful.Molarity{T} = upreferred(1/(NA*LDL0^3)) #T(1e+3)u"mol/m^3" # reference concentration [mol / m^3]
    KH0::Permeability{T} = upreferred(LP0^2) # reference hydraulic permeability
    Kϕ0::Permeability{T} = upreferred(LDL0^2) # reference electrokinetic permeability
    M0::MolarMass{T} = T(1.0)*u"kg/mol" # reference molar mass
    ρ0::Unitful.Density{T} = upreferred(C0*M0) # reference mass density [kg / m^3]
    λ0::HeatConductivity{T} = T(1.0)u"W/m/K" # reference heat conductivity [W / (m K)]
    S0::MolarEntropy{T} = upreferred(R) # reference molar entropy [J/(K mol)]
    V0::Unitful.Voltage{T} = upreferred(R*TEMP0/F) # reference voltage [V]
    μ0::Unitful.DynamicViscosity{T} = T(1e-3)u"Pa*s" # reference dynamic viscosity [kg/(m s) = Pa⋅s]
    ν0::Unitful.KinematicViscosity{T} = upreferred(μ0 / ρ0) # reference kinematic viscosity [m^2 / s]
    D0::Diffusivity{T} = upreferred(R*TEMP0 / (μ0*LDL0*NA)) #T(1e-9)u"m^2/s" # reference diffusion coefficient [m^2 / s]
    VEL0::Unitful.Velocity{T} = T(1e-3)u"m/s" #upreferred(KP0*P0/(μ0*L0)) # reference velocity [m/s]
    P0::Unitful.Pressure{T} = upreferred(L0*μ0*VEL0/KH0) #upreferred(R*TEMP0*C0) #T(1e+3)u"J/m^3" # reference pressure [J / m^3]
    T0::Unitful.Time{T} = upreferred(L0/VEL0) # reference time [s]
    σ0::Unitful.ElectricalConductivity{T} = upreferred((C0*D0*F)/(V0)) # reference conducitivy [Siemens/m]
    i0::CurrentDensity{T} = upreferred((C0*D0*F)/LP0) # reference current density [A/m^2]
    iv0::VolumetricCurrentDensity{T} = upreferred(i0/LP0) # reference volumetric current density [A/m^3]
    CP0::SpecificHeatCapacity{T} = upreferred(R/M0) # reference specific heat capacity [J / (kg K)]
    CPv0::VolumetricHeatCapacity{T} = upreferred(CP0*ρ0) # reference volumetric heat capacity [J / (m^3 K)]
    Ca0::ArealCapacitance{T} = upreferred(C0*F*LDL0/V0) # reference areal capacitance [F/m^2]
    Cv0::VolumetricCapacitance{T} = upreferred(Ca0/LP0) # reference volumetric capacitance [F/m^3]
    K0::HeterogeneousReactionRate{T} = upreferred(D0/LP0) # reference heterogeneous reaction rate [m/s]
end

@kwdef struct ScalingParameters{T<:AbstractFloat}
    ϵL0::T # ratio of pore-scale length to macroscopic length scale
    ϵl0::T # ratio of double layer length scale to pore-scale length
    PR0::T # Prandtl number
    SC0::T # Schmidt number
    RE0::T # Reynolds number
    KI0::T # Kinetic number
    PE0::T # Peclet number
    LE0::T # Lewis number
end

@kwdef struct DiscretizationParameters{T<:AbstractFloat}
    spatial_discretization::String
    temporal_discretization::String
    relative_tolerance::T
end

function ScalingParameters{T}(scales::CharacteristicScales{T}) where {T<:AbstractFloat}
    ϵL0::T = scales.LP0 / scales.L0 |> NoUnits # ratio of pore-scale length to macroscopic length
    ϵl0::T = scales.LDL0 / scales.LP0 |> NoUnits # ratio of pore-scale length to macroscopic length
    PR0::T = (scales.μ0 * scales.CP0) / (scales.λ0) |> NoUnits # Prandtl number
    SC0::T = scales.ν0 / scales.D0 |> NoUnits # Schmidt number
    RE0::T = scales.LP0*scales.VEL0/scales.ν0 |> NoUnits # Reynolds number
    KI0::T = scales.K0*scales.LP0 / scales.D0 |> NoUnits # Kinetic number
    PE0::T = scales.L0*scales.VEL0 / scales.D0 |> NoUnits # Peclet number
    LE0::T = SC0/PR0 # Lewis number
    ScalingParameters{T}(ϵL0=ϵL0, ϵl0=ϵl0, PR0=PR0, SC0=SC0, RE0=RE0,
                               KI0=KI0, PE0=PE0, LE0=LE0)
end

@kwdef struct ElectrolyteParams{T<:AbstractFloat, NSPEC, NATTR, NSPEC_TIMES_NATTR}
    # Electrolyte species without the counter species
    species::KeyedArray{T, 2, NamedDimsArray{(:row, :col), T, 2, MMatrix{NSPEC, NATTR, T, NSPEC_TIMES_NATTR}}, Tuple{Vector{String}, Vector{String}}}
    # Counter species, whose concentration is determined by the strong electroneurality condition
    counter::KeyedArray{T, 1, NamedDimsArray{(:col,), T, 1, MVector{NATTR, T}}, Base.RefValue{Vector{String}}}
    λₜ::T # thermal conductivity
    cpᵥ::T # thermal volumetric capacity
end

@kwdef mutable struct SystemResults{T<:AbstractFloat}
    current_density::Vector{T} = Vector{T}(undef,0)
end


@kwdef struct DomainVariables{T<:Real, ArrayType<:AbstractArray}
    p::T
    ϕₛ::T
    ϕₗ::T
    c::ArrayType
    temp::T
end

@kwdef struct ModelParameters{T<:AbstractFloat, NSPEC, NSPEC_SYSTEM, NATTR, NSPEC_TIMES_NATTR,
              NAMED_TUPLE_TYPE, DOMAIN_ID_TYPE, BOUNDARY_ID_TYPE}
    geom::FlowCellGeometry{T}
    mesh::Mesh2D{T}
    el_neg::ElectrodeParameters{T, NSPEC}
    el_pos::ElectrodeParameters{T, NSPEC}
    cc_neg::CurrentCollectorParameters{T}
    cc_pos::CurrentCollectorParameters{T}
    sep::SeparatorParameters{T}
    boundary::BoundaryConditions{T, NSPEC_SYSTEM}
    discr::DiscretizationParameters{T}
    study::StudyParameters{T}
    scales::CharacteristicScales{T}
    scaling_params::ScalingParameters{T}
    electrolyte::ElectrolyteParams{T, NSPEC_SYSTEM, NATTR, NSPEC_TIMES_NATTR}
    idx::NAMED_TUPLE_TYPE
    results::SystemResults{T}
    dom::DOMAIN_ID_TYPE
    bnd::BOUNDARY_ID_TYPE
end
