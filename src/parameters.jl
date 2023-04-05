@with_kw struct CurrentCollectorParameters{T}
    σₑ::T # electric conductivity
    hₜ::T # heat transfer coefficient
    λₜ::T # thermal conductivity
    cpᵥ::T # volumetric thermal capacity
end

@with_kw struct ReactionParameters{T}
    name::String = "" # name of reaction
    Δs::T = 0.0 # entropy change
    α::T = 0.5 # symmetry coefficient
    ki::T = 0.0 # kinetic number
    ∂Δϕ₀_∂temp::T = 0.0 # temperature dependence of ki
    temp_ref::T = 0.0 # reference temperature of ki and Δϕ₀
    ν_ox::Int64 = -1 # stoichiometric coefficient of the oxidized ion
    ν_red::Int64 = 1 # stoichiometric coefficient of the reduced ion
    ν_el::Int64 = -1 # stoichiometric coefficient of the transferred electrons
    Δϕ₀::T = 0.0 # (formal) standard reduction potential
end

@with_kw struct ElectrodeParameters{T}
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
    reactions::ReactionParameters{T}
end

@with_kw struct SeparatorParameters{T}
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

@with_kw struct BoundaryConditions{T}
    species_neg::AxisArray{T, 2, Matrix{T}, Tuple{AxisArrays.Axis{:row, Vector{String}}, AxisArrays.Axis{:col, Vector{String}}}} =
    AxisArray(Matrix{T}(undef,5,1);
            row=["solvent", "ox_neg", "red_neg", "ox_pos", "red_pos"],
            col=["concentration"])
    species_pos::AxisArray{T, 2, Matrix{T}, Tuple{AxisArrays.Axis{:row, Vector{String}}, AxisArrays.Axis{:col, Vector{String}}}} =
            AxisArray(Matrix{T}(undef,5,1);
                    row=["solvent", "ox_neg", "red_neg", "ox_pos", "red_pos"],
                    col=["concentration"])
    temp_amb::T # ambient temperature
    p_in_neg::T # electrolyte pressure at inlet in negative electrode
    p_in_pos::T # electrolyte pressure at inlet in positive electrode
    v_out_neg::T #TODO: add support for outlet pressure
    v_out_pos::T #TODO: add support for outlet pressure
    ϕₛ_neg::T # applied electrostatic potential at negative electrode
    ϕₛ_pos::T # applied electrostatic potential at positive electrode
end

@with_kw struct PolarizationParameters{T}
    voltage_start::T
    voltage_stop::T
    voltage_step::T
    output_folder::String
end
@with_kw struct StudyParameters{T}
    polarization::PolarizationParameters{T}
    non_isothermal::Bool
    migration::Bool
end

@with_kw struct CharacteristicScales{T}
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

@with_kw struct ScalingParameters{T}
    ϵL0::T # ratio of pore-scale length to macroscopic length scale
    ϵl0::T # ratio of double layer length scale to pore-scale length
    PR0::T # Prandtl number
    SC0::T # Schmidt number
    RE0::T # Reynolds number
    KI0::T # Kinetic number
    PE0::T # Peclet number
    LE0::T # Lewis number
end

@with_kw struct DiscretizationParameters{T}
    spatial_discretization::String
    temporal_discretization::String
    relative_tolerance::T
end

function ScalingParameters{T}(scales::CharacteristicScales{T}) where T
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


@with_kw struct ElectrolyteParams{T}
    species::AxisArray{T, 2, Matrix{T}, Tuple{AxisArrays.Axis{:row, Vector{String}}, AxisArrays.Axis{:col, Vector{String}}}} =
            AxisArray(Matrix{T}(undef, 6, 4);
                    row=["solvent", "ox_neg", "red_neg", "ox_pos", "red_pos", "counter"],
                    col=["charge", "molar_mass", "diffusivity"])
    λₜ::T # thermal conductivity
    cpᵥ::T # thermal volumetric capacity
end


@with_kw struct ModelParameters{T}
    geom::FlowCellGeometry2D{T}
    el_neg::ElectrodeParameters{T}
    el_pos::ElectrodeParameters{T}
    cc_neg::CurrentCollectorParameters{T}
    cc_pos::CurrentCollectorParameters{T}
    sep::SeparatorParameters{T}
    boundary::BoundaryConditions{T}
    discr::DiscretizationParameters{T}
    study::StudyParameters{T}
    scales::CharacteristicScales{T}
    scaling_params::ScalingParameters{T}
    electrolyte::ElectrolyteParams{T}
end
