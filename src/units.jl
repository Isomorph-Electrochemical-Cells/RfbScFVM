# Fundamental units:

# 𝐈: Current
# 𝐋: Length
# 𝐌: Mass
# 𝐓: Time
# 𝚯: Temperature
# 𝐍: Amount

# Define additional derived units


@derived_dimension Permeability  Unitful.𝐋^2    true

@derived_dimension Diffusivity  Unitful.𝐋^2*Unitful.𝐓^-1    true

@derived_dimension CurrentDensity  Unitful.𝐈*Unitful.𝐋^-2    true

@derived_dimension VolumetricCurrentDensity  Unitful.𝐈*Unitful.𝐋^-3    true

@derived_dimension VolumetricHeatCapacity  Unitful.𝐌*Unitful.𝐋^-1*Unitful.𝐓^-2*Unitful.𝚯^-1   true

@derived_dimension SpecificHeatCapacity  Unitful.𝐋^2*Unitful.𝐓^-2*Unitful.𝚯^-1   true

@derived_dimension HeatConductivity  Unitful.𝐌*Unitful.𝐋^1*Unitful.𝐓^-3*Unitful.𝚯^-1   true

@derived_dimension HeterogeneousReactionRate Unitful.𝐋*Unitful.𝐓^-1   true

@derived_dimension MolarEntropy Unitful.𝐌*Unitful.𝐋^2*Unitful.𝐓^-2*Unitful.𝚯^-1*Unitful.𝐍^-1   true

@derived_dimension ArealCapacitance Unitful.𝐌^(-1)*Unitful.𝐋^(-4)*Unitful.𝐓^4*Unitful.𝐈^2  true

@derived_dimension VolumetricCapacitance Unitful.𝐌^(-1)*Unitful.𝐋^(-5)*Unitful.𝐓^4*Unitful.𝐈^2   true

@derived_dimension MolarMass Unitful.𝐌*Unitful.𝐍^-1 true
