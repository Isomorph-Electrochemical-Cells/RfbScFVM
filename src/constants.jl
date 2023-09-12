# Define the float type to be used throughout the model
const RealType = Float64
const IntType = Int64
const SQRT_EPS = sqrt(eps(RealType))

# Define physical constants
const R = RealType(8.31446261815324)u"J/(K*mol)" # Ideal gas constant [J/(K mol)]
const F = RealType(96485.33212)u"C/mol" # Faraday constant [C/mol]*)
const NA = RealType(6.02214076e+23)u"1/mol" # Avogadro constant 1/Mol
const EV = RealType(1.602176634e-19)u"J" # Electron volt [J]
const KB = RealType(1.380649e-23)u"J/K" # Boltzmann constant [J/K]
const EC0 = (F/NA) # Electron charge [C]
