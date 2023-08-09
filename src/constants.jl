# Define the float type to be used throughout the model
const RealType = Float64
const IndexType = Int64

# Define physical constants
const R = RealType(8.31446261815324)u"J/(K*mol)" # Ideal gas constant [J/(K mol)]
const F = RealType(96485.33212)u"C/mol" # Faraday constant [C/mol]*)
const NA = RealType(6.02214076e+23)u"1/mol" # Avogadro constant 1/Mol
const EV = RealType(1.602176634e-19)u"J" # Electron volt [J]
const KB = RealType(1.380649e-23)u"J/K" # Boltzmann constant [J/K]
const EC0 = (F/NA) # Electron charge [C]

# @enum Variables begin
#     p_l = 1
#     p_r = 2
#     ϕₛ_l = 3
#     ϕₛ_r = 4
#     ϕₗ_l = 5
#     ϕₗ_r = 6
#     c_ox_neg_l = 7
#     c_ox_neg_r = 8
#     c_red_neg_l = 9
#     c_red_neg_r = 10
#     c_ox_pos_l = 11
#     c_ox_pos_r = 12
#     c_red_pos_l = 13
#     c_red_pos_r = 14
#     temp_l = 15
#     temp_r = 16
#     temp_i = 17
# end

# @enum Domains begin
#     cc_neg = 1
#     el_neg = 2
#     el_pos = 3
#     cc_pos = 4
#     sep = 5
# end

# @enum Boundaries begin
#     cc_neg_left = 6
#     cc_neg_el_neg = 7
#     el_pos_cc_pos = 8
#     cc_pos_right = 9
#     el_neg_inflow = 10
#     el_neg_outflow = 11
#     el_pos_inflow = 12
#     el_pos_outflow = 13
# end
