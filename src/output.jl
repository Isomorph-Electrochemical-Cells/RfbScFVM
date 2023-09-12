abstract type Results end

struct OcvResults{T} <: Results
    Δϕₛ::T
end

function Base.show(io::IO, result::OcvResults{T}) where T
    print(io, "Δϕₛ = ", result.Δϕₛ)
end

function dimensional(result::OcvResults{T}, data) where T<:AbstractFloat
    return OcvResults(uconvert(u"V", result.Δϕₛ * data.scales.V0))
end

function postprocess_results(result::OcvResults{T}, data) where T
    dim_ocv_result = dimensional(result, data)
    dict_ocv_result = unitful_to_dict(dim_ocv_result.Δϕₛ)
    return Dict("ocv"=>dict_ocv_result)
end

function save_results(dict_results, file_path)
    indentation = 2 # whitespace characters used for indentation
    open(file_path,"w") do io
        JSON.print(io, dict_results, indentation)
    end
end

struct PolarizationResults{VoltageType, CurrentDensityType, PowerDensityType, SolutionArray} <: Results
    vec_Δϕₛ::Vector{VoltageType} # cell voltage values
    vec_i::Vector{CurrentDensityType} # current density values
    vec_p::Vector{PowerDensityType} # power density values
    vec_solution::Vector{SolutionArray} # numerical solutions
end

function Base.show(io::IO, result::PolarizationResults)
    zipped_result = collect(zip(result.vec_Δϕₛ, result.vec_i, result.vec_p))
    print(io, "(cell voltage, current density, power density) = ", zipped_result)
end

function postprocess_results(result::PolarizationResults, data)
    dim_result = dimensional(result, data)

    dict_dim_vec_Δϕₛ = unitful_to_dict(dim_result.vec_Δϕₛ)
    dict_dim_vec_i = unitful_to_dict(dim_result.vec_i)
    dict_dim_vec_p = unitful_to_dict(dim_result.vec_p)

    dict_polarization = Dict("cell_voltage"=>dict_dim_vec_Δϕₛ,
                             "current_density"=>dict_dim_vec_i,
                             "power_density"=>dict_dim_vec_p)
    # dict_polarization = Dict("Δϕₛ"=>dict_dim_vec_Δϕₛ, "i"=>dict_dim_vec_i,
    #                          "solutions"=>dim_result.vec_solution)
    return Dict("polarization"=>dict_polarization)
end

function dimensional(result::PolarizationResults, data)
    return PolarizationResults(uconvert.(u"V", result.vec_Δϕₛ * data.scales.V0),
                               uconvert.(u"mA/cm^2", result.vec_i * data.scales.i0),
                               uconvert.(u"mW/cm^2", result.vec_p * data.scales.i0 * data.scales.V0),
                               result.vec_solution)
end
