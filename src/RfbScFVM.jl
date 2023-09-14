module RfbScFVM

using AxisKeys
import Base.show
using CSV
using ExtendableGrids
using FastGaussQuadrature
using FillArrays
using GridVisualize
using JSON
using LaTeXStrings
using LinearAlgebra
using Logging
using CairoMakie
using NonlinearSolve
using PrecompileTools
using Printf
using SparseArrays
using StaticArrays
using Tables
using VoronoiFVM
using Unitful

include("constants.jl")
include("units.jl")
include("utils.jl")
include("geometry.jl")
include("grid.jl")
include("parameters.jl")
include("kinetics.jl")
include("input.jl")
include("flux.jl")
include("subgrid_scale_models.jl")
include("physics_common.jl")
include("physics_1d.jl")
include("physics_2d.jl")
include("system.jl")
include("solver.jl")
include("plots.jl")
include("output.jl")
include("run.jl")


export read_config_file
export preprocess_parameters
export plot_grid
export create_grid_1d
export create_grid_2d
export save_results
export run_ocv
export run_polarization
export postprocess_results
export unitful_to_dict
export dimensional
export project_solutions
export plot_all_fields_1d
export plot_all_fields_2d
export plot_polarization_curve
export system

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)

    end
end


end # module RfbScFVM
