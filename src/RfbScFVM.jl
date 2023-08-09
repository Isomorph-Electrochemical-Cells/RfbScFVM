module RfbScFVM

using AxisKeys
using CSV
using DataFrames
using ExtendableGrids
using FillArrays
using GridVisualize
using JSON
using LaTeXStrings
using Makie, CairoMakie, GLMakie
using PrecompileTools
using Printf
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
include("study.jl")
include("postprocessing.jl")
include("output.jl")
include("run.jl")


export flow_cell_geometry
export create_grid_2d
export system
export initial_condition
export cell_thickness
export solve_steady_state_problem, solve_transient_problem
export plot_all_fields_2d, plot_all_fields_1d
export read_simulation_parameters
export ReactionParameters
export integrate_reaction_terms
export plot_polarization_curve
export volumetric_current_density
export ϕ_eq
export gridplot

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    # TODO
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        # TODO
    end
end


end # module RfbScFVM
