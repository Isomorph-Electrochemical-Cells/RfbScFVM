module RfbScFVM

using AxisArrays
using Printf
using JSON
using FillArrays
using Parameters
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LaTeXStrings
using Unitful
using DataFrames
using Tables
using CSV
using Makie, CairoMakie, GLMakie
using Infiltrator

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
include("physics.jl")
include("system.jl")
include("solver.jl")
include("study.jl")
include("postprocessing.jl")
include("run.jl")

export flow_cell_geometry_2D
export create_grid_2d
export physics
export system
export initial_condition
export cell_thickness
export solve_steady_state_problem, solve_transient_problem
export plot_all_fields_2d, plot_all_fields_1d
export read_simulation_parameters
export project_solution, project_solutions, combined_solution_fields_1d
export ReactionParameters
export update_physics_data
export integrate_reaction_terms
export plot_polarization_curve
export volumetric_current_density
export ϕ_eq

end # module RfbScFVM
