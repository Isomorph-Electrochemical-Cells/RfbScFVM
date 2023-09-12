module TestModel1Dvs2D

using RfbScFVM
using Test


function run_system(dict_input_data)
    # preprocess parameters
    data = preprocess_parameters(dict_input_data)

    # greate grid and generate visualization
    if data.discr.spatial_discr == "fvm_2d"
        (grid, subgrids) = create_grid_2d(data.geom, data.mesh, data.dom, data.bnd)
    elseif data.discr.spatial_discr == "fvm_1d"
        (grid, subgrids) = create_grid_1d(data.geom, data.mesh, data.dom, data.bnd)
    end

    # vis = plot_grid(grid; plotter=CairoMakie);
    # output_path_figures = ""
    # save(joinpath(output_path_figures, "grid.png"), vis)

    ocv = run_ocv(data)

    sys = system(grid, data)
    polarization = run_polarization(sys, grid, subgrids, data, ocv.Δϕₛ)

    ocv = postprocess_results(ocv, data)
    polarization = postprocess_results(polarization, data)

    # dict_results = Dict()
    # merge!(dict_results, ocv)
    # merge!(dict_results, polarization)

    return polarization["polarization"]
end

function test_1d_vs_2d_systems(;non_isothermal=non_isothermal)

    # read configuration parameters from file
    input_file_path = "test/test_data/mv_temptma_polarization_soc50_isothermal_1d.json"

    dict_input_data = read_config_file(input_file_path)
    dict_input_data["study_parameters"]["non_isothermal"] = non_isothermal

    dict_input_data["discretization_parameters"]["spatial_discretization"] = "fvm_1d"
    results_1d = run_system(dict_input_data)

    dict_input_data["discretization_parameters"]["spatial_discretization"] = "fvm_2d"
    results_2d = run_system(dict_input_data)

    @test results_1d["cell_voltage"]["value"] ≈ results_2d["cell_voltage"]["value"] rtol=1e-6 atol=1e-6
    @test results_1d["current_density"]["value"] ≈ results_2d["current_density"]["value"] rtol=1e-2 atol=1e-2
end

function test()
    test_1d_vs_2d_systems(;non_isothermal=false)
    test_1d_vs_2d_systems(;non_isothermal=true)
end

end
