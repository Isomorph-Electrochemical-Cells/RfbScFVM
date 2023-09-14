module TestHomogeneousReactions

using RfbScFVM
using Test
using Unitful


function run_system(dict_input_data)
    # preprocess parameters
    data = preprocess_parameters(dict_input_data)

    # greate grid and generate visualization
    if data.discr.spatial_discr == "fvm_2d"
        (grid, subgrids) = create_grid_2d(data.geom, data.mesh, data.dom, data.bnd)
    elseif data.discr.spatial_discr == "fvm_1d"
        (grid, subgrids) = create_grid_1d(data.geom, data.mesh, data.dom, data.bnd)
    end

    ocv = run_ocv(data)

    sys = system(grid, data)
    polarization = run_polarization(sys, grid, subgrids, data, ocv.Δϕₛ)

    return (polarization.vec_solution[1], data, subgrids)
end

function homogeneous_reaction_1d()
    # read configuration parameters from file

    package_path = pkgdir(RfbScFVM)
    input_file_path = "test/test_data/homogeneous_reactions.json"

    dict_input_data = read_config_file(joinpath(package_path, input_file_path))
    dict_input_data["discretization_parameters"]["spatial_discretization"] = "fvm_1d"
    (solution, data, subgrids) = run_system(dict_input_data)

    species_names = ["A", "B"]
    # numerical reference solution evaluated with Mathematica
    ref_solutions = [0.875394u"mol/l", 2.12461u"mol/l"]

    max_rel_err = 0.0
    for idx_species_name in eachindex(species_names)
        species_name = species_names[idx_species_name]
        ref_sol = ref_solutions[idx_species_name]

        species_idx = findfirst(x->x==species_name, data.electrolyte.species.row)

        idx_c_sol = data.idx[data.dom[:el_neg]].c[species_idx]
        c_sol = view(solution[idx_c_sol, :], subgrids.subgrid_el_neg) * data.scales.C0
        max_rel_err = max(maximum(abs.((c_sol .- ref_sol)/ref_sol)), max_rel_err)
    end
    @test max_rel_err < 0.01

end

function test()
    homogeneous_reaction_1d()
end

end
