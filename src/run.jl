function initial_condition(data)
    domain_sym = [:cc_neg, :el_neg, :sep, :el_pos, :cc_pos]
    num_species = length(data.electrolyte.species[col=1])

    Δϕₛ_eq_neg = ϕ_eq(data.el_neg.reactions, idx -> data.boundary.species_neg[idx],
                        data.boundary.temp_amb)
    Δϕₛ_eq_pos = ϕ_eq(data.el_pos.reactions, idx -> data.boundary.species_pos[idx],
                        data.boundary.temp_amb)

    ϕₗ_neg = data.boundary.ϕₛ_neg[] - Δϕₛ_eq_neg
    ϕₗ_pos = data.boundary.ϕₛ_pos[] - Δϕₛ_eq_pos

    vec_c_zero = SVector{num_species, RealType}(zeros(RealType, num_species))
    cc_neg_init = DomainVariables(p=zero(RealType),
                                    ϕₛ=data.boundary.ϕₛ_neg[],
                                    ϕₗ=zero(RealType),
                                    c=vec_c_zero,
                                    temp=data.boundary.temp_amb)
    cc_pos_init = DomainVariables(p=zero(RealType),
                                    ϕₛ=data.boundary.ϕₛ_pos[],
                                    ϕₗ=zero(RealType),
                                    c=zeros(RealType, num_species),
                                    temp=data.boundary.temp_amb)
    @assert length(data.electrolyte.species[col=1])-1 == length(data.boundary.species_neg)
    c_el_neg = copy(data.boundary.species_neg)
    p_neg = data.boundary.p_in_neg
    if occursin("1d", data.discr.spatial_discr)
        p_neg -= data.el_neg.μ * data.geom.ly_cell / (2 * data.el_neg.kₕ) * data.boundary.v_out_neg
    end
    el_neg_init = DomainVariables(p=p_neg,
                                    ϕₛ=data.boundary.ϕₛ_neg[],
                                    ϕₗ=ϕₗ_neg,
                                    c=c_el_neg,
                                    temp=data.boundary.temp_amb)
    @assert length(data.electrolyte.species[col=1])-1 == length(data.boundary.species_pos)
    c_el_pos = copy(data.boundary.species_pos)
    p_pos = data.boundary.p_in_pos
    if occursin("1d", data.discr.spatial_discr)
        p_pos -= data.el_pos.μ * data.geom.ly_cell / (2 * data.el_pos.kₕ) * data.boundary.v_out_pos
    end
    el_pos_init = DomainVariables(p=p_pos,
                                    ϕₛ=data.boundary.ϕₛ_pos[],
                                    ϕₗ=ϕₗ_pos,
                                    c=c_el_pos,
                                    temp=data.boundary.temp_amb)
    sep_init = DomainVariables(p=zero(RealType),
                                ϕₛ=zero(RealType),
                                ϕₗ=zero(RealType),
                                c=vec_c_zero,
                                temp=data.boundary.temp_amb)

    initial_data = (cc_neg_init, el_neg_init, el_pos_init, cc_pos_init, sep_init)

    return initial_data
end

function run_ocv(data)
    Δϕₛ_eq_neg = 0.0
    for reaction in data.el_neg.reactions
        Δϕₛ_eq_neg += ϕ_eq(reaction, idx -> data.boundary.species_neg[idx],
                           data.boundary.temp_amb)
    end
    Δϕₛ_eq_pos = 0.0
    for reaction in data.el_pos.reactions
        Δϕₛ_eq_pos += ϕ_eq(reaction, idx -> data.boundary.species_pos[idx],
                           data.boundary.temp_amb)
    end

    cc_pos = cc_neg = zero(RealType)
    species = data.electrolyte.species
    num_species = size(species)[1] - 1
    species_neg = data.boundary.species_neg
    for idx in 1:num_species
        cc_neg -= species[row=idx, col=Key("charge")] * species_neg[idx]
    end
    species_pos = data.boundary.species_pos
    for idx in 1:num_species
        cc_pos -= species[row=idx, col=Key("charge")] * species_pos[idx]
    end
    z_c = RealType(data.electrolyte.species[row=end, col=Key("charge")])
    cc_neg /= z_c
    cc_pos /= z_c

    z = data.electrolyte.species[col=Key("charge")].data.data
    cf = data.sep.cf
    zf = data.sep.zf
    temp = data.boundary.temp_amb

    c_neg = vcat(species_neg, cc_neg)
    ϕₗ_l = 0.0
    ϕₗ_sep_l = donnan_equlibrium_potential(ϕₗ_l, c_neg, z, cf, zf, temp)

    c_pos = vcat(species_pos, cc_pos)
    ϕₗ_r = 0.0
    ϕₗ_sep_r = donnan_equlibrium_potential(ϕₗ_r, c_pos, z, cf, zf, temp)
    ###########################

    Δϕₛ_eq = Δϕₛ_eq_pos - Δϕₛ_eq_neg + ϕₗ_sep_l - ϕₗ_sep_r

    return OcvResults(Δϕₛ_eq)
end

function run_polarization(sys, grid, subgrids, data, Δϕₛ_eq; save_figures=false)
    params_polarization = data.study.polarization
    voltage_start = RealType(params_polarization.voltage_start)
    voltage_start = RealType(voltage_start === NaN ? Δϕₛ_eq : voltage_start)
    voltage_stop = RealType(params_polarization.voltage_stop)
    voltage_stop = RealType(voltage_stop === NaN ? Δϕₛ_eq : voltage_stop)
    voltage_step = RealType(params_polarization.voltage_step)

    data.boundary.ϕₛ_neg[] = -voltage_start/2
    data.boundary.ϕₛ_pos[] = voltage_start/2

    initial_data = initial_condition(data)

    vec_Δϕₛ = [Δϕₛ_value for Δϕₛ_value in
                  range(voltage_start, voltage_stop, step=voltage_step)]
    println("Running polarization study...")
    println("Equilibrium voltage: ", ustrip(Δϕₛ_eq * data.scales.V0), " [V]")

    num_Δϕₛ_values = length(vec_Δϕₛ)
    init = get_initial_condition(sys, grid, initial_data, data)

    indices = unknown_indices(init)
    data = ModelParameters(; data..., indices=indices)
    if dim_space(grid) == 1
        physics!(sys, physics_1d(data))
    elseif dim_space(grid) == 2
        physics!(sys, physics_2d(data))
    end

    if data.study.generate_figures
        if occursin("1d", data.discr.spatial_discr)
            dict_figures = plot_all_fields_1d(init, grid, subgrids, data)
        else
            dict_figures = plot_all_fields_2d(init, grid, subgrids, data)
        end

        if save_figures
            relative_folder_path = data.study.output_folder
            relative_path = joinpath(relative_folder_path, "polarization", "initial_condition")
            save_plots(dict_figures; relative_path=relative_path)
        end
    end

    solution = solve_steady_state_problem(sys, init, vec_Δϕₛ)

    @assert num_Δϕₛ_values == length(solution.t)

    vec_i = sys.physics.data.results.current_density

    return PolarizationResults(vec_Δϕₛ, vec_i, -vec_Δϕₛ .* vec_i, solution.u)
end

function save_plots(dict_figures::Dict{String, Makie.Figure};
                    relative_path, file_format="pdf", file_name_suffix="")
    dir_path = dirname(dirname(pathof(@__MODULE__)))
    output_dir_path = joinpath(dir_path, relative_path)
    mkpath(output_dir_path)

    for (name, fig) in dict_figures
        file_name = name * file_name_suffix * "." * file_format
        Makie.save(joinpath(output_dir_path, file_name), fig)
    end
end

function run(input_file_path)

    dict_params = read_config_file(input_file_path)
    data = preprocess_parameters(dict_params)

    if data.discr.spatial_discr == "fvm_2d"
        (grid, subgrids) = create_grid_2d(data.geom, data.mesh, data.dom, data.bnd)
    elseif data.discr.spatial_discr == "fvm_1d"
        (grid, subgrids) = create_grid_1d(data.geom, data.mesh, data.dom, data.bnd)
    else
        @assert false "Unknown spatial discretization scheme"
    end

    if data.study.generate_figures
        output_path = data.study.output_folder
        mkpath(output_path)
        output_path_figures = joinpath(pwd(), output_path)
        vis = plot_grid(grid; plotter=CairoMakie);
        save(joinpath(output_path_figures, "grid.png"), vis)
    end

    results = Dict()

    sys = system(grid, data)
    max_node_id = maximum(grid.components[ExtendableGrids.CellNodes])
    min_node_id = minimum(grid.components[ExtendableGrids.CellNodes])
    num_nodes = max_node_id - min_node_id + 1
    @info "\nNumber of grid nodes: " * string(num_nodes)
    @info "\nTotal degrees of freedom: " * string(VoronoiFVM.num_dof(sys))

    # run simulations
    ocv_results = run_ocv(data)
    polarization_results = run_polarization(sys, grid, subgrids, data, ocv_results.Δϕₛ;
                                             save_figures=true)

    # post-process results
    dict_ocv_results = postprocess_results(ocv_results, data)
    dict_polarization_results = postprocess_results(polarization_results, data)

    # save results
    dict_results = Dict()
    merge!(dict_results, dict_ocv_results)
    merge!(dict_results, dict_polarization_results)
    filepath = joinpath(pwd(), data.study.output_folder, data.study.output_file_name)
    save_results(dict_results, filepath)

    if data.study.generate_figures
        fig_polarization = plot_polarization_curve(dict_polarization_results["polarization"])
        output_path_figures = joinpath(pwd(), data.study.output_folder)
        mkpath(output_path_figures)
        Makie.save(joinpath(output_path_figures, "polarization.pdf"), fig_polarization)

        solution = polarization_results.vec_solution
        cell_values = dict_polarization_results["polarization"]["cell_voltage"]["value"]

        if occursin("2d", data.discr.spatial_discr)
            grid_proj_1d, subgrids_proj_1d = create_grid_1d(data.geom, Mesh1D(data.mesh), data.dom, data.bnd)
        end

        polarization_folder = joinpath(data.study.output_folder , "polarization")
        for idx in eachindex(solution)
            str_voltage = @sprintf "%0.3f" cell_values[idx]
            subfolder_name = "voltage_" * str_voltage * "_[V]_" * string(idx)
            if occursin("1d", data.discr.spatial_discr)
                dict_figures = plot_all_fields_1d(solution[idx], grid, subgrids, data)
                relative_path = joinpath(polarization_folder, subfolder_name)
                save_plots(dict_figures; relative_path=relative_path)
            else
                dict_figures_2d = plot_all_fields_2d(solution[idx], grid, subgrids, data)
                relative_path = joinpath(polarization_folder, subfolder_name)
                save_plots(dict_figures_2d; relative_path=relative_path,
                           file_name_suffix = "_2d")

                solution_proj = project_solutions(solution[idx], sys, grid, subgrids, data)
                dict_figures = plot_all_fields_1d(solution_proj[2], grid_proj_1d, subgrids_proj_1d, data)

                relative_path = joinpath(polarization_folder, "avg_x")
                save_plots(dict_figures; relative_path=relative_path)
            end
        end
    end
end
