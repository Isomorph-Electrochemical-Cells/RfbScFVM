function run(input_file_path)

    params = read_simulation_parameters(input_file_path)

    (grid, subgrids) = create_grid_2d(params.geom, reduced_membrane_model=true)

    initial_data = (p_l = params.boundary.p_in_neg,
                    p_r = params.boundary.p_in_pos,
                    ϕₛ_l = params.boundary.ϕₛ_neg,
                    ϕₛ_r = params.boundary.ϕₛ_pos,
                    ϕₗ_l = (params.boundary.ϕₛ_neg+params.boundary.ϕₛ_pos)/2,
                    ϕₗ_r = (params.boundary.ϕₛ_neg+params.boundary.ϕₛ_pos)/2,
                    c_ox_neg_l = params.boundary.species_neg[row="ox_neg",col="concentration"],
                    c_ox_neg_r = params.boundary.species_pos[row="ox_neg",col="concentration"],
                    c_red_neg_l = params.boundary.species_neg[row="red_neg",col="concentration"],
                    c_red_neg_r = params.boundary.species_pos[row="red_neg",col="concentration"],
                    c_ox_pos_l = params.boundary.species_neg[row="ox_pos",col="concentration"],
                    c_ox_pos_r = params.boundary.species_pos[row="ox_pos",col="concentration"],
                    c_red_pos_l = params.boundary.species_neg[row="red_pos",col="concentration"],
                    c_red_pos_r = params.boundary.species_pos[row="red_pos",col="concentration"],
                    temp_l = params.boundary.temp_amb,
                    temp_r = params.boundary.temp_amb,
                    temp_i = params.boundary.temp_amb)

    physics_data = (geom=params.geom,
                    el_neg=params.el_neg,
                    el_pos=params.el_pos,
                    cc_neg=params.cc_neg,
                    cc_pos=params.cc_pos,
                    sep=params.sep,
                    boundary=params.boundary,
                    discr=params.discr,
                    study=params.study,
                    scales=params.scales,
                    scaling_params=params.scaling_params,
                    electrolyte=params.electrolyte,
                    pressure_boundary_type = :p_in_v_out)

    sys = system(grid, physics_data)
    println("\nTotal degrees of freedom: ", VoronoiFVM.num_dof(sys))

    inival = initial_condition(sys, grid, initial_data, physics_data.study.non_isothermal)

    solution = solve_steady_state_problem(sys, inival)

    ####
    # Polarization study

    # Evaluate equilibrium potential # TODO: Move this to a separate function
    c_ox_neg_l = params.boundary.species_neg[row="ox_neg", col="concentration"]
    c_red_neg_l = params.boundary.species_neg[row="red_neg", col="concentration"]
    Δϕₛ_eq_neg = ϕ_eq(params.el_neg.reactions,
                      c_ox_neg_l, c_red_neg_l,
                      params.boundary.temp_amb)
    c_ox_pos_r = params.boundary.species_pos[row="ox_pos", col="concentration"]
    c_red_pos_r = params.boundary.species_pos[row="red_pos", col="concentration"]
    Δϕₛ_eq_pos = ϕ_eq(params.el_pos.reactions,
                      c_ox_pos_r, c_red_pos_r,
                      params.boundary.temp_amb)
    # Membrane Potential = Total potential difference due to Donnan Potentials
    z_c = params.electrolyte.species[row="counter", col="charge"]
    #TODO: Generalize to arbitrary electrolytes
    c_counter_l = (2*c_ox_neg_l + c_red_neg_l) / z_c
    c_counter_r = (2*c_ox_pos_r + c_red_pos_r) / z_c
    Δϕₘ_eq =  1/z_c * log(c_counter_l / c_counter_r)
    Δϕₛ_eq = Δϕₛ_eq_pos - Δϕₛ_eq_neg + Δϕₘ_eq

    params_polarization = params.study.polarization
    voltage_start = params_polarization.voltage_start
    voltage_stop = params_polarization.voltage_stop
    voltag_step = params_polarization.voltage_step
    if voltage_start === NaN
        voltage_start = Δϕₛ_eq
    end

    fontsize_theme = Theme(fontsize = 22)
    set_theme!(fontsize_theme)

    Δϕₛ_values = [Δϕₛ_value for Δϕₛ_value in range(voltage_start, voltage_stop, step=voltag_step)]
    println("Running polarization study...")
    println("Equilibrium voltage: ", Δϕₛ_eq, " [-], ",
                                      ustrip(Δϕₛ_eq * params.scales.V0), " [V]")
    num_Δϕₛ_values = length(Δϕₛ_values)
    integrated_reaction_terms = zeros(num_Δϕₛ_values, 2)
    solution = inival
    for idx_iter in eachindex(Δϕₛ_values)
        Δϕₛ_value = Δϕₛ_values[idx_iter]
        Δϕₛ_value_dim = ustrip(Δϕₛ_value * params.scales.V0)
        println("Applied voltage: ", Δϕₛ_value, " [-], ", Δϕₛ_value_dim, " [V]")
        solution_old = solution
        physics_data_updated = update_physics_data(physics_data, Δϕₛ_value)

        sys = system(grid, physics_data_updated)
        solution = solve_steady_state_problem(sys, solution_old)
        str_Δϕₛ_value = @sprintf "%0.4f_V" Δϕₛ_value_dim
        # plot_all_fields_2d(solution, grid, subgrids, physics_data_updated,
        #                    joinpath("polarization", "voltage_" * str_Δϕₛ_value))

        integrated_reaction_terms[idx_iter,:] .= integrate_reaction_terms(sys, solution)

        @assert isapprox(integrated_reaction_terms[idx_iter,1],
                        -integrated_reaction_terms[idx_iter,2], atol=1e-6)
    end

    current_density_values_nd = integrated_reaction_terms[:,2]/(params.geom.ly_cell*params.scaling_params.ϵL0)
    current_density_values = ustrip.(uconvert.(u"mA/cm^2", current_density_values_nd*params.scales.i0))
    voltage_values = ustrip.(Δϕₛ_values * params.scales.V0)


    output_path = params.study.polarization.output_folder
    mkpath(output_path)

    polarization_data = hcat(current_density_values, voltage_values)
    CSV.write(joinpath(output_path, "polarization.csv"), Tables.table(polarization_data);
              header=["CurrentDensity[mA/cm^2]", "CellVoltage[V]"])

    fig_polarization = plot_polarization_curve(current_density_values, voltage_values)
    output_path_figures = joinpath(pwd(), "output/figures")
    mkpath(output_path_figures)
    Makie.save(joinpath(output_path_figures, "polarization.pdf"), fig_polarization)

end
