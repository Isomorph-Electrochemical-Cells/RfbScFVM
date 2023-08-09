function initial_condition(prms)
    domain_sym = [:cc_neg, :el_neg, :sep, :el_pos, :cc_pos]
    num_species = length(prms.electrolyte.species[col=1])

    ϕₗ_neg = prms.boundary.ϕₛ_neg[] - prms.el_neg.reactions.Δϕ₀
    ϕₗ_pos = prms.boundary.ϕₛ_pos[] - prms.el_pos.reactions.Δϕ₀

    vec_c_zero = SVector{num_species, RealType}(zeros(RealType, num_species))
    cc_neg_init = DomainVariables(p=zero(RealType),
                                    ϕₛ=prms.boundary.ϕₛ_neg[],
                                    ϕₗ=zero(RealType),
                                    c=vec_c_zero,
                                    temp=prms.boundary.temp_amb)
    cc_pos_init = DomainVariables(p=zero(RealType),
                                    ϕₛ=prms.boundary.ϕₛ_pos[],
                                    ϕₗ=zero(RealType),
                                    c=zeros(RealType, num_species),
                                    temp=prms.boundary.temp_amb)
    c_el_neg = copy(prms.boundary.species_neg)
    el_neg_init = DomainVariables(p=prms.boundary.p_in_neg,
                                    ϕₛ=prms.boundary.ϕₛ_neg[],
                                    ϕₗ=ϕₗ_neg,
                                    c=c_el_neg,
                                    temp=prms.boundary.temp_amb)
    c_el_pos = copy(prms.boundary.species_pos)
    el_pos_init = DomainVariables(p=prms.boundary.p_in_pos,
                                    ϕₛ=prms.boundary.ϕₛ_pos[],
                                    ϕₗ=ϕₗ_pos,
                                    c=c_el_pos,
                                    temp=prms.boundary.temp_amb)
    sep_init = DomainVariables(p=zero(RealType),
                                ϕₛ=zero(RealType),
                                ϕₗ=zero(RealType),
                                c=vec_c_zero,
                                temp=prms.boundary.temp_amb)

    initial_data = (cc_neg_init, el_neg_init, el_pos_init, cc_pos_init, sep_init)

    return initial_data
end

function run_ocv(prms)
    idx_c2u = collect(1:length(prms.boundary.species_neg))
    Δϕₛ_eq_neg = ϕ_eq(prms.el_neg.reactions,
                      prms.boundary.species_neg, idx_c2u,
                      prms.boundary.temp_amb)
    Δϕₛ_eq_pos = ϕ_eq(prms.el_pos.reactions,
                      prms.boundary.species_pos, idx_c2u,
                      prms.boundary.temp_amb)

    c_counter_neg = c_counter_pos = zero(RealType)
    species = prms.electrolyte.species
    num_species = size(species)[1]
    species_neg = prms.boundary.species_neg
    for idx in 1:num_species
        c_counter_neg -= species[row=idx, col=Key("charge")] * species_neg[idx]
    end
    species_pos = prms.boundary.species_pos
    for idx in 1:num_species
        c_counter_pos -= species[row=idx, col=Key("charge")] * species_pos[idx]
    end
    z_c = RealType(prms.electrolyte.counter(col="charge"))
    c_counter_neg /= z_c
    c_counter_pos /= z_c

    # Membrane Potential = Total potential difference due to Donnan Potentials
    Δϕₘ_eq =  1/z_c * log(c_counter_neg / c_counter_pos)
    Δϕₛ_eq = Δϕₛ_eq_pos - Δϕₛ_eq_neg + Δϕₘ_eq

    return Δϕₛ_eq
end

function run_polarization(sys, grid, prms, Δϕₛ_eq)

    params_polarization = prms.study.polarization
    voltage_start = params_polarization.voltage_start
    voltage_start = voltage_start === NaN ? Δϕₛ_eq : voltage_start

    voltage_stop = params_polarization.voltage_stop
    voltage_stop = voltage_stop === NaN ? Δϕₛ_eq : voltage_stop

    voltag_step = params_polarization.voltage_step

    prms.boundary.ϕₛ_neg[] = -voltage_start/2
    prms.boundary.ϕₛ_pos[] = voltage_start/2

    initial_data = initial_condition(prms)

    Δϕₛ_values = [Δϕₛ_value for Δϕₛ_value in
                  range(voltage_start, voltage_stop, step=voltag_step)]
    println("Running polarization study...")
    println("Equilibrium voltage: ", Δϕₛ_eq, " [-], ",
                                      ustrip(Δϕₛ_eq * prms.scales.V0), " [V]")

    num_Δϕₛ_values = length(Δϕₛ_values)
    integrated_reaction_terms = zeros(num_Δϕₛ_values, 2)
    init = initial_condition(sys, grid, initial_data, prms)
    solution = solve_steady_state_problem(sys, init, Δϕₛ_values)

    @assert num_Δϕₛ_values == length(solution.t)

    return (Δϕₛ_values, solution, sys)
end

function run(input_file_path)

    prms = read_simulation_parameters(input_file_path)

    (grid, subgrids) = create_grid_2d(prms.geom, prms.mesh, prms.dom, prms.bnd)
    results = Dict()

    sys = system(grid, prms)
    println("\nTotal degrees of freedom: ", VoronoiFVM.num_dof(sys))

    Δϕₛ_eq = run_ocv(prms)
    results["ocv"] = unitful_to_dict(Δϕₛ_eq * prms.scales.V0, "V")

    (Δϕₛ_values, solution, sys) = run_polarization(sys, grid, prms, Δϕₛ_eq)

    current_density_values_nd = sys.physics.data.results.current_density
    current_density_values = ustrip.(uconvert.(u"mA/cm^2", current_density_values_nd*prms.scales.i0))
    voltage_values = ustrip.(Δϕₛ_values * prms.scales.V0)
    power_density_values = abs.(current_density_values .* voltage_values)

    dict_voltage = Dict("value"=>voltage_values, "unit"=>"V")
    dict_current_density = Dict("value"=>current_density_values, "unit"=>"mA/cm^2")
    dict_power_density = Dict("value"=>power_density_values, "unit"=>"mW/cm^2")

    results["polarization"] = Dict("voltage"=>dict_voltage,
                                   "current_density"=>dict_current_density,
                                   "power_density"=>dict_power_density)

    output_path = prms.study.output_folder
    mkpath(output_path)

    filepath = joinpath(pwd(), prms.study.output_folder, prms.study.output_file_name)
    output_results(results, filepath)

    polarization_data = hcat(current_density_values, voltage_values)
    CSV.write(joinpath(output_path, "polarization.csv"), Tables.table(polarization_data);
              header=["CurrentDensity[mA/cm^2]", "CellVoltage[V]"])

    if prms.study.generate_figures
        fig_polarization = plot_polarization_curve(current_density_values, voltage_values)
        output_path_figures = joinpath(pwd(), prms.study.output_folder)
        mkpath(output_path_figures)
        Makie.save(joinpath(output_path_figures, "polarization.pdf"), fig_polarization)

        for idx in eachindex(solution)
            str_voltage = @sprintf "%0.3f" voltage_values[idx]
            filename = "voltage_" * str_voltage * "_[V]_" * string(idx) * ".pdf"
            plot_all_fields_2d(solution[idx], grid, subgrids, prms,
                            joinpath("polarization", filename))
        end
    end
end
