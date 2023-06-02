function parse_json_file(input_file_path)
    dict_input_parameters = JSON.parsefile(input_file_path; dicttype=Dict, inttype=Int64,
                                           use_mmap=true)
    return dict_input_parameters
end

function read_simulation_parameters(input_file_path)
    dict_params = parse_json_file(input_file_path)
    return preprocess_model_parameters(dict_params)
end

function dict_to_unitful(dict_quantity)
    RealType(dict_quantity["value"]) * uparse(dict_quantity["unit"])
end

function preprocess_model_parameters(dict_params)
    dict_model_params = dict_params["model_parameters"]

    # geometry
    dict_geom = dict_model_params["geometry"]

    lx_cc_neg = dict_to_unitful(dict_geom["current_collector_neg_thickness"])
    lx_el_neg = dict_to_unitful(dict_geom["electrode_neg_thickness"])
    lx_sep = dict_to_unitful(dict_geom["separator_thickness"])
    lx_el_pos = dict_to_unitful(dict_geom["electrode_pos_thickness"])
    lx_cc_pos = dict_to_unitful(dict_geom["current_collector_pos_thickness"])
    ly_cell = dict_to_unitful(dict_geom["cell_height"])

    lx = lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos + lx_cc_pos
    lx_cc_neg = lx_cc_neg / lx |> NoUnits
    lx_el_neg = lx_el_neg / lx |> NoUnits
    lx_sep = lx_sep / lx |> NoUnits
    lx_el_pos = lx_el_pos / lx |> NoUnits
    lx_cc_pos = lx_cc_pos / lx |> NoUnits
    ly_cell = ly_cell / lx |> NoUnits

    mesh_params = dict_params["discretization_parameters"]["mesh"]
    @assert mesh_params["type"] == "rectilinear"
    mesh_sizes = mesh_params["relative_mesh_sizes"]
    hy_rel_cell = mesh_sizes["hy_cell"]
    hx_rel_cc_neg = mesh_sizes["hx_cc_neg"]
    hx_rel_el_neg = mesh_sizes["hx_el_neg"]
    hx_rel_el_pos = mesh_sizes["hx_el_pos"]
    hx_rel_cc_pos = mesh_sizes["hx_cc_pos"]

    geom = flow_cell_geometry_2D(lx_cc_neg=lx_cc_neg, lx_el_neg=lx_el_neg,
        lx_sep=lx_sep, lx_el_pos=lx_el_pos, lx_cc_pos=lx_cc_pos, ly_cell=ly_cell,
        hy_cell=hy_rel_cell*ly_cell,
        hx_cc_neg=hx_rel_cc_neg*lx_cc_neg, hx_el_neg=hx_rel_el_neg*lx_el_neg,
        hx_cc_pos=hx_rel_cc_pos*lx_cc_pos, hx_el_pos=hx_rel_el_pos*lx_el_pos)

    dict_boundary_conditions = dict_params["boundary_conditions"]
    temp_amb = dict_to_unitful(dict_boundary_conditions["ambient_temperature"])
    # define characteristic scales
    scales = CharacteristicScales{RealType}(L0=lx, TEMP0=temp_amb)

    # negative and positive current collectors
    dict_cc_neg_pos = [dict_model_params["current_collector_neg"],
                       dict_model_params["current_collector_pos"]]

    cc_neg_pos = Array{CurrentCollectorParameters{RealType}}(undef,2)
    for (index, dict_cc) in enumerate(dict_cc_neg_pos)
        σₑ = dict_to_unitful(dict_cc["electric_conductivity"])/scales.σ0 |> NoUnits
        hₜ = dict_to_unitful(dict_cc["heat_transfer_coefficient"])/(scales.λ0/scales.L0) |> NoUnits
        λₜ = dict_to_unitful(dict_cc["thermal_conductivity"])/scales.λ0 |> NoUnits
        cpᵥ = dict_to_unitful(dict_cc["thermal_vol_capacity"])/scales.CPv0 |> NoUnits

        cc_neg_pos[index] = CurrentCollectorParameters{RealType}(σₑ=σₑ, hₜ=hₜ, λₜ=λₜ, cpᵥ=cpᵥ)
    end

    # # negative and postive electrodes
    dict_el_neg_pos = [dict_model_params["electrode_neg"],
                       dict_model_params["electrode_pos"]]

    el_neg_pos = Array{ElectrodeParameters{RealType}}(undef,2)
    for (index, dict_el) in enumerate(dict_el_neg_pos)
        λₜ = dict_to_unitful(dict_el["thermal_conductivity"])/scales.λ0 |> NoUnits
        cpᵥ = dict_to_unitful(dict_el["thermal_vol_capacity"])/scales.CPv0 |> NoUnits
        εₗ = dict_to_unitful(dict_el["porosity"]) |> NoUnits
        kₕ = dict_to_unitful(dict_el["hydraulic_permeability"])/scales.KH0 |> NoUnits
        dict_mass_transfer_model = dict_el["mass_transfer_model"]
        mass_transfer_factor = dict_to_unitful(dict_mass_transfer_model["mass_transfer_factor"])
        mass_transfer_exp = NoUnits(dict_to_unitful(dict_mass_transfer_model["mass_transfer_exponent"]))
        mass_transfer_exp = Rational{Int32}(mass_transfer_exp) # convert float to rational
        sh_factor = NoUnits(mass_transfer_factor*(scales.LP0/scales.D0)*scales.VEL0^mass_transfer_exp)
        sh_exponent = mass_transfer_exp
        σₑ = dict_to_unitful(dict_el["electric_conductivity"])/scales.σ0 |> NoUnits
        aᵥ = dict_to_unitful(dict_el["specific_surface_area"])*scales.LP0 |> NoUnits
        μ = dict_to_unitful(dict_el["dynamic_viscosity"])/scales.μ0 |> NoUnits
        Cₐ = dict_to_unitful(dict_el["areal_capacitance"])/scales.Ca0 |> NoUnits
        Cᵥ = Cₐ * aᵥ # volumetric capacitance
        diffusivity_model = dict_el["effective_diffusivity_model"]
        Deff_const = dict_to_unitful(diffusivity_model["constant_coefficient"]) |> NoUnits
        Deff_linear = dict_to_unitful(diffusivity_model["linear_coefficient"])/(scales.L0/scales.D0) |> NoUnits
        Deff_quadratic = dict_to_unitful(diffusivity_model["quadratic_coefficient"])/(scales.L0^2/scales.D0^2) |> NoUnits

        reaction = dict_el["reactions"][1] # TODO: Support multiple reactions
        reaction_name = reaction["name"]
        reduction_potential = reaction["standard_reduction_potential"]
        Δϕ₀ = dict_to_unitful(reduction_potential["standard_value"])/scales.V0 |> NoUnits
        ∂Δϕ₀_∂temp = dict_to_unitful(reduction_potential["temperature_coefficient"])/(scales.V0/scales.TEMP0) |> NoUnits

        Δs = dict_to_unitful(reaction["entropy_change"])/scales.S0 |> NoUnits
        dict_kinetics = reaction["kinetics"]
        α = dict_to_unitful(dict_kinetics["transfer_coefficient"]) |> NoUnits
        ki = dict_to_unitful(dict_kinetics["rate_constant"])/scales.K0 |> NoUnits # kinetic number
        temp_ref = dict_to_unitful(reaction["reference_temperature"])/scales.TEMP0 |> NoUnits


        coeffs = reaction["stoichiometric_coefficients"]
        # TODO: Initialize stoichiometric coefficients from inputs
        ν_ox = -1
        ν_red = 1
        ν_el = -1

        reactions = ReactionParameters(name=reaction_name,
                                      Δϕ₀=Δϕ₀, ∂Δϕ₀_∂temp=∂Δϕ₀_∂temp,
                                      Δs=Δs, α=α, ki=ki, temp_ref=temp_ref,
                                      ν_ox=ν_ox, ν_red=ν_red, ν_el=ν_el)

        el_neg_pos[index] = ElectrodeParameters{RealType}(λₜ=λₜ, cpᵥ=cpᵥ, εₗ=εₗ, kₕ=kₕ,
                                                          μ=μ, σₑ=σₑ, aᵥ=aᵥ, Cᵥ=Cᵥ,
                                                          Deff_const=Deff_const,
                                                          Deff_linear=Deff_linear,
                                                          Deff_quadratic=Deff_quadratic,
                                                          reactions=reactions,
                                                          sh_factor = sh_factor,
                                                          sh_exponent = sh_exponent)
    end

    # # separator
    dict_sep = dict_model_params["separator"]
    λₜ = dict_to_unitful(dict_sep["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = dict_to_unitful(dict_sep["thermal_vol_capacity"])/scales.CPv0 |> NoUnits
    kh = dict_to_unitful(dict_sep["hydraulic_permeability"])/scales.KH0 |> NoUnits
    kϕ = dict_to_unitful(dict_sep["electrokinetic_permeability"])/scales.Kϕ0 |> NoUnits
    μ = dict_to_unitful(dict_sep["dynamic_viscosity"])/scales.μ0 |> NoUnits
    conductivity = dict_sep["electrical_conductivity"]
    σₑ = dict_to_unitful(conductivity["reference_value"])/scales.σ0 |> NoUnits
    ∂σₑ∂temp = dict_to_unitful(conductivity["temperature_coefficient"])/(scales.σ0/scales.TEMP0) |> NoUnits
    temp_ref = dict_to_unitful(conductivity["reference_temperature"])/(scales.TEMP0) |> NoUnits

    dict_fixed_ionic_groups = dict_sep["fixed_ionic_groups"]
    cf = dict_to_unitful(dict_fixed_ionic_groups["concentration"])/scales.C0 |> NoUnits
    zf = dict_to_unitful(dict_fixed_ionic_groups["charge"]) |> NoUnits

    sep = SeparatorParameters{RealType}(λₜ=λₜ, cpᵥ=cpᵥ, kh=kh, kϕ=kϕ, μ=μ,
                                        σₑ=σₑ, ∂σₑ∂temp=∂σₑ∂temp, temp_ref=temp_ref,
                                        cf=cf, zf=zf)

    dict_boundary_conditions = dict_params["boundary_conditions"]

    dict_el_neg_pos = [dict_boundary_conditions["electrode_neg"],
                       dict_boundary_conditions["electrode_pos"]]

    species_neg_pos = Array{Any}(undef,2)
    for (index, dict_el) in enumerate(dict_el_neg_pos)
        species_names = [species["name"] for species in dict_el["species"]]
        num_species = length(species_names)
        concentrations = reshape([dict_to_unitful(species["concentration"]) for species in dict_el["species"]], (num_species, 1)) ./scales.C0 .|> NoUnits
        species_neg_pos[index] = AxisArray(concentrations; row=species_names, col=["concentration"])
    end

    temp_amb = dict_to_unitful(dict_boundary_conditions["ambient_temperature"])/scales.TEMP0 |> NoUnits
    p_in_neg = dict_to_unitful(dict_boundary_conditions["pressure_inlet_neg"])/scales.P0 |> NoUnits
    v_out_neg = dict_to_unitful(dict_boundary_conditions["velocity_outlet_neg"])/scales.VEL0 |> NoUnits
    #p_out_neg = dict_boundary_conditions["pressure_outlet_neg"] #FIXME
    p_in_pos = dict_to_unitful(dict_boundary_conditions["pressure_inlet_pos"])/scales.P0 |> NoUnits
    v_out_pos = dict_to_unitful(dict_boundary_conditions["velocity_outlet_pos"])/scales.VEL0 |> NoUnits
    #p_out_pos = dict_boundary_conditions["pressure_outlet_pos"] #FIXME

    ϕₛ_neg = dict_to_unitful(dict_boundary_conditions["voltage_neg"])/scales.V0 |> NoUnits
    ϕₛ_pos = dict_to_unitful(dict_boundary_conditions["voltage_pos"])/scales.V0 |> NoUnits

    boundary_conditions = BoundaryConditions{RealType}(species_neg=species_neg_pos[1],
                                          species_pos=species_neg_pos[2],
                                          temp_amb=temp_amb,
                                          p_in_neg=p_in_neg, p_in_pos=p_in_pos,
                                          #p_out_neg=p_out_neg, p_out_pos=p_out_pos, #FIXME
                                          v_out_neg=v_out_neg, v_out_pos=v_out_pos,
                                          ϕₛ_neg=ϕₛ_neg, ϕₛ_pos=ϕₛ_pos)

    dict_discr_params = dict_params["discretization_parameters"]
    spatial_discr = dict_discr_params["spatial_discretization"]
    temporal_discr = dict_discr_params["temporal_discretization"]
    relative_tol = dict_discr_params["relative_tolerance"]
    discr_params = DiscretizationParameters{RealType}(
                                        spatial_discretization=spatial_discr,
                                        temporal_discretization=temporal_discr,
                                        relative_tolerance=relative_tol)
    dict_study_params = dict_params["study_parameters"]

    dict_polarization_params = dict_study_params["polarization"]
    voltage_start = dict_polarization_params["voltage_start"]
    if voltage_start == "ocv"
        voltage_start = NaN # TODO: use flags to indicate type of voltage start / step / end
    else
        voltage_start = dict_to_unitful(dict_polarization_params["voltage_start"])/scales.V0 |> NoUnits
    end
    voltage_stop = dict_to_unitful(dict_polarization_params["voltage_stop"])/scales.V0 |> NoUnits
    voltage_step = dict_to_unitful(dict_polarization_params["voltage_step"])/scales.V0 |> NoUnits
    output_folder = dict_polarization_params["output_folder"]
    polarization_params = PolarizationParameters{RealType}(voltage_start=voltage_start,
                                                          voltage_stop=voltage_stop,
                                                          voltage_step=voltage_step,
                                                          output_folder=output_folder)

    non_isothermal = dict_study_params["non_isothermal"]
    migration = dict_study_params["migration"]

    study_params = StudyParameters{RealType}(polarization=polarization_params,
                                             non_isothermal=non_isothermal,
                                             migration=migration)

    scaling_params = ScalingParameters{RealType}(scales)

    dict_electrolyte_params = dict_model_params["electrolyte"]
    dict_species = dict_electrolyte_params["species"]
    species_names = [species["name"] for species in dict_species]
    num_species = length(species_names) # consider electroneutrality
    colnames = ["charge", "molar_mass", "diffusivity"]
    reference_units = [1.0*NoUnits, scales.ρ0/scales.C0, scales.D0]
    num_cols = length(colnames)
    species = AxisArray(Matrix{RealType}(undef, num_species, num_cols);
                                  row=species_names, col=colnames)

    for idx_species in eachindex(species_names)
        for idx_col in eachindex(colnames)
            dict_quantity = dict_species[idx_species][colnames[idx_col]]
            species[row=idx_species, col=idx_col] =
                    NoUnits(dict_to_unitful(dict_quantity) / reference_units[idx_col])
        end
    end

    λₜ = dict_to_unitful(dict_electrolyte_params["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = dict_to_unitful(dict_electrolyte_params["thermal_vol_capacity"])/scales.CPv0 |> NoUnits

    electrolyte = ElectrolyteParams{RealType}(species=species, λₜ=λₜ, cpᵥ=cpᵥ)

    system_variables = init_system_variables(num_species - 1;
                                        non_isothermal=study_params.non_isothermal)

    model_params = ModelParameters(geom=geom,
                                    el_neg=el_neg_pos[1],
                                    el_pos=el_neg_pos[2],
                                    cc_neg=cc_neg_pos[1],
                                    cc_pos=cc_neg_pos[2],
                                    sep=sep,
                                    boundary=boundary_conditions,
                                    discr=discr_params,
                                    study=study_params,
                                    scales=scales,
                                    scaling_params=scaling_params,
                                    electrolyte=electrolyte,
                                    var=system_variables)
    return model_params
end


function init_system_variables(num_species; non_isothermal=true)
    domains = domain_symbols()
    num_domains = length(domains)
    domain_def = domain_definitions_table(num_species, non_isothermal)
    var_indices = subdomain_variable_indices(domain_def, non_isothermal)
    #subdomain_var_id_to_domain_ids = subdomain_variable_id_to_domain_ids(var_indices)

    var_symbols = variable_symbols()
    var_ext_symbols = variable_symbols(num_species)

    vec_domain_variables = Vector{DomainVariables{Int64, Vector{Int64}}}(undef, num_domains)

    for idx_subdomain in eachindex(vec_domain_variables)
        indices = var_indices[cols=idx_subdomain]

        vec = Vector{Union{Int64,Vector{Int64}}}(undef, length(var_symbols))
        for idx in eachindex(var_symbols)
            ind = collect(indices[var_ext_symbols .== var_symbols[idx]])
            vec[idx] = length(ind)==1 ? ind[1] : ind
        end
        var_symbols_to_indices = NamedTuple{tuple(var_symbols...)}(tuple(vec...))
        vec_domain_variables[idx_subdomain] = DomainVariables(var_symbols_to_indices...)
    end

    @infiltrate

    return NamedTuple{tuple(domains...)}(vec_domain_variables)

    # return (NamedTuple{tuple(domains...)}(vec_domain_variables),
    #         subdomain_var_id_to_domain_ids)
end
