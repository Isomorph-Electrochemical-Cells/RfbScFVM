function read_config_file(input_file_path)
    dict_input_parameters = JSON.parsefile(input_file_path; dicttype=Dict, inttype=IntType,
                                           use_mmap=true)
    return dict_input_parameters
end

function unitful(dict_quantity)
    RealType(dict_quantity["value"]) * uparse(dict_quantity["unit"])
end

function unitful_to_dict(quantity)
    Dict("value"=>ustrip(quantity), "unit"=>string((unit(quantity))))
end

function unitful_to_dict(quantity::S) where {S <: AbstractVector}
    Dict("value"=>ustrip(quantity), "unit"=>string((unit(first(quantity)))))
end

function unitful_to_dict(quantity, unit)
    Dict("value"=>ustrip(quantity), "unit"=>unit)
end

function set_global_logging_level(logging_level::String)
    level = nothing
    if logging_level == "debug"
        level = Logging.Debug
    elseif logging_level == "info"
        level = Logging.Info
    elseif logging_level == "warning"
        level = Logging.Warn
    elseif logging_level == "error"
        level = Logging.Error
    else
        @error "unknown logging level: " logging_level
    end

    if level !== nothing
        logger = ConsoleLogger(stdout, level)
        global_logger(logger)
    end
end

function preprocess_parameters(dict_params)
    dict_model_params = dict_params["model_parameters"]

    # geometry parameters
    geom, lx = geometry_parameters(dict_model_params["geometry"])

    # mesh parameters
    mesh = mesh_parameters(dict_params["discretization_parameters"], geom)

    # define characteristic scales
    dict_boundary_conditions = dict_params["boundary_conditions"]
    temp_amb = unitful(dict_boundary_conditions["ambient_temperature"])
    scales = CharacteristicScales{RealType}(L0=lx, TEMP0=temp_amb)

    # negative and positive current collectors
    cc_neg = cc_parameters(dict_model_params["current_collector_neg"], scales)
    cc_pos = cc_parameters(dict_model_params["current_collector_pos"], scales)




    dict_discr_params = dict_params["discretization_parameters"]
    spatial_discr = dict_discr_params["spatial_discretization"]
    temporal_discr = dict_discr_params["temporal_discretization"]
    maxiters = dict_discr_params["max_iterations"]
    abstol = dict_discr_params["absolute_tolerance"]
    reltol = dict_discr_params["relative_tolerance"]

    discr_params = DiscretizationParameters{RealType}(spatial_discr=spatial_discr,
                                                      temporal_discr=temporal_discr,
                                                      maxiters=maxiters,
                                                      abstol=abstol,
                                                      reltol=reltol)
    dict_study_params = dict_params["study_parameters"]

    dict_polarization_params = dict_study_params["polarization"]
    voltage_start = dict_polarization_params["voltage_start"]
    if voltage_start == "ocv"
        voltage_start = NaN # TODO: use flags to indicate type of voltage start / step / end
    else
        voltage_start = unitful(dict_polarization_params["voltage_start"])/scales.V0 |> NoUnits
    end
    voltage_stop = unitful(dict_polarization_params["voltage_stop"])/scales.V0 |> NoUnits
    voltage_step = unitful(dict_polarization_params["voltage_step"])/scales.V0 |> NoUnits
    polarization_params = PolarizationParameters{RealType}(voltage_start=voltage_start,
                                                          voltage_stop=voltage_stop,
                                                          voltage_step=voltage_step)
    non_isothermal = dict_study_params["non_isothermal"]
    logging_level = dict_study_params["logging_level"]
    set_global_logging_level(logging_level)

    output_folder = dict_study_params["output_folder"]
    output_file_name = dict_study_params["output_file_name"]
    generate_figures = dict_study_params["generate_figures"]

    study_params = StudyParameters{RealType}(polarization=polarization_params,
                                             non_isothermal=non_isothermal,
                                             output_folder=output_folder,
                                             output_file_name=output_file_name,
                                             generate_figures=generate_figures,
                                             logging_level=logging_level)

    scaling_params = ScalingParameters{RealType}(scales)

    # electrolyte parameters
    electrolyte = electrolyte_parameters(dict_model_params["electrolyte"], scales)

    # boundary parameters
    boundary = boundary_parameters(dict_params["boundary_conditions"], electrolyte, scales)

    # electrode Parameters
    el_neg = electrode_parameters(dict_model_params["electrode_neg"], electrolyte, scales)
    el_pos = electrode_parameters(dict_model_params["electrode_pos"], electrolyte, scales)

    # separator Parameters
    sep = separator_parameters(dict_model_params["separator"], electrolyte, scales)

    num_free_species = count_free_species(electrolyte.species)
    system_variables = init_system_variables(num_free_species;
                                             non_isothermal=study_params.non_isothermal)

    dom = (cc_neg = 1, el_neg = 2, el_pos = 3, cc_pos = 4, sep = 5)
    bnd = (sep = 5, cc_neg_left = 6, cc_neg_el_neg = 7, el_pos_cc_pos = 8,
               cc_pos_right = 9, el_neg_inflow = 10, el_neg_outflow = 11,
               el_pos_inflow = 12, el_pos_outflow = 13 )

    nodes, weights = init_quadrature_rule([0.0, geom.ly_cell]; num_nodes=2)
    quadrature = Quadrature(nodes=nodes, weights=weights)
    results = SystemResults{RealType}()

    model_params = (geom=geom,
                    mesh=mesh,
                    el_neg=el_neg,
                    el_pos=el_pos,
                    cc_neg=cc_neg,
                    cc_pos=cc_pos,
                    sep=sep,
                    boundary=boundary,
                    discr=discr_params,
                    study=study_params,
                    scales=scales,
                    scaling_params=scaling_params,
                    electrolyte=electrolyte,
                    idx=system_variables,
                    results=results,
                    dom=dom,
                    bnd=bnd,
                    quadrature=quadrature)

    return model_params
end

count_all_species(electrolyte_species) = size(electrolyte_species)[1]
count_free_species(electrolyte_species) = size(electrolyte_species)[1] - 1

function boundary_species_concentrations(dict_species, electrolyte, scales)
    species_names = [species["name"] for species in dict_species]
    num_free_species = length(species_names)
    @assert num_free_species == count_free_species(electrolyte.species)

    idx_perm = indexin(species_names, electrolyte.species.row)
    species_names = species_names[idx_perm]

    concentrations = [unitful(species["concentration"]) for species in dict_species]
    concentrations = NoUnits.(concentrations ./ scales.C0)[idx_perm]

    KeyedArray(SVector{num_free_species, RealType}(concentrations); row=species_names)
end

function boundary_parameters(dict_boundary, electrolyte, scales)
    species_neg = boundary_species_concentrations(dict_boundary["species_inlet_neg"],
                                                  electrolyte, scales)
    species_pos = boundary_species_concentrations(dict_boundary["species_inlet_pos"],
                                                  electrolyte, scales)

    temp_amb = unitful(dict_boundary["ambient_temperature"])/scales.TEMP0 |> NoUnits
    p_in_neg = unitful(dict_boundary["pressure_inlet_neg"])/scales.P0 |> NoUnits
    v_out_neg = unitful(dict_boundary["velocity_outlet_neg"])/scales.VEL0 |> NoUnits
    p_in_pos = unitful(dict_boundary["pressure_inlet_pos"])/scales.P0 |> NoUnits
    v_out_pos = unitful(dict_boundary["velocity_outlet_pos"])/scales.VEL0 |> NoUnits

    ϕₛ_neg = unitful(dict_boundary["voltage_neg"])/scales.V0 |> NoUnits
    ϕₛ_pos = unitful(dict_boundary["voltage_pos"])/scales.V0 |> NoUnits

    BoundaryConditions(species_neg=species_neg, species_pos=species_pos,
                       temp_amb=temp_amb,
                       p_in_neg=p_in_neg, p_in_pos=p_in_pos,
                       v_out_neg=v_out_neg, v_out_pos=v_out_pos,
                       ϕₛ_neg=Base.RefValue(ϕₛ_neg),
                       ϕₛ_pos=Base.RefValue(ϕₛ_pos))
end

function cc_parameters(dict_cc, scales)
    σₑ = unitful(dict_cc["electric_conductivity"])/scales.σ0 |> NoUnits
    hₜ = unitful(dict_cc["heat_transfer_coefficient"])/(scales.λ0/scales.L0) |> NoUnits
    λₜ = unitful(dict_cc["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = unitful(dict_cc["thermal_vol_capacity"])/scales.CPv0 |> NoUnits

    return CurrentCollectorParameters{RealType}(σₑ=σₑ, hₜ=hₜ, λₜ=λₜ, cpᵥ=cpᵥ)
end

function mesh_parameters(discr_params, geom)
    if discr_params["spatial_discretization"] == "fvm_1d"
        mesh_params = discr_params["mesh"]
        mesh_sizes = mesh_params["relative_mesh_sizes"]
        hx_rel_cc_neg = mesh_sizes["hx_cc_neg"]
        hx_rel_el_neg = mesh_sizes["hx_el_neg"]
        hx_rel_el_pos = mesh_sizes["hx_el_pos"]
        hx_rel_cc_pos = mesh_sizes["hx_cc_pos"]

        mesh = mesh_1d(geom;
                       hx_cc_neg=hx_rel_cc_neg*geom.lx_cc_neg,
                       hx_el_neg=hx_rel_el_neg*geom.lx_el_neg,
                       hx_cc_pos=hx_rel_cc_pos*geom.lx_cc_pos,
                       hx_el_pos=hx_rel_el_pos*geom.lx_el_pos)
    elseif discr_params["spatial_discretization"] == "fvm_2d"
        mesh_params = discr_params["mesh"]
        mesh_sizes = mesh_params["relative_mesh_sizes"]
        hx_rel_cc_neg = mesh_sizes["hx_cc_neg"]
        hx_rel_el_neg = mesh_sizes["hx_el_neg"]
        hx_rel_el_pos = mesh_sizes["hx_el_pos"]
        hx_rel_cc_pos = mesh_sizes["hx_cc_pos"]
        hy_rel_cell = mesh_sizes["hy_cell"]

        mesh = mesh_2d(geom;
                       hx_cc_neg=hx_rel_cc_neg*geom.lx_cc_neg,
                       hx_el_neg=hx_rel_el_neg*geom.lx_el_neg,
                       hx_cc_pos=hx_rel_cc_pos*geom.lx_cc_pos,
                       hx_el_pos=hx_rel_el_pos*geom.lx_el_pos,
                       hy_cell=hy_rel_cell*geom.ly_cell)
    end
end

function geometry_parameters(dict_geom)
    lx_cc_neg = unitful(dict_geom["current_collector_neg_thickness"])
    lx_el_neg = unitful(dict_geom["electrode_neg_thickness"])
    lx_sep = unitful(dict_geom["separator_thickness"])
    lx_el_pos = unitful(dict_geom["electrode_pos_thickness"])
    lx_cc_pos = unitful(dict_geom["current_collector_pos_thickness"])
    ly_cell = unitful(dict_geom["cell_height"])

    lx = lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos + lx_cc_pos

    lx_cc_neg = lx_cc_neg / lx |> NoUnits
    lx_el_neg = lx_el_neg / lx |> NoUnits
    lx_sep = lx_sep / lx |> NoUnits
    lx_el_pos = lx_el_pos / lx |> NoUnits
    lx_cc_pos = lx_cc_pos / lx |> NoUnits
    ly_cell = ly_cell / lx |> NoUnits

    (flow_cell_geometry(lx_cc_neg=lx_cc_neg, lx_el_neg=lx_el_neg, lx_sep=lx_sep,
                       lx_el_pos=lx_el_pos, lx_cc_pos=lx_cc_pos, ly_cell=ly_cell), lx)
end

function separator_parameters(dict_sep, electrolyte, scales)
    λₜ = unitful(dict_sep["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = unitful(dict_sep["thermal_vol_capacity"])/scales.CPv0 |> NoUnits
    kh = unitful(dict_sep["hydraulic_permeability"])/scales.KH0 |> NoUnits
    kϕ = unitful(dict_sep["electrokinetic_permeability"])/scales.Kϕ0 |> NoUnits
    μ = unitful(dict_sep["dynamic_viscosity"])/scales.μ0 |> NoUnits

    dict_fixed_ionic_groups = dict_sep["fixed_ionic_groups"]
    cf = unitful(dict_fixed_ionic_groups["concentration"])/scales.C0 |> NoUnits
    zf = unitful(dict_fixed_ionic_groups["charge"]) |> NoUnits

    dict_species_sep = dict_sep["species"]

    colnames = ["reference_diffusivity", "reference_temperature", "temperature_coefficient"]
    ref_units = [scales.D0, scales.TEMP0, 1/scales.TEMP0]
    num_cols = length(colnames)
    species_names = electrolyte.species.row
    num_species = length(species_names)
    species = KeyedArray(MMatrix{num_species, num_cols, RealType}(undef),
                        row=species_names,col=colnames)

    sep_species_names = [species["name"] for species in dict_species_sep]
    idx_species_names = indexin(sep_species_names, species_names)
    # iterate over species names in the order of species_names
    for idx_species in idx_species_names
        for idx_col in eachindex(colnames)
            dict_quantity = dict_species_sep[idx_species][colnames[idx_col]]
            species[idx_species, idx_col] = NoUnits(unitful(dict_quantity) /
                                            ref_units[idx_col])
        end
    end

    SeparatorParameters(λₜ=λₜ, cpᵥ=cpᵥ, kh=kh, kϕ=kϕ, μ=μ, cf=cf, zf=zf, species=species)
end

function electrolyte_parameters(dict_electrolyte, scales)
    dict_species = dict_electrolyte["species"]
    species_names = [species["name"] for species in dict_species]
    num_species = length(species_names)

    # identify eliminated species and check that it is contained in species_names
    eliminated_species_name = dict_electrolyte["electroneutrality"]["name"]
    indices_eliminated_species = findall(x->x==eliminated_species_name, species_names)
    error_msg = "species to be determined by the electroneutrality condition not found"
    @assert length(indices_eliminated_species)==1 error_msg
    idx_eliminated_species = indices_eliminated_species[1]

    idx_perm = collect(eachindex(species_names))
    idx_tmp = idx_perm[idx_eliminated_species]
    idx_perm[idx_eliminated_species] = idx_perm[end]
    idx_perm[end] = idx_tmp
    # update order of species names, so that the eliminated species is in the last entry
    species_names = species_names[idx_perm]

    colnames = ["charge", "reference_diffusivity", "reference_temperature",
                "temperature_coefficient"]
    reference_units = [1.0*NoUnits, scales.D0, scales.TEMP0, 1/scales.TEMP0]
    num_cols = length(colnames)
    species = KeyedArray(MMatrix{num_species, num_cols, RealType}(undef),
                         row=species_names, col=colnames)
    for idx_species in eachindex(species_names)
        for idx_col in eachindex(colnames)
            dict_quantity = dict_species[idx_perm[idx_species]][colnames[idx_col]]
            species[idx_species, idx_col] =
                    NoUnits(unitful(dict_quantity) / reference_units[idx_col])
        end
    end

    λₜ = unitful(dict_electrolyte["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = unitful(dict_electrolyte["thermal_vol_capacity"])/scales.CPv0 |> NoUnits

    dict_hom_reactions = dict_electrolyte["homogeneous_reactions"]
    hom_reactions = homogeneous_reaction_params(dict_hom_reactions, species_names, scales)

    ElectrolyteParams(species=species, λₜ=λₜ, cpᵥ=cpᵥ, hom_reactions=hom_reactions)
end

function electrode_parameters(dict_electrode, electrolyte, scales)

    λₜ = unitful(dict_electrode["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = unitful(dict_electrode["thermal_vol_capacity"])/scales.CPv0 |> NoUnits
    εₗ = unitful(dict_electrode["porosity"]) |> NoUnits
    kₕ = unitful(dict_electrode["hydraulic_permeability"])/scales.KH0 |> NoUnits
    dict_mass_transfer_model = dict_electrode["mass_transfer_model"]
    mass_transfer_factor = unitful(dict_mass_transfer_model["mass_transfer_factor"])
    mass_transfer_exp = NoUnits(unitful(dict_mass_transfer_model["mass_transfer_exponent"]))
    mass_transfer_exp = Rational{Int32}(mass_transfer_exp) # convert float to rational
    sh_factor = RealType(NoUnits(mass_transfer_factor*(scales.LP0/scales.D0)*scales.VEL0^mass_transfer_exp))
    sh_exponent = RealType(mass_transfer_exp)
    σₑ = unitful(dict_electrode["electric_conductivity"])/scales.σ0 |> NoUnits
    aᵥ = unitful(dict_electrode["specific_surface_area"])*scales.LP0 |> NoUnits
    μ = unitful(dict_electrode["dynamic_viscosity"])/scales.μ0 |> NoUnits
    Cₐ = unitful(dict_electrode["areal_capacitance"])/scales.Ca0 |> NoUnits
    Cᵥ = Cₐ * aᵥ # volumetric capacitance
    diff_model = dict_electrode["effective_diffusivity_model"]
    Deff_const = unitful(diff_model["constant_coefficient"]) |> NoUnits
    Deff_linear = unitful(diff_model["linear_coefficient"])/(scales.L0/scales.D0) |> NoUnits
    Deff_quadratic = unitful(diff_model["quadratic_coefficient"])/(scales.L0^2/scales.D0^2) |> NoUnits

    vec_reactions = dict_electrode["reactions"]
    num_free_species = count_free_species(electrolyte.species)
    vec_reaction_params = ReactionParameters{RealType, num_free_species}[]
    for idx_reaction in eachindex(vec_reactions)

        reaction = vec_reactions[idx_reaction]
        reaction_name = reaction["name"]
        reduction_potential = reaction["standard_reduction_potential"]
        Δϕ₀ = unitful(reduction_potential["standard_value"])/scales.V0 |> NoUnits
        ∂Δϕ₀_∂temp = unitful(reduction_potential["temperature_coefficient"])/(scales.V0/scales.TEMP0) |> NoUnits

        Δs = unitful(reaction["entropy_change"])/scales.S0 |> NoUnits
        dict_kinetics = reaction["kinetics"]
        α = unitful(dict_kinetics["transfer_coefficient"]) |> NoUnits
        ki = unitful(dict_kinetics["rate_constant"])/scales.K0 |> NoUnits # kinetic number
        temp_ref = unitful(reaction["reference_temperature"])/scales.TEMP0 |> NoUnits

        species_names = electrolyte.species.row
        eliminated_species_name = species_names[end] # name of eliminated species

        vec_coeffs = reaction["stoichiometric_coefficients"]
        ν_coeff = @MVector zeros(IntType, num_free_species)
        ν_el = idx_red = idx_ox = zero(IntType)
        for idx_species in eachindex(vec_coeffs)
            name = vec_coeffs[idx_species]["name"]
            if name == "e"
                ν_el = vec_coeffs[idx_species]["nu"]
                continue
            end
            vec_idx_el = IntType.(findall(==(name), species_names))
            @assert length(vec_idx_el) == 1 "no unique corresponding electrolyte species could be determined for " * name
            idx_el = vec_idx_el[1]

            if vec_coeffs[idx_species]["name"] != eliminated_species_name
                ν_coeff[idx_el] += vec_coeffs[idx_species]["nu"]
            end

            if haskey(vec_coeffs[idx_species], "state")
                if vec_coeffs[idx_species]["state"] == "ox"
                    idx_ox = idx_el
                elseif vec_coeffs[idx_species]["state"] == "red"
                    idx_red = idx_el
                end
            end
        end
        @assert ν_el != 0 "number of exchanged electrons must be non-zero"

        reaction_params = ReactionParameters(name=reaction_name,
                                Δϕ₀=Δϕ₀, ∂Δϕ₀_∂temp=∂Δϕ₀_∂temp,
                                Δs=Δs, α=α, ki=ki, temp_ref=temp_ref,
                                ν_coeff=SVector{num_free_species, IntType}(ν_coeff),
                                ν_el=ν_el, idx_ox=idx_ox, idx_red=idx_red)

        push!(vec_reaction_params, reaction_params)
    end

    ElectrodeParameters(λₜ=λₜ, cpᵥ=cpᵥ, εₗ=εₗ, kₕ=kₕ,
                        μ=μ, σₑ=σₑ, aᵥ=aᵥ, Cᵥ=Cᵥ,
                        Deff_const=Deff_const,
                        Deff_linear=Deff_linear,
                        Deff_quadratic=Deff_quadratic,
                        reactions=vec_reaction_params,
                        sh_factor=sh_factor,
                        sh_exponent=sh_exponent)
end

function homogeneous_reaction_params(homogeneous_reactions, species_names, scales)
    num_species = length(species_names)
    num_reactions = length(homogeneous_reactions)

    # matrix with non-negative stoichiometric coefficients of reactants
    mat_a = @MMatrix zeros(IntType, num_reactions, num_species)
    # matrix with non-negative stoichiometric coefficients of products
    mat_b = @MMatrix zeros(IntType, num_reactions, num_species)
    # vector storing the reaction constants
    vec_k = @MVector zeros(RealType, num_reactions)

    for idx_reaction in eachindex(homogeneous_reactions)
        dict_reaction = homogeneous_reactions[idx_reaction]
        dict_species = dict_reaction["species"]

        for idx_coeff in eachindex(dict_species)
            name = dict_species[idx_coeff]["name"]
            react_coeff = dict_species[idx_coeff]["reactant_coefficient"]
            prod_coeff = dict_species[idx_coeff]["product_coefficient"]

            found_species_idx = findall(x->x==name, species_names)
            @assert length(found_species_idx) == 1 "species " *
                    name * " not known in reaction " * dict_reaction["name"]
            idx_species = found_species_idx[1]
            mat_a[idx_reaction, idx_species] += react_coeff
            mat_b[idx_reaction, idx_species] += prod_coeff
        end
    end

    for idx_reaction in eachindex(homogeneous_reactions)
        dict_reaction = homogeneous_reactions[idx_reaction]
        order = sum(mat_a[idx_reaction,:]) # reaction order
         # reference scale of reaction constant k for a reaction of a given order
        dim_rate_const = scales.D0/(scales.C0^(order-1)*scales.L0^2)
        dim_rate_const_value = unitful(dict_reaction["rate_constant"])
        try
            vec_k[idx_reaction] = dim_rate_const_value / dim_rate_const |> NoUnits
        catch e
            @error "provided unit of rate constant in reaction " * dict_reaction["name"] *
                    " is " * string(unit(dim_rate_const_value)) * ", but should be " *
                    string(unit(dim_rate_const))
            throw(e)
        end
    end

    # eliminate species due to electroneutrality constraint
    idx_counter_species = num_species
    mat_a = SMatrix{num_reactions, num_species-1, IntType}(mat_a[:, (1:end) .∉ idx_counter_species])
    mat_b = SMatrix{num_reactions, num_species-1, IntType}(mat_b[:, (1:end) .∉ idx_counter_species])
    mat_c = SMatrix(transpose(mat_b-mat_a) * diagm(vec_k))
    HomogeneousReactionParams{num_species-1, num_reactions, RealType}(mat_a=mat_a, mat_c=mat_c)
end

function init_quadrature_rule(interval; num_nodes=1)
    nodes, weights = gausslegendre(num_nodes)
    n_scaled = (interval[1] + interval[2]) / 2 .+ nodes .* (interval[2] - interval[1]) / 2
    w_scaled = weights .* (interval[2] - interval[1]) / 2
    (SVector(n_scaled...), SVector(w_scaled...))
end


function init_system_variables(num_species; non_isothermal=true)
    domains = domain_symbols()
    num_domains = length(domains)
    domain_def = domain_definitions_table(num_species, non_isothermal)
    var_indices = subdomain_variable_indices(domain_def, non_isothermal)

    var_symbols = variable_symbols()
    var_ext_symbols = variable_symbols(num_species)

    DomainVariablesType = DomainVariables{IntType, SVector{num_species, IntType}}
    vec_domain_variables = Vector{DomainVariablesType}(undef, num_domains)

    for idx_subdomain in eachindex(vec_domain_variables)
        indices = var_indices[cols=idx_subdomain]

        vec = Vector{Union{IntType,Vector{IntType}}}(undef, length(var_symbols))
        for idx in eachindex(var_symbols)
            ind = collect(indices[var_ext_symbols .== var_symbols[idx]])
            vec[idx] = length(ind)==1 ? ind[1] : ind
        end
        vec_domain_variables[idx_subdomain] = DomainVariablesType(vec...)
    end

    return vec_domain_variables
end
