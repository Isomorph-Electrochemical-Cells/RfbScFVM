function parse_json_file(input_file_path)
    dict_input_parameters = JSON.parsefile(input_file_path; dicttype=Dict, inttype=IndexType,
                                           use_mmap=true)
    return dict_input_parameters
end

function read_simulation_parameters(input_file_path)
    dict_params = parse_json_file(input_file_path)
    return preprocess_model_parameters(dict_params)
end

function unitful(dict_quantity)
    RealType(dict_quantity["value"]) * uparse(dict_quantity["unit"])
end

function unitful_to_dict(quantity)
    Dict("value"=>ustrip(quantity), "unit"=>(unit(quantity)))
end

function unitful_to_dict(quantity, unit)
    Dict("value"=>ustrip(quantity), "unit"=>unit)
end

function preprocess_model_parameters(dict_params)
    dict_model_params = dict_params["model_parameters"]

    # geometry
    dict_geom = dict_model_params["geometry"]

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

    mesh_params = dict_params["discretization_parameters"]["mesh"]
    mesh_sizes = mesh_params["relative_mesh_sizes"]
    hy_rel_cell = mesh_sizes["hy_cell"]
    hx_rel_cc_neg = mesh_sizes["hx_cc_neg"]
    hx_rel_el_neg = mesh_sizes["hx_el_neg"]
    hx_rel_el_pos = mesh_sizes["hx_el_pos"]
    hx_rel_cc_pos = mesh_sizes["hx_cc_pos"]

    geom = flow_cell_geometry(lx_cc_neg=lx_cc_neg, lx_el_neg=lx_el_neg,
        lx_sep=lx_sep, lx_el_pos=lx_el_pos, lx_cc_pos=lx_cc_pos, ly_cell=ly_cell)

    mesh = mesh_2d(geom,
                   hx_cc_neg=hx_rel_cc_neg*lx_cc_neg,
                   hx_el_neg=hx_rel_el_neg*lx_el_neg,
                   hx_cc_pos=hx_rel_cc_pos*lx_cc_pos,
                   hx_el_pos=hx_rel_el_pos*lx_el_pos,
                   hy_cell=hy_rel_cell*ly_cell)

    dict_boundary_conditions = dict_params["boundary_conditions"]
    temp_amb = unitful(dict_boundary_conditions["ambient_temperature"])
    # define characteristic scales
    scales = CharacteristicScales{RealType}(L0=lx, TEMP0=temp_amb)

    # negative and positive current collectors
    dict_cc_neg_pos = [dict_model_params["current_collector_neg"],
                       dict_model_params["current_collector_pos"]]

    cc_neg_pos = Array{CurrentCollectorParameters{RealType}}(undef,2)
    for (index, dict_cc) in enumerate(dict_cc_neg_pos)
        σₑ = unitful(dict_cc["electric_conductivity"])/scales.σ0 |> NoUnits
        hₜ = unitful(dict_cc["heat_transfer_coefficient"])/(scales.λ0/scales.L0) |> NoUnits
        λₜ = unitful(dict_cc["thermal_conductivity"])/scales.λ0 |> NoUnits
        cpᵥ = unitful(dict_cc["thermal_vol_capacity"])/scales.CPv0 |> NoUnits

        cc_neg_pos[index] = CurrentCollectorParameters{RealType}(σₑ=σₑ, hₜ=hₜ, λₜ=λₜ, cpᵥ=cpᵥ)
    end

    # # negative and postive electrodes
    dict_el_neg_pos = [dict_model_params["electrode_neg"],
                       dict_model_params["electrode_pos"]]

    el_neg_pos = Array{Any}(undef,2)
    for (index, dict_el) in enumerate(dict_el_neg_pos)
        λₜ = unitful(dict_el["thermal_conductivity"])/scales.λ0 |> NoUnits
        cpᵥ = unitful(dict_el["thermal_vol_capacity"])/scales.CPv0 |> NoUnits
        εₗ = unitful(dict_el["porosity"]) |> NoUnits
        kₕ = unitful(dict_el["hydraulic_permeability"])/scales.KH0 |> NoUnits
        dict_mass_transfer_model = dict_el["mass_transfer_model"]
        mass_transfer_factor = unitful(dict_mass_transfer_model["mass_transfer_factor"])
        mass_transfer_exp = NoUnits(unitful(dict_mass_transfer_model["mass_transfer_exponent"]))
        mass_transfer_exp = Rational{Int32}(mass_transfer_exp) # convert float to rational
        sh_factor = RealType(NoUnits(mass_transfer_factor*(scales.LP0/scales.D0)*scales.VEL0^mass_transfer_exp))
        sh_exponent = RealType(mass_transfer_exp)
        σₑ = unitful(dict_el["electric_conductivity"])/scales.σ0 |> NoUnits
        aᵥ = unitful(dict_el["specific_surface_area"])*scales.LP0 |> NoUnits
        μ = unitful(dict_el["dynamic_viscosity"])/scales.μ0 |> NoUnits
        Cₐ = unitful(dict_el["areal_capacitance"])/scales.Ca0 |> NoUnits
        Cᵥ = Cₐ * aᵥ # volumetric capacitance
        diff_model = dict_el["effective_diffusivity_model"]
        Deff_const = unitful(diff_model["constant_coefficient"]) |> NoUnits
        Deff_linear = unitful(diff_model["linear_coefficient"])/(scales.L0/scales.D0) |> NoUnits
        Deff_quadratic = unitful(diff_model["quadratic_coefficient"])/(scales.L0^2/scales.D0^2) |> NoUnits

        reaction = dict_el["reactions"][1] # TODO: Support multiple reactions
        reaction_name = reaction["name"]
        reduction_potential = reaction["standard_reduction_potential"]
        Δϕ₀ = unitful(reduction_potential["standard_value"])/scales.V0 |> NoUnits
        ∂Δϕ₀_∂temp = unitful(reduction_potential["temperature_coefficient"])/(scales.V0/scales.TEMP0) |> NoUnits

        Δs = unitful(reaction["entropy_change"])/scales.S0 |> NoUnits
        dict_kinetics = reaction["kinetics"]
        α = unitful(dict_kinetics["transfer_coefficient"]) |> NoUnits
        ki = unitful(dict_kinetics["rate_constant"])/scales.K0 |> NoUnits # kinetic number
        temp_ref = unitful(reaction["reference_temperature"])/scales.TEMP0 |> NoUnits


        dict_electrolyte_params = dict_model_params["electrolyte"]
        dict_species = dict_electrolyte_params["species"]
        species_names = [species["name"] for species in dict_species]
        num_species = length(species_names)

        vec_coeffs = reaction["stoichiometric_coefficients"]
        ν_coeff = @MVector zeros(IndexType, num_species)
        ν_el = idx_red = idx_ox = zero(IndexType)
        for idx_species in eachindex(vec_coeffs)
            name = vec_coeffs[idx_species]["name"]
            if name == "e"
                ν_el = vec_coeffs[idx_species]["nu"]
                continue
            end
            vec_idx = IndexType.(findall(==(name), species_names))
            @assert length(vec_idx) == 1 "no unique corresponding electrolyte species could be determined for " * name
            ν_coeff[vec_idx[1]] = vec_coeffs[idx_species]["nu"]

            if haskey(vec_coeffs[idx_species], "state")
                if vec_coeffs[idx_species]["state"] == "ox"
                    idx_ox = vec_idx[1]
                elseif vec_coeffs[idx_species]["state"] == "red"
                    idx_red = vec_idx[1]
                end
            end
        end
        @assert ν_el != 0 "number of exchanged electrons must be non-zero"

        reactions = ReactionParameters(name=reaction_name,
                                      Δϕ₀=Δϕ₀, ∂Δϕ₀_∂temp=∂Δϕ₀_∂temp,
                                      Δs=Δs, α=α, ki=ki, temp_ref=temp_ref,
                                      ν_coeff=SVector{num_species,IndexType}(ν_coeff),
                                      ν_el=ν_el, idx_ox=idx_ox, idx_red=idx_red)

        el_neg_pos[index] = ElectrodeParameters(λₜ=λₜ, cpᵥ=cpᵥ, εₗ=εₗ, kₕ=kₕ,
                                                μ=μ, σₑ=σₑ, aᵥ=aᵥ, Cᵥ=Cᵥ,
                                                Deff_const=Deff_const,
                                                Deff_linear=Deff_linear,
                                                Deff_quadratic=Deff_quadratic,
                                                reactions=reactions,
                                                sh_factor=sh_factor,
                                                sh_exponent=sh_exponent)
    end

    # # separator
    dict_sep = dict_model_params["separator"]
    λₜ = unitful(dict_sep["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = unitful(dict_sep["thermal_vol_capacity"])/scales.CPv0 |> NoUnits
    kh = unitful(dict_sep["hydraulic_permeability"])/scales.KH0 |> NoUnits
    kϕ = unitful(dict_sep["electrokinetic_permeability"])/scales.Kϕ0 |> NoUnits
    μ = unitful(dict_sep["dynamic_viscosity"])/scales.μ0 |> NoUnits
    conductivity = dict_sep["electrical_conductivity"]
    σₑ = unitful(conductivity["reference_value"])/scales.σ0 |> NoUnits
    ∂σₑ∂temp = unitful(conductivity["temperature_coefficient"])/(scales.σ0/scales.TEMP0) |> NoUnits
    temp_ref = unitful(conductivity["reference_temperature"])/(scales.TEMP0) |> NoUnits

    dict_fixed_ionic_groups = dict_sep["fixed_ionic_groups"]
    cf = unitful(dict_fixed_ionic_groups["concentration"])/scales.C0 |> NoUnits
    zf = unitful(dict_fixed_ionic_groups["charge"]) |> NoUnits

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
        concentrations = [unitful(species["concentration"]) for species in dict_el["species"]] ./scales.C0 .|> NoUnits
        species_neg_pos[index] = KeyedArray(SVector{num_species, RealType}(concentrations); row=species_names)
    end

    temp_amb = unitful(dict_boundary_conditions["ambient_temperature"])/scales.TEMP0 |> NoUnits
    p_in_neg = unitful(dict_boundary_conditions["pressure_inlet_neg"])/scales.P0 |> NoUnits
    v_out_neg = unitful(dict_boundary_conditions["velocity_outlet_neg"])/scales.VEL0 |> NoUnits
    #p_out_neg = dict_boundary_conditions["pressure_outlet_neg"] #FIXME
    p_in_pos = unitful(dict_boundary_conditions["pressure_inlet_pos"])/scales.P0 |> NoUnits
    v_out_pos = unitful(dict_boundary_conditions["velocity_outlet_pos"])/scales.VEL0 |> NoUnits
    #p_out_pos = dict_boundary_conditions["pressure_outlet_pos"] #FIXME

    ϕₛ_neg = unitful(dict_boundary_conditions["voltage_neg"])/scales.V0 |> NoUnits
    ϕₛ_pos = unitful(dict_boundary_conditions["voltage_pos"])/scales.V0 |> NoUnits

    boundary_conditions = BoundaryConditions(species_neg=species_neg_pos[1],
                                          species_pos=species_neg_pos[2],
                                          temp_amb=temp_amb,
                                          p_in_neg=p_in_neg, p_in_pos=p_in_pos,
                                          #p_out_neg=p_out_neg, p_out_pos=p_out_pos, #FIXME
                                          v_out_neg=v_out_neg, v_out_pos=v_out_pos,
                                          ϕₛ_neg=Base.RefValue(ϕₛ_neg),
                                          ϕₛ_pos=Base.RefValue(ϕₛ_pos)
                                          )

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
    migration = dict_study_params["migration"]
    output_folder = dict_study_params["output_folder"]
    output_file_name = dict_study_params["output_file_name"]
    generate_figures = dict_study_params["generate_figures"]

    study_params = StudyParameters{RealType}(polarization=polarization_params,
                                             non_isothermal=non_isothermal,
                                             migration=migration,
                                             output_folder=output_folder,
                                             output_file_name=output_file_name,
                                             generate_figures=generate_figures)

    scaling_params = ScalingParameters{RealType}(scales)

    dict_electrolyte_params = dict_model_params["electrolyte"]
    dict_species = dict_electrolyte_params["species"]
    species_names = [species["name"] for species in dict_species]
    num_species = length(species_names)
    eliminated_species_name = dict_electrolyte_params["electroneutrality"]["name"]
    @assert eliminated_species_name in species_names "species to be determined by the electroneutrality condition not found"
    species_variable_names = species_names[species_names .!= eliminated_species_name]
    colnames = ["charge", "molar_mass", "diffusivity"]
    reference_units = [1.0*NoUnits, scales.ρ0/scales.C0, scales.D0]
    num_cols = length(colnames)
    species = KeyedArray(MMatrix{num_species-1,num_cols,RealType}(undef),row=species_variable_names,col=colnames)

    for idx_species in eachindex(species_variable_names)
        for idx_col in eachindex(colnames)
            dict_quantity = dict_species[idx_species][colnames[idx_col]]
            species[idx_species, idx_col] =
                    NoUnits(unitful(dict_quantity) / reference_units[idx_col])
        end
    end
    counter = KeyedArray(MVector{num_cols,RealType}(undef),col=colnames)

    idx_counter_species = findfirst(species_names .== eliminated_species_name)
    for idx_col in eachindex(colnames)
        dict_quantity = dict_species[idx_counter_species][colnames[idx_col]]
        counter[idx_col] = NoUnits(unitful(dict_quantity) / reference_units[idx_col])
    end

    λₜ = unitful(dict_electrolyte_params["thermal_conductivity"])/scales.λ0 |> NoUnits
    cpᵥ = unitful(dict_electrolyte_params["thermal_vol_capacity"])/scales.CPv0 |> NoUnits

    electrolyte = ElectrolyteParams(species=species, counter=counter, λₜ=λₜ, cpᵥ=cpᵥ)

    system_variables = init_system_variables(num_species - 1;
                                        non_isothermal=study_params.non_isothermal)

    dom = (cc_neg = 1, el_neg = 2, el_pos = 3, cc_pos = 4, sep = 5)
    bnd = (sep = 5, cc_neg_left = 6, cc_neg_el_neg = 7, el_pos_cc_pos = 8,
               cc_pos_right = 9, el_neg_inflow = 10, el_neg_outflow = 11,
               el_pos_inflow = 12, el_pos_outflow = 13 )

    model_params = ModelParameters(geom=geom,
                                   mesh=mesh,
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
                                   idx=system_variables,
                                   results=SystemResults{RealType}(),
                                   dom=dom,
                                   bnd=bnd)

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

    DomainVariablesType = DomainVariables{IndexType, SVector{num_species, IndexType}}
    vec_domain_variables = Vector{DomainVariablesType}(undef, num_domains)

    for idx_subdomain in eachindex(vec_domain_variables)
        indices = var_indices[cols=idx_subdomain]

        vec = Vector{Union{IndexType,Vector{IndexType}}}(undef, length(var_symbols))
        for idx in eachindex(var_symbols)
            ind = collect(indices[var_ext_symbols .== var_symbols[idx]])
            vec[idx] = length(ind)==1 ? ind[1] : ind
        end
        vec_domain_variables[idx_subdomain] = DomainVariablesType(vec...)

        #var_symbols_to_indices = NamedTuple{tuple(var_symbols...)}(tuple(vec...))
        #vec_domain_variables[idx_subdomain] = DomainVariables{IndexType, MVector{num_species, IndexType}}(var_symbols_to_indices...)
    end

    return vec_domain_variables #NamedTuple{tuple(domains...)}(vec_domain_variables)

    # return (NamedTuple{tuple(domains...)}(vec_domain_variables),
    #         subdomain_var_id_to_domain_ids)
end
