# function system(grid, physics_data)
#     unknown_storage=:sparse
#     sys = VoronoiFVM.System(grid, physics_2d(physics_data), unknown_storage=unknown_storage)

#     var_id_to_domain_ids = physics_data.subdomain_var_id_to_domain_ids
#     for var_id in keys(var_id_to_domain_ids)
#         enable_species!(sys, var_id, var_id_to_domain_ids[var_id])
#     end

#     system_idx = unknown_indices(unknowns(sys))
#     physics_data = (physics_data...,system_idx=system_idx)
#     physics!(sys, physics_2d(physics_data))

#     return sys
# end

function system(grid, physics_data)
    unknown_storage=:sparse
    sys = VoronoiFVM.System(grid, physics_2d(physics_data), unknown_storage=unknown_storage)

    # specify the domains of definition of the variables
    # in the reduced interface modelling approach each variable x is split into
    # a variable x_l and x_r defined over the domain left and right of the interface, respectively
    enable_species!(sys, Int(p_l), [Int(domain_id_el_neg)])
    enable_species!(sys, Int(p_r), [Int(domain_id_el_pos)])
    enable_species!(sys, Int(ϕₛ_l), [Int(domain_id_cc_neg), Int(domain_id_el_neg)])
    enable_species!(sys, Int(ϕₛ_r), [Int(domain_id_el_pos), Int(domain_id_cc_pos)])
    enable_species!(sys, Int(ϕₗ_l), [Int(domain_id_el_neg)])
    enable_species!(sys, Int(ϕₗ_r), [Int(domain_id_el_pos)])

    enable_species!(sys, Int(c_ox_neg_l), [Int(domain_id_el_neg)])
    enable_species!(sys, Int(c_ox_neg_r), [Int(domain_id_el_pos)])
    enable_species!(sys, Int(c_red_neg_l), [Int(domain_id_el_neg)])
    enable_species!(sys, Int(c_red_neg_r), [Int(domain_id_el_pos)])

    enable_species!(sys, Int(c_ox_pos_l), [Int(domain_id_el_neg)])
    enable_species!(sys, Int(c_ox_pos_r), [Int(domain_id_el_pos)])
    enable_species!(sys, Int(c_red_pos_l), [Int(domain_id_el_neg)])
    enable_species!(sys, Int(c_red_pos_r), [Int(domain_id_el_pos)])

    if physics_data.study.non_isothermal
        println("Enable Heat Balance Equation")
        enable_species!(sys, Int(temp_l), [Int(domain_id_cc_neg), Int(domain_id_el_neg)])
        enable_species!(sys, Int(temp_r), [Int(domain_id_el_pos), Int(domain_id_cc_pos)])
        enable_boundary_species!(sys, Int(temp_i), [Int(boundary_id_el_sep_neg)])
    end

    system_idx = unknown_indices(unknowns(sys))
    physics_data = (physics_data...,system_idx=system_idx)
    physics!(sys, physics_2d(physics_data))

    return sys
end

function initial_condition(system, grid, initial_data, non_isothermal=false)
    # Create a solution array
    inival = unknowns(system)
    inival .= 0.0

    node_coords = grid.components[ExtendableGrids.Coordinates]

    # for idx_coords in size(node_coords)[2]
    #     coords = vec(node_coords[:, idx_coords])

    #     inival[1, idx_coords] .= voltage_neg
    #     inival[2, idx_coords] .= voltage_pos
    #     inival[3, idx_coords] .= 0.0
    #     inival[4, idx_coords] .= initial_species_concentration(coords)
    #     inival[5, idx_coords] .= initial_species_concentration(coords)
    #     inival[6, idx_coords] .= initial_species_concentration(coords)
    #     inival[7, idx_coords] .= initial_species_concentration(coords)
    # end

    cell_nodes = grid.components[ExtendableGrids.CellNodes]
    cell_regions = grid.components[ExtendableGrids.CellRegions]
    num_variables = non_isothermal ? length(instances(Variables)) : length(instances(Variables))-3
    values = zeros(num_variables)
    for idx_cell in 1:size(cell_nodes)[2]
        for idx_node in 1:size(cell_nodes)[1]
            node = cell_nodes[idx_node, idx_cell]
            coords = node_coords[:,node]
            region = cell_regions[node]
            initial_condition!(values, region, coords, initial_data, non_isothermal)
            inival[:, node] .= values
        end
    end
    return inival
end

function initial_condition!(values, region,
                            coords, initial_data,
                            non_isothermal=false)
    if region == Int(domain_id_cc_neg) # negative current collector
        values[3] = initial_data.ϕₛ_l
        if non_isothermal
            values[15] = initial_data.temp_l
        end
    end
    if region == Int(domain_id_el_neg) # negative electrode
        values[1] = initial_data.p_l
        values[3] = initial_data.ϕₛ_l
        values[5] = initial_data.ϕₗ_l
        values[7] = initial_data.c_ox_neg_l
        values[9] = initial_data.c_red_neg_l
        values[11] = initial_data.c_ox_pos_l
        values[13] = initial_data.c_red_pos_l
        if non_isothermal
            values[15] = initial_data.temp_l
        end
    end
    if region == Int(boundary_id_el_sep_neg) && non_isothermal # separator interface boundary
        values[17] = initial_data.temp_i
    end
    if region == Int(domain_id_el_pos) # positive electrode
        values[2] = initial_data.p_r
        values[4] = initial_data.ϕₛ_r
        values[6] = initial_data.ϕₗ_r
        values[8] = initial_data.c_ox_neg_r
        values[10] = initial_data.c_red_neg_r
        values[12] = initial_data.c_ox_pos_r
        values[14] = initial_data.c_red_pos_r
        if non_isothermal
            values[16] = initial_data.temp_r
        end
    end
    if region == Int(domain_id_cc_pos) # positive current collector
        values[4] = initial_data.ϕₛ_r
        if non_isothermal
            values[16] = initial_data.temp_r
        end
    end
end
