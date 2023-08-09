function system(grid, data)
    unknown_storage=:sparse
    sys = VoronoiFVM.System(grid, physics_2d(data), unknown_storage=unknown_storage)

    # specify the domains of definition of the variables
    enable_species!(sys, data.idx[data.dom[:el_neg]].p, [data.dom[:el_neg]])
    enable_species!(sys, data.idx[data.dom[:el_pos]].p, [data.dom[:el_pos]])
    enable_species!(sys, data.idx[data.dom[:el_neg]].ϕₛ, [data.dom[:cc_neg],
                                                             data.dom[:el_neg]])
    enable_species!(sys, data.idx[data.dom[:el_pos]].ϕₛ, [data.dom[:el_pos],
                                                              data.dom[:cc_pos]])
    enable_species!(sys, data.idx[data.dom[:el_neg]].ϕₗ, [data.dom[:el_neg]])
    enable_species!(sys, data.idx[data.dom[:el_pos]].ϕₗ, [data.dom[:el_pos]])

    for idx_species in data.idx[data.dom[:el_neg]].c
        enable_species!(sys, idx_species, [data.dom[:el_neg]])
    end

    for idx_species in data.idx[data.dom[:el_pos]].c
        enable_species!(sys, idx_species, [data.dom[:el_pos]])
    end

    if data.study.non_isothermal
        println("Enable Heat Balance Equation")
        temp_l = data.idx[data.dom[:el_neg]].temp
        temp_r = data.idx[data.dom[:el_pos]].temp
        temp_i = data.idx[data.dom[:sep]].temp
        enable_species!(sys, temp_l, [data.dom[:cc_neg], data.dom[:el_neg]])
        enable_species!(sys, temp_r, [data.dom[:cc_pos], data.dom[:el_pos]])
        enable_boundary_species!(sys, temp_i, [data.bnd[:sep]])
    end

    physics!(sys, physics_2d(data))

    return sys
end

function initial_condition(system, grid, initial_data, data)
    # Create a solution array
    inival = unknowns(system)
    inival .= 0.0

    node_coords = grid.components[ExtendableGrids.Coordinates]

    cell_nodes = grid.components[ExtendableGrids.CellNodes]
    cell_regions = grid.components[ExtendableGrids.CellRegions]
    non_isothermal = data.study.non_isothermal
    num_species = 2*length(data.electrolyte.species[col=1])
    num_variables = non_isothermal ? 9+num_species : 6+num_species

    values = zeros(num_variables)
    for idx_cell in 1:size(cell_nodes)[2]
        for idx_node in 1:size(cell_nodes)[1]
            node = cell_nodes[idx_node, idx_cell]
            coords = node_coords[:,node]
            region = cell_regions[node]
            initial_condition!(values, region, coords, initial_data, data)
            inival[:, node] .= values
        end
    end
    return inival
end

function initial_condition!(values, region,
                            coords, initial_data,
                            data)
    non_isothermal = data.study.non_isothermal

    # negative or positive current collector
    if region == data.dom[:cc_neg] || region == data.dom[:cc_pos]
        values[data.idx[region].ϕₛ] = initial_data[region].ϕₛ
        if non_isothermal
            values[data.idx[region].temp] = initial_data[region].temp
        end
    end
    # negative or positive electrode
    if region == data.dom[:el_neg] || region == data.dom[:el_pos]
        values[data.idx[region].p] = initial_data[region].p
        values[data.idx[region].ϕₛ] = initial_data[region].ϕₛ
        values[data.idx[region].ϕₗ] = initial_data[region].ϕₗ

        for idx_species in eachindex(data.idx[region].c)
            values[data.idx[region].c[idx_species]] = initial_data[region].c[idx_species]
        end
        if non_isothermal
            values[data.idx[region].temp] = initial_data[region].temp
        end
    end
    # separator interface boundary
    if region == data.bnd[:sep] && non_isothermal
        values[data.idx[region].temp] = initial_data[region].temp
    end
end
