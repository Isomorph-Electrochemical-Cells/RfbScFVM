function system(grid, data)
    unknown_storage=:sparse

    if dim_space(grid) == 1
        sys = VoronoiFVM.System(grid, physics_1d(data), unknown_storage=unknown_storage)
    elseif dim_space(grid) == 2
        sys = VoronoiFVM.System(grid, physics_2d(data), unknown_storage=unknown_storage)
    end

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
        temp_l = data.idx[data.dom[:el_neg]].temp
        temp_r = data.idx[data.dom[:el_pos]].temp
        temp_i = data.idx[data.dom[:sep]].temp
        enable_species!(sys, temp_l, [data.dom[:cc_neg], data.dom[:el_neg]])
        enable_species!(sys, temp_r, [data.dom[:cc_pos], data.dom[:el_pos]])
        enable_boundary_species!(sys, temp_i, [data.bnd[:sep]])
    end

    if dim_space(grid) == 1
        physics!(sys, physics_1d(data))
    elseif dim_space(grid) == 2
        physics!(sys, physics_2d(data))
    end

    return sys
end

function get_initial_condition(system, grid, initial_data, data)
    # Create a solution array
    inival = unknowns(system)
    inival .= 0.0

    node_coords = grid.components[ExtendableGrids.Coordinates]
    cell_nodes = grid.components[ExtendableGrids.CellNodes]
    cell_regions = grid.components[ExtendableGrids.CellRegions]
    non_isothermal = data.study.non_isothermal
    num_species = 2*(length(data.electrolyte.species[col=1])-1)
    num_variables = non_isothermal ? 9+num_species : 6+num_species

    values = zeros(num_variables)
    for idx_cell in 1:size(cell_nodes)[2]
        for idx_node in 1:size(cell_nodes)[1]
            node = cell_nodes[idx_node, idx_cell]
            coords = node_coords[:,node]
            region = cell_regions[idx_cell] #cell_regions[node] # TODO: IS THIS CORRECT FOR THE 2D CASE, TOO???
            initial_condition!(values, region, coords, initial_data, data)
            inival[:, node] .= values
        end
    end

############## # FIXME
    # num_cells = size(inival)[2]
    #inival[end, Int(ceil(num_cells/2))] = 1.0
    if data.study.non_isothermal
        inival[end, :] .= initial_data[data.dom[:sep]].temp
    end
##############

    return inival
end

function initial_pressure_field(y; p_in=p_in, v=v, μ=μ, kₕ=kₕ)
    p_in - y * μ * v / kₕ
end

function initial_condition!(values, region, coords, initial_data, data)
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
        if occursin("2d", data.discr.spatial_discr)
            v = region == data.dom[:el_neg] ? data.boundary.v_out_neg : data.boundary.v_out_pos
            μ = region == data.dom[:el_neg] ? data.el_neg.μ : data.el_pos.μ
            kₕ = region == data.dom[:el_neg] ? data.el_neg.kₕ : data.el_pos.kₕ
            values[data.idx[region].p] = initial_pressure_field(coords[2];
                                            p_in=initial_data[region].p, v=v, μ=μ, kₕ=kₕ)
        else
            values[data.idx[region].p] = initial_data[region].p
        end

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
    if region == data.dom[:sep] && non_isothermal
        values[data.idx[region].temp] = initial_data[region].temp
    end
end


function project_solutions(solution, sys, grid, subgrids, data, dim=1)
    # TODO: Generalize this implementation to allow for projections of arbitrary subdomains
    # Currently, only the projection on the x-axis (dim=1) is supported!

    coords_x = grid.components[ExtendableGrids.Coordinates][dim, :]
    sp = sortperm(coords_x)

    coords_x_sorted = coords_x[sp]
    solution_sorted = solution[:, sp]
    replace!(solution_sorted, NaN=>0.0)

    coords_x_sorted_unique = unique(coords_x_sorted)
    solution_x = zeros(size(solution_sorted)[1], length(coords_x_sorted_unique))
    count_x = zeros(length(coords_x_sorted_unique))

    x_last = -Inf
    idx_solution_x = 0
    for idx in 1:size(solution_sorted)[2]
        x_current = coords_x_sorted[idx]
        if x_current > x_last + eps(x_current)
            x_last = x_current
            idx_solution_x += 1
        end
        count_x[idx_solution_x] += 1
        solution_x[:, idx_solution_x] .+= solution_sorted[:, idx]
    end
    @assert idx_solution_x == size(solution_x)[2]
    @assert isapprox(sum(count_x), length(coords_x))
    @assert isapprox(sum(count_x), size(solution_sorted)[2])

    for idx_row in 1:size(solution_x)[1]
        solution_x[idx_row, :] ./= count_x
    end

    solutions_on_subgrids = Dict{Symbol, Tuple{Vector{Float64}, Matrix{Float64}}}()

    for key_subgrid in keys(subgrids)
        subgrid = subgrids[key_subgrid]
        coords_x_subgrid = vec(subgrid.components[ExtendableGrids.Coordinates][dim, :])
        coords_x_subgrid_sorted_unique = unique(sort(coords_x_subgrid))
        idx_coords_x_subgrid = 1
        solution_x_subgrid = zeros(size(solution_sorted)[1], length(coords_x_subgrid_sorted_unique))
        for idx_coords_x in eachindex(coords_x_sorted)
            if idx_coords_x_subgrid <= length(coords_x_subgrid_sorted_unique) &&
                isapprox(coords_x[idx_coords_x], coords_x_subgrid_sorted_unique[idx_coords_x_subgrid])
                solution_x_subgrid[:, idx_coords_x_subgrid] = solution_x[:, idx_coords_x]
                idx_coords_x_subgrid += 1
            end
        end

        solutions_on_subgrids[key_subgrid] = (coords_x_subgrid_sorted_unique,
                                              solution_x_subgrid)
    end

    return (coords_x_sorted_unique, solution_x, solutions_on_subgrids)
end
