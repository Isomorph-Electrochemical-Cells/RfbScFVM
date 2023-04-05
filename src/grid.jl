@enum SubDomains begin
    domain_id_cc_neg = 1
    domain_id_el_neg = 2
    domain_id_el_pos = 3
    domain_id_cc_pos = 4
    boundary_id_cc_neg = 5
    boundary_id_cc_el_neg = 6
    boundary_id_el_sep_neg = 7
    boundary_id_el_cc_pos = 8
    boundary_id_cc_pos = 9
    boundary_id_el_neg_inflow = 10
    boundary_id_el_neg_outflow = 11
    boundary_id_el_pos_inflow = 12
    boundary_id_el_pos_outflow = 13
end

function generate_grid_nodes(domain, h)
    num_subdomain_nodes = length(h)
    subdomain_nodes = collect(range(domain..., length=num_subdomain_nodes))
    nodes = []
    for idx_subdomain in 1:num_subdomain_nodes-1
        nodes_subdomain = geomspace(subdomain_nodes[idx_subdomain],
                                    subdomain_nodes[idx_subdomain+1],
                                    h[idx_subdomain], h[idx_subdomain+1])
        nodes = isempty(nodes) ? nodes_subdomain : glue(nodes, nodes_subdomain)
    end
    return nodes
end

function create_grid_2d(cell_geometry::FlowCellGeometry2D{T};
                        reduced_membrane_model = false) where {T<:AbstractFloat}

    lx_cc_neg = cell_geometry.lx_cc_neg
    lx_el_neg = cell_geometry.lx_el_neg
    lx_sep = reduced_membrane_model ? 0.0 : cell_geometry.lx_sep
    lx_sep_physical = cell_geometry.lx_sep
    lx_el_pos = cell_geometry.lx_el_pos
    lx_cc_pos = cell_geometry.lx_cc_pos
    ly = cell_geometry.ly_cell

    x_left = 0.0
    x_cc_el_neg = lx_cc_neg
    x_el_sep_neg = lx_cc_neg + lx_el_neg
    x_sep_el_pos = lx_cc_neg + lx_el_neg + lx_sep
    x_el_cc_pos = lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos
    x_right =  lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos + lx_cc_pos

    y_bottom = 0.0
    y_top = ly

    domain_cc_neg = (x_left, x_cc_el_neg)
    domain_el_neg = (x_cc_el_neg, x_el_sep_neg)
    domain_sep = (x_el_sep_neg, x_sep_el_pos)
    domain_el_pos = (x_sep_el_pos, x_el_cc_pos)
    domain_cc_pos = (x_el_cc_pos, x_right)
    domain_x = (x_left, x_right)
    domain_y = (y_bottom, y_top)

    # num_cells_y = cell_geometry.num_cells_y
    # h_cc_neg = lx_cc_neg/cell_geometry.num_cells_cc_neg_x
    # h_sep = lx_sep_physical/cell_geometry.num_cells_sep_x
    # h_cc_pos = lx_cc_pos/cell_geometry.num_cells_cc_pos_x
    # h_y = ly/num_cells_y
    # h_y_max = h_y
    # h_y_min = h_y/10

    nodes_cc_neg = generate_grid_nodes(domain_cc_neg, cell_geometry.hx_cc_neg)
    nodes_el_neg = generate_grid_nodes(domain_el_neg, cell_geometry.hx_el_neg)
    nodes_el_pos = generate_grid_nodes(domain_el_pos, cell_geometry.hx_el_pos)
    nodes_cc_pos = generate_grid_nodes(domain_cc_pos, cell_geometry.hx_cc_pos)
    nodes_y = generate_grid_nodes(domain_y, cell_geometry.hy_cell)


    # nodes_cc_neg = collect(range(domain_cc_neg..., step=h_cc_neg))
    # nodes_el_neg = geomspace(domain_el_neg..., h_cc_neg, h_sep)
    # if !reduced_membrane_model
    #     nodes_sep = collect(range(domain_sep..., step=h_sep))
    # end
    # nodes_el_pos = geomspace(domain_el_pos..., h_sep, h_cc_neg)
    # nodes_cc_pos = collect(range(domain_cc_pos..., step=h_cc_pos))

    all_nodes_x = glue(nodes_cc_neg, nodes_el_neg)
    all_nodes_x = glue(all_nodes_x, nodes_el_pos)
    all_nodes_x = glue(all_nodes_x, nodes_cc_pos)

    # equidistant discretization in y direction
    # nodes_y = collect(range(domain_y..., step=h_y))
    # y_middle = (y_bottom+y_top)/2
    # domain_y_lower = (y_bottom, y_middle-h_y_max/2)
    # domain_y_upper = (y_middle+h_y_max/2, y_top)
    # nodes_y_lower = geomspace(domain_y_lower..., h_y_min, h_y_max)
    # nodes_y_upper = geomspace(domain_y_upper..., h_y_max, h_y_min)
    # nodes_y = vcat(nodes_y_lower, nodes_y_upper)


    # generate grid
    grid = simplexgrid(all_nodes_x, nodes_y)

    # assign IDs to boundaries
    bfacemask!(grid,
                [0.0, 0.0],
                [0.0, ly],
                Int(boundary_id_cc_neg))
    bfacemask!(grid,
                [lx_cc_neg, 0.0],
                [lx_cc_neg, ly],
                Int(boundary_id_cc_el_neg))
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg, 0.0],
                [lx_cc_neg+lx_el_neg, ly],
                Int(boundary_id_el_sep_neg))
    if !reduced_membrane_model
        bfacemask!(grid,
                    [lx_cc_neg+lx_el_neg+lx_sep, 0.0],
                    [lx_cc_neg+lx_el_neg+lx_sep, ly],
                    Int(boundary_id_sep_el_pos))
    end
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, ly],
                Int(boundary_id_el_cc_pos))
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos+lx_cc_pos, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos+lx_cc_pos, ly],
                Int(boundary_id_cc_pos))

    # assign IDs to inflow and outflow boundaries
    bfacemask!(grid,
                [lx_cc_neg, 0.0],
                [lx_cc_neg+lx_el_neg, 0.0],
                Int(boundary_id_el_neg_inflow))
    bfacemask!(grid,
                [lx_cc_neg, ly],
                [lx_cc_neg+lx_el_neg, ly],
                Int(boundary_id_el_neg_outflow))
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, 0.0],
                Int(boundary_id_el_pos_inflow))
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep, ly],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, ly],
                Int(boundary_id_el_pos_outflow))

    # assign domain IDs to subdomains
    cellmask!(grid,
                [0.0, 0.0],
                [lx_cc_neg, ly],
                Int(domain_id_cc_neg))
    cellmask!(grid,
                [lx_cc_neg, 0.0],
                [lx_cc_neg+lx_el_neg, ly],
                Int(domain_id_el_neg))
    if !reduced_membrane_model
        cellmask!(grid,
                    [lx_cc_neg+lx_el_neg, 0.0],
                    [lx_cc_neg+lx_el_neg+lx_sep, ly],
                    Int(domain_id_sep))
    end
    cellmask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, ly],
                Int(domain_id_el_pos))
    cellmask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos+lx_cc_pos, ly],
                Int(domain_id_cc_pos))

    # generate subgrids over different subdomains (and unions of subdomains)
    subgrid_neg = subgrid(grid, [Int(domain_id_cc_neg), Int(domain_id_el_neg)])
    subgrid_el_neg = subgrid(grid, [Int(domain_id_el_neg)])
    if reduced_membrane_model
        subgrid_sep = subgrid(grid, [Int(boundary_id_el_sep_neg)]; boundary=true, project=false)
    else
        subgrid_sep = subgrid(grid, [Int(domain_id_sep)])
    end
    subgrid_pos = subgrid(grid, [Int(domain_id_el_pos), Int(domain_id_cc_pos)])
    subgrid_el_pos = subgrid(grid, [Int(domain_id_el_pos)])
    if !reduced_membrane_model
        subgrid_el = subgrid(grid,
                     [Int(domain_id_el_neg), Int(domain_id_sep), Int(domain_id_el_pos)])
        subgrids = (subgrid_neg = subgrid_neg,
                    subgrid_el_neg = subgrid_el_neg,
                    subgrid_sep = subgrid_sep,
                    subgrid_pos = subgrid_pos,
                    subgrid_el_pos = subgrid_el_pos,
                    subgrid_el = subgrid_el)
    else
        subgrid_el = subgrid(grid, [Int(domain_id_el_neg), Int(domain_id_el_pos)])
        subgrids = (subgrid_neg = subgrid_neg,
                    subgrid_el_neg = subgrid_el_neg,
                    subgrid_sep = subgrid_sep,
                    subgrid_pos = subgrid_pos,
                    subgrid_el_pos = subgrid_el_pos,
                    subgrid_el = subgrid_el)
    end

    return (grid=grid, subgrids=subgrids)
end
