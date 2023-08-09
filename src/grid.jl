"""
    generate_grid_nodes(interval, h)

Generate grid nodes for a one-dimensional domain interval. The array h specifies the
    spatial spacings at uniformly distributed locations. The spacings in between are
    calculated according to a geometric sequence.

# Arguments
- `interval`: two element array containing the end points of the interval
- `h`: array with spatial spacings

# Returns
- array containing the grid nodes
"""
function generate_grid_nodes(interval, h)
    num_subdomain_nodes = length(h)
    subdomain_nodes = collect(range(interval..., length=num_subdomain_nodes))
    nodes = []
    for idx_subdomain in 1:num_subdomain_nodes-1
        nodes_subdomain = geomspace(subdomain_nodes[idx_subdomain],
                                    subdomain_nodes[idx_subdomain+1],
                                    h[idx_subdomain], h[idx_subdomain+1])
        nodes = isempty(nodes) ? nodes_subdomain : glue(nodes, nodes_subdomain)
    end
    return nodes
end

function create_grid_1d(cell_geometry::FlowCellGeometry{T},
    mesh::Mesh1D{T}, dom_ids, bnd_ids) where {T<:AbstractFloat}

    lx_cc_neg = cell_geometry.lx_cc_neg
    lx_el_neg = cell_geometry.lx_el_neg
    lx_sep = 0.0
    lx_el_pos = cell_geometry.lx_el_pos
    lx_cc_pos = cell_geometry.lx_cc_pos
    ly = cell_geometry.ly_cell

    x_left = 0.0
    x_cc_el_neg = lx_cc_neg
    x_el_sep_neg = lx_cc_neg + lx_el_neg
    x_sep_el_pos = lx_cc_neg + lx_el_neg + lx_sep
    x_el_cc_pos = lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos
    x_right = lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos + lx_cc_pos

    domain_cc_neg = (x_left, x_cc_el_neg)
    domain_el_neg = (x_cc_el_neg, x_el_sep_neg)
    domain_sep = (x_el_sep_neg, x_sep_el_pos)
    domain_el_pos = (x_sep_el_pos, x_el_cc_pos)
    domain_cc_pos = (x_el_cc_pos, x_right)
    domain_x = (x_left, x_right)


    nodes_cc_neg = generate_grid_nodes(domain_cc_neg, mesh.hx_cc_neg)
    nodes_el_neg = generate_grid_nodes(domain_el_neg, mesh.hx_el_neg)
    nodes_el_pos = generate_grid_nodes(domain_el_pos, mesh.hx_el_pos)
    nodes_cc_pos = generate_grid_nodes(domain_cc_pos, mesh.hx_cc_pos)

    all_nodes_x = glue(nodes_cc_neg, nodes_el_neg)
    all_nodes_x = glue(all_nodes_x, nodes_el_pos)
    all_nodes_x = glue(all_nodes_x, nodes_cc_pos)

    # generate grid
    grid = simplexgrid(all_nodes_x)

    # assign IDs to boundaries
    bfacemask!(grid, [0.0], [0.0], bnd_ids[:cc_neg_left])
    bfacemask!(grid, [lx_cc_neg], [lx_cc_neg], bnd_ids[:cc_neg_el_neg])
    bfacemask!(grid, [lx_cc_neg + lx_el_neg], [lx_cc_neg + lx_el_neg], bnd_ids[:sep])
    bfacemask!(grid, [lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos],
                     [lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos],
                      bnd_ids[:el_pos_cc_pos])
    bfacemask!(grid, [lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos + lx_cc_pos],
                     [lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos + lx_cc_pos],
                      bnd_ids[:cc_pos_right])

                      # assign domain IDs to subdomains
    cellmask!(grid, [0.0], [lx_cc_neg], dom_ids[:cc_neg])
    cellmask!(grid, [lx_cc_neg], [lx_cc_neg + lx_el_neg], dom_ids[:el_neg])
    cellmask!(grid, [lx_cc_neg + lx_el_neg + lx_sep],
                    [lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos],
                    dom_ids[:el_pos])
    cellmask!(grid, [lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos],
                    [lx_cc_neg + lx_el_neg + lx_sep + lx_el_pos + lx_cc_pos],
                    dom_ids[:cc_pos])

    # generate subgrids over different subdomains (and unions of subdomains)
    subgrid_neg = subgrid(grid, [dom_ids[:cc_neg], dom_ids[:el_neg]])
    subgrid_el_neg = subgrid(grid, [dom_ids[:el_neg]])
    subgrid_sep = subgrid(grid, [dom_ids[:sep]]; boundary=true, project=false)
    subgrid_pos = subgrid(grid, [dom_ids[:el_pos], dom_ids[:cc_pos]])
    subgrid_el_pos = subgrid(grid, [dom_ids[:el_pos]])
    subgrid_el = subgrid(grid, [dom_ids[:el_neg], dom_ids[:el_pos]])
    subgrids = (subgrid_neg=subgrid_neg,
        subgrid_el_neg=subgrid_el_neg,
        subgrid_sep=subgrid_sep,
        subgrid_pos=subgrid_pos,
        subgrid_el_pos=subgrid_el_pos,
        subgrid_el=subgrid_el)

    return (grid=grid, subgrids=subgrids)
end



function create_grid_2d(cell_geometry::FlowCellGeometry{T},
                        mesh::Mesh2D{T}, dom_ids, bnd_ids) where {T<:AbstractFloat}

    lx_cc_neg = cell_geometry.lx_cc_neg
    lx_el_neg = cell_geometry.lx_el_neg
    lx_sep = 0.0
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

    nodes_cc_neg = generate_grid_nodes(domain_cc_neg, mesh.hx_cc_neg)
    nodes_el_neg = generate_grid_nodes(domain_el_neg, mesh.hx_el_neg)
    nodes_el_pos = generate_grid_nodes(domain_el_pos, mesh.hx_el_pos)
    nodes_cc_pos = generate_grid_nodes(domain_cc_pos, mesh.hx_cc_pos)
    nodes_y = generate_grid_nodes(domain_y, mesh.hy_cell)

    all_nodes_x = glue(nodes_cc_neg, nodes_el_neg)
    all_nodes_x = glue(all_nodes_x, nodes_el_pos)
    all_nodes_x = glue(all_nodes_x, nodes_cc_pos)

    # generate grid
    grid = simplexgrid(all_nodes_x, nodes_y)

    # assign IDs to boundaries
    bfacemask!(grid,
                [0.0, 0.0],
                [0.0, ly],
                bnd_ids[:cc_neg_left])
    bfacemask!(grid,
                [lx_cc_neg, 0.0],
                [lx_cc_neg, ly],
                bnd_ids[:cc_neg_el_neg])
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg, 0.0],
                [lx_cc_neg+lx_el_neg, ly],
                bnd_ids[:sep])
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, ly],
                bnd_ids[:el_pos_cc_pos])
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos+lx_cc_pos, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos+lx_cc_pos, ly],
                bnd_ids[:cc_pos_right])

    # assign IDs to inflow and outflow boundaries
    bfacemask!(grid,
                [lx_cc_neg, 0.0],
                [lx_cc_neg+lx_el_neg, 0.0],
                bnd_ids[:el_neg_inflow])
    bfacemask!(grid,
                [lx_cc_neg, ly],
                [lx_cc_neg+lx_el_neg, ly],
                bnd_ids[:el_neg_outflow])
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, 0.0],
                bnd_ids[:el_pos_inflow])
    bfacemask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep, ly],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, ly],
                bnd_ids[:el_pos_outflow])

    # assign domain IDs to subdomains
    cellmask!(grid,
                [0.0, 0.0],
                [lx_cc_neg, ly],
                dom_ids[:cc_neg])
    cellmask!(grid,
                [lx_cc_neg, 0.0],
                [lx_cc_neg+lx_el_neg, ly],
                dom_ids[:el_neg])
    cellmask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, ly],
                dom_ids[:el_pos])
    cellmask!(grid,
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos, 0.0],
                [lx_cc_neg+lx_el_neg+lx_sep+lx_el_pos+lx_cc_pos, ly],
                dom_ids[:cc_pos])

    # generate subgrids over different subdomains (and unions of subdomains)
    subgrid_neg = subgrid(grid, [dom_ids[:cc_neg], dom_ids[:el_neg]])
    subgrid_el_neg = subgrid(grid, [ dom_ids[:el_neg]])
    subgrid_sep = subgrid(grid, [dom_ids[:sep]]; boundary=true, project=false)
    subgrid_pos = subgrid(grid, [dom_ids[:el_pos], dom_ids[:cc_pos]])
    subgrid_el_pos = subgrid(grid, [dom_ids[:el_pos]])
    subgrid_el = subgrid(grid, [dom_ids[:el_neg], dom_ids[:el_pos]])
    subgrids = (subgrid_neg = subgrid_neg,
                subgrid_el_neg = subgrid_el_neg,
                subgrid_sep = subgrid_sep,
                subgrid_pos = subgrid_pos,
                subgrid_el_pos = subgrid_el_pos,
                subgrid_el = subgrid_el)

    return (grid=grid, subgrids=subgrids)
end
