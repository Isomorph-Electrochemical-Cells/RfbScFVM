function system_bcondition_1d!(f, u, bnode, data)

    system_internal_interface_flux!(f, u, bnode, data)

    ly_cell = data.geom.ly_cell

    v_neg = data.boundary.v_out_neg
    p_in_neg = data.boundary.p_in_neg
    p_neg = p_in_neg - data.el_neg.μ * ly_cell / (2 * data.el_neg.kₕ) * v_neg
    boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_neg]].p,
                        data.bnd[:cc_neg_el_neg], p_neg)

    v_pos = data.boundary.v_out_pos
    p_in_pos = data.boundary.p_in_pos
    p_pos = p_in_pos - data.el_pos.μ * ly_cell / (2 * data.el_pos.kₕ) * v_pos
    boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_pos]].p,
                        data.bnd[:el_pos_cc_pos], p_pos)

    # ϕₛ_l
    boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_neg]].ϕₛ,
                        data.bnd[:cc_neg_left], data.boundary.ϕₛ_neg[])
    # ϕₛ_r
    boundary_dirichlet!(f, u, bnode,  data.idx[data.dom[:el_pos]].ϕₛ,
                        data.bnd[:cc_pos_right], data.boundary.ϕₛ_pos[])

    if data.study.non_isothermal
        heat_exchange_boundary!(f, u, bnode, data)
    end
end

function system_reaction_1d!(f, u, node, data)

    system_reaction!(f, u, node, data)

    region = node.region
    if region == data.dom[:el_neg] || region == data.dom[:el_pos]
        temp = data.study.non_isothermal ? u[data.idx[region].temp] : 1.0
        PE0 = data.scaling_params.PE0
        ly = data.geom.ly_cell
        idx_c2u = data.idx[region].c


        if region == data.dom[:el_neg]
            vy = data.boundary.v_out_neg # FIXME
            c_in = data.boundary.species_neg
        elseif region == data.dom[:el_pos]
            vy = data.boundary.v_out_pos # FIXME
            c_in = data.boundary.species_pos
        end

        for idx in eachindex(idx_c2u)
            f[idx_c2u[idx]] -= 2 * PE0 * vy * (c_in[idx] - u[idx_c2u[idx]]) / ly
        end

        if data.study.non_isothermal
            temp_in = data.boundary.temp_amb
            cpᵥ = data.electrolyte.cpᵥ
            f[data.idx[region].temp] = -2 * PE0 * cpᵥ * vy * (temp_in - temp) / ly
        end
    end

end

function system_edgereaction_1d!(f, u, edge, data)
    # evaluate joule heating term
    if data.study.non_isothermal #FIXME: UNCOMMENT
        region = edge.region
        temp_l = data.idx[data.dom[:el_neg]].temp
        temp_r = data.idx[data.dom[:el_pos]].temp
        ϕₛ_l = data.idx[data.dom[:el_neg]].ϕₛ
        ϕₛ_r = data.idx[data.dom[:el_pos]].ϕₛ
        ϕₗ_l = data.idx[data.dom[:el_neg]].ϕₗ
        ϕₗ_r = data.idx[data.dom[:el_pos]].ϕₗ
        deff_factor = one(eltype(u)) # TODO
        if region == data.dom[:cc_neg]
            σₑ = data.cc_neg.σₑ
            s_ohmₛ = σₑ * (u[ϕₛ_l,1]-u[ϕₛ_l,2]) * (u[ϕₛ_l,1]-u[ϕₛ_l,2])
            f[temp_l] = -s_ohmₛ
        elseif region == data.dom[:cc_pos]
            σₑ = data.cc_pos.σₑ
            s_ohmₛ = σₑ * (u[ϕₛ_r,1]-u[ϕₛ_r,2]) * (u[ϕₛ_r,1]-u[ϕₛ_r,2])
            f[temp_r] = -s_ohmₛ
        elseif region == data.dom[:el_neg]
            σₗ = (σl_eff(u, region, data, deff_factor, 1) +
                 σl_eff(u, region, data, deff_factor, 2)) / 2
            s_ohmₗ = σₗ * (u[ϕₗ_l,1]-u[ϕₗ_l,2]) * (u[ϕₗ_l,1]-u[ϕₗ_l,2])
            σₑ = data.el_neg.σₑ
            s_ohmₛ = σₑ * (u[ϕₛ_l,1]-u[ϕₛ_l,2]) * (u[ϕₛ_l,1]-u[ϕₛ_l,2])
            f[temp_l] = -(s_ohmₗ +  s_ohmₛ)
        elseif edge.region == data.dom[:el_pos]
            σₗ = (σl_eff(u, region, data, deff_factor, 1) +
                 σl_eff(u, region, data, deff_factor, 2)) / 2
            s_ohmₗ = σₗ * (u[ϕₗ_r,1]-u[ϕₗ_r,2]) * (u[ϕₗ_r,1]-u[ϕₗ_r,2])
            σₑ = data.el_pos.σₑ
            s_ohmₛ = σₑ * (u[ϕₛ_r,1]-u[ϕₛ_r,2]) * (u[ϕₛ_r,1]-u[ϕₛ_r,2])
            f[temp_r] = -(s_ohmₗ + s_ohmₛ)
        end
    end
end

function system_breaction_1d!(f, u, bnode, data)
    f .= 0.0

    system_bcondition_1d!(f, u, bnode, data)

    joule_heating_separator!(f, u, bnode, data)

end

function generic_operator_1d!(f, u, sys)
    data = VoronoiFVM.data(sys)
    xcoords = sys.grid.components[ExtendableGrids.XCoordinates]
    num_nodes = length(xcoords)
    indices = data.indices

    idx_p_neg = data.idx[data.dom[:el_neg]].p
    idx_p_pos = data.idx[data.dom[:el_pos]].p

    lx_sep = data.geom.lx_sep #thickness of the separator
    lx_el_neg = data.geom.lx_el_neg
    lx_el_pos = data.geom.lx_el_pos
    idx_sep = Int(ceil(num_nodes/2))
    x_sep = xcoords[idx_sep]
    @assert data.geom.lx_cc_neg + data.geom.lx_el_neg ≈ x_sep

    vxₕ = -data.sep.kh / data.sep.μ * (u[indices[idx_p_pos, idx_sep]] -
                                      u[indices[idx_p_neg, idx_sep]]) / lx_sep

    idx_ϕₗ_neg = data.idx[data.dom[:el_neg]].ϕₗ
    idx_ϕₗ_pos = data.idx[data.dom[:el_pos]].ϕₗ
    vxₑ = data.sep.kϕ / data.sep.μ * data.sep.cf * data.sep.zf *
         (u[indices[idx_ϕₗ_pos, idx_sep]] - u[indices[idx_ϕₗ_neg, idx_sep]]) / lx_sep
    vxₑ /= (data.scaling_params.PE0)

    vx = vxₕ + vxₑ

    f .= 0.0
    for idx_node in 1:(num_nodes)
        idx_neg = indices[idx_p_neg, idx_node]
        if idx_neg > 0
            x_neg = xcoords[idx_node] - data.geom.lx_cc_neg
            if idx_node > 1 && idx_node < num_nodes
                dx = (xcoords[idx_node+1] - xcoords[idx_node-1])/2
            elseif idx_node == 1
                dx = (xcoords[idx_node+1] - xcoords[idx_node])/2
            else
                dx = (xcoords[idx_node] - xcoords[idx_node-1])/2
            end
            @assert x_neg >= 0.0
            f[idx_neg] -= dx * vx / (lx_el_neg)
        end
        idx_pos = indices[idx_p_pos, idx_node]
        if idx_pos > 0
            x_pos_ref = cell_thickness(data.geom) - data.geom.lx_cc_pos
            x_pos = xcoords[idx_node] + data.geom.lx_sep - x_pos_ref
            if idx_node > 1 && idx_node < num_nodes
                dx = (xcoords[idx_node+1] - xcoords[idx_node-1])/2
            elseif idx_node == 1
                dx = (xcoords[idx_node+1] - xcoords[idx_node])/2
            else
                dx = (xcoords[idx_node] - xcoords[idx_node-1])/2
            end
            @assert x_pos <= 0.0
            f[idx_pos] += dx * vx / (lx_el_pos)
        end
    end
end

function get_velocity(u, region, data)
    u[data.idx[region].p]
end

function integrated_volumetric_current_density(r, u, temp, region, data)
    nodes = data.quadrature.nodes
    weights = data.quadrature.weights
    ly = data.geom.ly_cell

    Δϕ = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
    idx_c2u = data.idx[region].c
    c_avg = u[idx_c2u]
    c_in = region == :el_neg ? data.boundary.species_neg.data : data.boundary.species_pos.data
    c_fun(y, idx) = c_in[idx] + 2 * y * (c_avg[idx] - c_in[idx]) / ly

    v = (data.boundary.v_out_neg + data.boundary.v_out_pos) / 2 #region == :el_neg ? data.boundary.v_out_neg : data.boundary.v_out_pos # FIXME

    integral = zero(eltype(u))
    for quad_idx in eachindex(nodes)
        n = nodes[quad_idx]
        w = weights[quad_idx]

        c_ox =  c_fun(n, r.idx_ox)
        c_red = c_fun(n, r.idx_red)

        η = Δϕ - ϕ_eq(r, idx -> c_fun(n, idx), temp) # overpotential
        iᵥ = volumetric_current_density(r, η, c_ox, c_red, v, temp, region, data)

        integral += iᵥ*w
    end
    return integral /= ly
end

function physics_1d(data)
    num_eqs = num_variables(data.idx)
    VoronoiFVM.Physics(num_species  = num_eqs,
                       storage      = system_storage!,
                       bstorage     = system_bstorage!,
                       flux         = system_flux!,
                       reaction     = system_reaction_1d!,
                       edgereaction = system_edgereaction_1d!,
                       breaction    = system_breaction_1d!,
                       generic      = generic_operator_1d!,
                       data         = data)
end
