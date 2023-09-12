function darcy_outflow!(f, u, node, data)
    if node.region == data.bnd[:el_neg_outflow]
        f[data.idx[data.dom[:el_neg]].p] = data.boundary.v_out_neg
    elseif node.region == data.bnd[:el_pos_outflow]
        f[data.idx[data.dom[:el_pos]].p] = data.boundary.v_out_pos
    end
end

function species_outflow!(f, u, node, data)
    PE0 = data.scaling_params.PE0
    if node.region == data.bnd[:el_neg_outflow]
        v_out = PE0 * data.boundary.v_out_neg
        species_outflow_region!(f, u, v_out, data.idx[data.dom[:el_neg]].c)
    elseif node.region == data.bnd[:el_pos_outflow]
        v_out = PE0 * data.boundary.v_out_pos
        species_outflow_region!(f, u, v_out, data.idx[data.dom[:el_pos]].c)
    end
end

function species_outflow_region!(f, u, v_out, idx_species)
    for idx in idx_species
        f[idx] = v_out*u[idx]
    end
end

function heat_outflow!(f, u, bnode, data)
    @assert data.study.non_isothermal

    if bnode.region == data.bnd[:el_neg_outflow] ||
       bnode.region == data.bnd[:el_pos_outflow]

        PE0 = data.scaling_params.PE0
        if bnode.region == data.bnd[:el_neg_outflow]
            v_out = PE0 * data.boundary.v_out_neg
            cpᵥ = data.electrolyte.cpᵥ
            idx_temp = data.idx[data.dom[:el_neg]].temp
        else
            v_out = PE0 * data.boundary.v_out_pos
            cpᵥ = data.electrolyte.cpᵥ
            idx_temp = data.idx[data.dom[:el_pos]].temp
        end
        f[idx_temp] = v_out * cpᵥ * u[idx_temp]
    end
end

function system_bcondition_2d!(f, u, bnode, data)

    darcy_outflow!(f, u, bnode, data)

    species_outflow!(f, u, bnode, data)

    system_internal_interface_flux!(f, u, bnode, data)
    if data.study.non_isothermal
        heat_outflow!(f, u, bnode, data)
    end

    # # p_l
    boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_neg]].p,
                        data.bnd[:el_neg_inflow], data.boundary.p_in_neg) # pressure inlet

    # p_r
    boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_pos]].p,
                        data.bnd[:el_pos_inflow], data.boundary.p_in_pos) # pressure inlet

    # ϕₛ_l
    boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_neg]].ϕₛ,
                                     data.bnd[:cc_neg_left],
                                     data.boundary.ϕₛ_neg[])
    # ϕₛ_r
    boundary_dirichlet!(f, u, bnode,  data.idx[data.dom[:el_pos]].ϕₛ,
                                      data.bnd[:cc_pos_right],
                                      data.boundary.ϕₛ_pos[])

    for idx in eachindex(data.boundary.species_neg)
        idx_var = data.idx[data.dom[:el_neg]].c[idx]
        bdr_val = data.boundary.species_neg[idx]
        boundary_dirichlet!(f, u, bnode, idx_var, data.bnd[:el_neg_inflow], bdr_val)
    end

    for idx in eachindex(data.boundary.species_pos)
        idx_var = data.idx[data.dom[:el_pos]].c[idx]
        bdr_val = data.boundary.species_pos[idx]
        boundary_dirichlet!(f, u, bnode, idx_var, data.bnd[:el_pos_inflow], bdr_val)
    end

    if data.study.non_isothermal
        # temp_l
        boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_neg]].temp,
            data.bnd[:el_neg_inflow], data.boundary.temp_amb)
        # temp_r
        boundary_dirichlet!(f, u, bnode,  data.idx[data.dom[:el_pos]].temp,
            data.bnd[:el_pos_inflow], data.boundary.temp_amb)
        heat_exchange_boundary!(f, u, bnode, data)
    end
end

# Boundary flux (i.e. transport processes along boundaries)
function system_bflux_2d!(f, u, edge, data)
    #f .= 0.0

    if edge.region == data.dom[:sep] && data.study.non_isothermal
        # Heat flux along separator interface
        LE0 = data.scaling_params.LE0
        λₜ = data.sep.λₜ

        temp_i = data.idx[data.dom[:sep]].temp
        f[temp_i] = fvc_flux_Δx(u[temp_i, 1], u[temp_i, 2], 0.0, LE0 * λₜ)
    end
end


function system_edgereaction_2d!(f, u, edge, data)
    # evaluate joule heating term
    if data.study.non_isothermal
        region = edge.region
        idx_temp = data.idx[region].temp
        idx_ϕₛ = data.idx[region].ϕₛ
        idx_ϕₗ = data.idx[region].ϕₗ
        deff_factor = one(eltype(u)) # TODO
        if region == data.dom[:cc_neg] || region == data.dom[:cc_pos]
            σₑ = region == data.dom[:cc_neg] ? data.cc_neg.σₑ : data.cc_pos.σₑ
            s_ohmₛ = σₑ * (u[idx_ϕₛ,1]-u[idx_ϕₛ,2]) * (u[idx_ϕₛ,1]-u[idx_ϕₛ,2])
            f[idx_temp] = -s_ohmₛ
        elseif region == data.dom[:el_neg] || region == data.dom[:el_pos]
            σₑ = region == data.dom[:el_neg] ? data.el_neg.σₑ : data.el_pos.σₑ
            σₗ = (σl_eff(u, region, data, deff_factor, 1) +
                 σl_eff(u, region, data, deff_factor, 2)) / 2
            s_ohmₗ = σₗ * (u[idx_ϕₗ,1]-u[idx_ϕₗ,2])*(u[idx_ϕₗ,1]-u[idx_ϕₗ,2])
            s_ohmₛ = σₑ * (u[idx_ϕₛ,1]-u[idx_ϕₛ,2]) * (u[idx_ϕₛ,1]-u[idx_ϕₛ,2])
            f[idx_temp] = -(s_ohmₗ +  s_ohmₛ)
        end
    end
end

function system_breaction_2d!(f, u, bnode, data)
    system_bcondition_2d!(f, u, bnode, data)

    # if bnode.region == data.bnd[:sep]
    #     # joule heating over separator domain integrated along the through-plane direction
    #     if data.study.non_isothermal # FIXME: USE GENERAL IMPLEMENTATION OF DONNAN PORENTIALS
    #         h = data.geom.lx_sep
    #         # evaluate Donnan potentials
    #         z_c = data.electrolyte.counter[col=Key("charge")]
    #         c_counter_neg = concentration_eliminated_species(u, data.dom[:el_neg], data, 0)
    #         c_counter_neg += eps(eltype(u))
    #         @assert c_counter_neg > 0
    #         c_counter_pos = concentration_eliminated_species(u, data.dom[:el_pos], data, 0)
    #         c_counter_pos += eps(eltype(u))
    #         @assert c_counter_pos > 0
    #         c_counter_sep = -data.sep.cf * data.sep.zf / z_c
    #         @assert c_counter_sep > 0
    #         ϕₗ_sep_neg = u[data.idx[data.dom[:el_neg]].ϕₗ] #+ 1/z_c * log(c_counter_neg / c_counter_sep) # FIXME
    #         ϕₗ_sep_pos = u[data.idx[data.dom[:el_pos]].ϕₗ] #+ 1/z_c * log(c_counter_pos / c_counter_sep) # FIXME

    #         σₑ_sep = data.sep.σₑ
    #         idx_temp_sep = data.idx[data.dom[:sep]].temp
    #         σₑ_sep += data.sep.∂σₑ∂temp * (u[idx_temp_sep] - data.sep.temp_ref) # FIXME: MAKE THIS CONSISTENT WITH THE HANDLING IN THE SEPARATOR

    #         f[idx_temp_sep] -= σₑ_sep * ((ϕₗ_sep_pos - ϕₗ_sep_neg) / h)^2
    #     end
    # end

    joule_heating_separator!(f, u, bnode, data)

end


function physics_2d(data)
    num_eqs = num_variables(data.idx)
    VoronoiFVM.Physics(num_species  = num_eqs,
                       storage      = system_storage!,
                       bstorage     = system_bstorage!,
                       flux         = system_flux!,
                       bflux        = system_bflux_2d!,
                       reaction     = system_reaction!,
                       edgereaction = system_edgereaction_2d!,
                       breaction    = system_breaction_2d!,
                       data         = data)
end
