
function σl_eff_1d(u, region, data, deff_factor, node=0)
    if region == data.dom[:el_neg] || region == data.dom[:el_pos]
        indices = data.idx[region].c
        temp = data.study.non_isothermal ? u[data.idx[region].temp] : one(eltype(u))
        species = data.electrolyte.species
        σl_eff_value = zero(eltype(u))
        for idx_species in eachindex(indices)
            d = species[row=idx_species, col=Key("diffusivity")]
            z = species[row=idx_species, col=Key("charge")]
            if node == 0
                c = u[indices[idx_species]]
            else
                c = u[indices[idx_species], node]
            end
            σl_eff_value += d * z^2 * c
        end
        # evaluate contribution of counter species
        counter = data.electrolyte.counter
        d = counter[col=Key("diffusivity")]
        z = counter[col=Key("charge")]
        @inline c = concentration_counter_species(u, region, data, node)
        σl_eff_value += d * z^2 * c
        return σl_eff_value * deff_factor / temp
    elseif region == data.dom[:sep]
        σₑ_sep = data.sep.σₑ
        if data.study.non_isothermal
            σₑ_sep += data.sep.∂σₑ∂temp * (u[data.idx[region].temp] - data.sep.temp_ref)
        end
        return σₑ_sep
    end
    return 0.0
end



function system_flux_1d!(f, u, edge, data)
    # Flux of Darcy's law
    if edge.region == data.dom[:el_neg]
        kₕ = data.el_neg.kₕ
        μ = data.el_neg.μ
        idx_p = data.idx[edge.region].p
        vh = kₕ/μ * (u[idx_p, 1] - u[idx_p, 2])
        f[idx_p] = vh
    elseif edge.region == data.dom[:el_pos]
        kₕ = data.el_pos.kₕ
        μ = data.el_pos.μ
        idx_p = data.idx[edge.region].p
        vh = kₕ/μ * (u[idx_p, 1] - u[idx_p, 2])
        f[idx_p] = vh
    else
        vh = 0
    end

    # Flux of the electrostatic potential ϕₛ
    if edge.region == data.dom[:cc_neg]
        σs_eff = data.cc_neg.σₑ
    elseif edge.region == data.dom[:cc_pos]
        σs_eff = data.cc_pos.σₑ
    elseif edge.region == data.dom[:el_neg]
        σs_eff = data.el_neg.σₑ
    elseif edge.region == data.dom[:el_pos]
        σs_eff = data.el_pos.σₑ
    end
    if edge.region == data.dom[:cc_neg] || edge.region == data.dom[:el_neg]
        idx_ϕₛ = data.idx[edge.region].ϕₛ
        f[idx_ϕₛ] = σs_eff * (u[idx_ϕₛ, 1] - u[idx_ϕₛ, 2])
    elseif edge.region == data.dom[:cc_pos] || edge.region == data.dom[:el_pos]
        idx_ϕₛ = data.idx[edge.region].ϕₛ
        f[idx_ϕₛ] = σs_eff * (u[idx_ϕₛ, 1] - u[idx_ϕₛ, 2])
    end
    PE0 = data.scaling_params.PE0

    # Flux of the electrostatic potential ϕₗ
    if edge.region == data.dom[:el_neg] || edge.region == data.dom[:el_pos]
        deff_factor = effective_diffusivity_factor(0.0, 0.0, edge.region, data)
        v_out = edge.region == data.dom[:el_neg] ? data.boundary.v_out_neg : data.boundary.v_out_pos

        dstar_factor = effective_diffusivity_factor(abs(v_out), 0.0, edge.region, data)
        σl_eff_avg = (σl_eff_1d(u, edge.region, data, deff_factor, 1) +
                      σl_eff_1d(u, edge.region, data, deff_factor, 2))/2
        idx_ϕₗ = data.idx[edge.region].ϕₗ
        f[idx_ϕₗ] = σl_eff_avg * (u[idx_ϕₗ, 1] - u[idx_ϕₗ, 2])
        indices = data.idx[edge.region].c
        species = data.electrolyte.species
        for idx_species in eachindex(indices)
            d = species[row=idx_species, col=Key("diffusivity")]
            deff = d * dstar_factor
            z = species[row=idx_species, col=Key("charge")]
            Δc = u[indices[idx_species], 1] - u[indices[idx_species], 2]
            f[idx_ϕₗ] += deff * z * Δc
        end
        counter = data.electrolyte.counter
        d = counter[col=Key("diffusivity")]
        deff = d * dstar_factor
        z = counter[col=Key("charge")]
        Δc = concentration_counter_species(u, edge.region, data, 1) -
                concentration_counter_species(u, edge.region, data, 2)
        f[idx_ϕₗ] += deff * z * Δc

        for idx in eachindex(indices)
            d = data.electrolyte.species[row=idx, col=Key("diffusivity")]
            deff_factor = effective_diffusivity_factor(0.0, 0.0, edge.region, data)
            v_out = edge.region == data.dom[:el_neg] ? data.boundary.v_out_neg : data.boundary.v_out_pos
            dstar_factor = effective_diffusivity_factor(abs(v_out), 0.0, edge.region, data)

            # Migration term
            migration = 0.0
            if data.study.migration
                if data.study.non_isothermal
                    idx_temp = data.idx[edge.region].temp
                    temp = (u[idx_temp, 1] + u[idx_temp, 2]) / 2
                else
                    temp = 1.0
                end

                idx_ϕₗ = data.idx[edge.region].ϕₗ
                ∂ϕ∂n =  (u[idx_ϕₗ, 1] - u[idx_ϕₗ, 2])

                z = data.electrolyte.species[row=idx, col=Key("charge")]
                migration = -z * (d * deff_factor) * ∂ϕ∂n / temp
            end

            idx_c = indices[idx]
            f[idx_c] = fvc_flux_Δx(u[idx_c, 1], u[idx_c, 2],
                                PE0*vh - migration, d * dstar_factor)
        end
    end

    if data.study.non_isothermal
        # Heat flux
        LE0 = data.scaling_params.LE0
        if edge.region == data.dom[:cc_neg]
            cpᵥ = data.cc_neg.cpᵥ
            λₜ = data.cc_neg.λₜ
        elseif edge.region == data.dom[:cc_pos]
            cpᵥ = data.cc_pos.cpᵥ
            λₜ = data.cc_pos.λₜ
        elseif edge.region == data.dom[:el_neg]
            εₗ = data.el_neg.εₗ
            cpᵥ = data.electrolyte.cpᵥ # εₗ * data.electrolyte.cpᵥ + (1-εₗ) * data.el_neg.cpᵥ #FIXME
            λₜ = data.electrolyte.λₜ # εₗ * data.electrolyte.λₜ + (1-εₗ) * data.el_neg.λₜ #FIXME
        elseif edge.region == data.dom[:el_pos]
            εₗ = data.el_pos.εₗ
            cpᵥ = data.electrolyte.cpᵥ #εₗ * data.electrolyte.cpᵥ + (1-εₗ) * data.el_pos.cpᵥ #FIXME
            λₜ = data.electrolyte.λₜ #εₗ * data.electrolyte.λₜ + (1-εₗ) * data.el_pos.λₜ #FIXME
        else
            @assert false
        end

        idx_temp = data.idx[edge.region].temp
        f[idx_temp] = fvc_flux_Δx(u[idx_temp, 1], u[idx_temp, 2],
                                  PE0 * cpᵥ * vh, LE0 * λₜ)
    end
end


function system_bcondition_1d!(f, u, bnode, data)

    system_internal_interface_flux!(f, u, bnode, data)

    # # p_l
    boundary_dirichlet!(f, u, bnode, 1, data.bnd[:el_neg_inflow],
                                        data.boundary.p_in_neg) # pressure inlet

    # p_r
    boundary_dirichlet!(f, u, bnode, 2, data.bnd[:el_pos_inflow],
                                        data.boundary.p_in_pos) # pressure inlet

    # ϕₛ_l
    boundary_dirichlet!(f, u, bnode, data.idx[data.dom[:el_neg]].ϕₛ,
                                              data.bnd[:cc_neg_left],
                                              data.boundary.ϕₛ_neg[])
    # ϕₛ_r
    boundary_dirichlet!(f, u, bnode,  data.idx[data.dom[:el_pos]].ϕₛ,
                                      data.bnd[:cc_pos_right],
                                        data.boundary.ϕₛ_pos[])

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
function system_bflux_1d!(f, u, edge, data)
    #f .= 0.0

    if edge.region == data.dom[:sep] && data.study.non_isothermal
        # Heat flux along separator interface
        LE0 = data.scaling_params.LE0
        λₜ = data.sep.λₜ
        # the y-component of the velocity is assumed to be negligible within the seprator

        temp_i = data.idx[data.dom[:sep]].temp
        f[temp_i] = fvc_flux_Δx(u[temp_i, 1], u[temp_i, 2], 0.0, LE0 * λₜ)
    end
end

# Note: The reaction term r appears on the left-hand-side of the transport eq., i.e.
#       ∂ₜs(u) + ∇⋅j + r(u) - f = 0
function system_reaction_1d!(f, u, node, data)
    f .= 0
    ϵL0 = data.scaling_params.ϵL0
    if node.region == data.dom[:el_neg]
        r = data.el_neg.reactions
        temp_l = data.idx[data.dom[:el_neg]].temp
        temp = data.study.non_isothermal ? u[temp_l] : 1.0

        Δϕ_l = u[data.idx[data.dom[:el_neg]].ϕₛ] - u[data.idx[data.dom[:el_neg]].ϕₗ]
        idx_ox_neg_l = data.idx[data.dom[:el_neg]].c[1] # TODO: Determine index of ox_neg_l in a more general way
        idx_red_neg_l = data.idx[data.dom[:el_neg]].c[2] # TODO: Determine index of ox_neg_l in a more general way
        η_l = Δϕ_l - ϕ_eq(r, u[idx_ox_neg_l], u[idx_red_neg_l], temp)

        @inline iᵥ = volumetric_current_density(r,
                                       η_l,
                                       u[idx_ox_neg_l], u[idx_red_neg_l],
                                       data.boundary.v_out_neg, #TODO: Use local velocity field
                                       temp,
                                       node.region, data)
        ϵL0_inv_squared = 1.0/(ϵL0*ϵL0)
        f[idx_ox_neg_l] = -ϵL0_inv_squared * r.ν_ox/r.ν_el * iᵥ
        f[idx_red_neg_l] = -ϵL0_inv_squared * r.ν_red/r.ν_el * iᵥ

        f[data.idx[data.dom[:el_neg]].ϕₛ] = ϵL0_inv_squared * iᵥ
        f[data.idx[data.dom[:el_neg]].ϕₗ] = -ϵL0_inv_squared * iᵥ

        # if iᵥ*η_l < -eps(typeof(iᵥ))
        #     @warn "iᵥ*η_r = $(iᵥ*η_l)"
        # end

        if data.study.non_isothermal
            f[temp_l] = -ϵL0_inv_squared * iᵥ * (η_l + u[temp_l]/r.ν_el * r.Δs)
        end

    elseif node.region == data.dom[:el_pos]
        r = data.el_pos.reactions
        temp_r = data.idx[data.dom[:el_pos]].temp
        temp = data.study.non_isothermal ? u[temp_r] : 1.0

        idx_ox_pos_r = data.idx[data.dom[:el_pos]].c[3] # TODO: Determine index of ox_neg_l in a more general way
        idx_red_pos_r = data.idx[data.dom[:el_pos]].c[4] # TODO: Determine index of ox_neg_l in a more general way

        Δϕ_r = u[data.idx[data.dom[:el_pos]].ϕₛ] - u[data.idx[data.dom[:el_pos]].ϕₗ]
        η_r = Δϕ_r - ϕ_eq(r, u[idx_ox_pos_r], u[idx_red_pos_r], temp)
        @inline iᵥ = volumetric_current_density(r,
                                       η_r,
                                       u[idx_ox_pos_r], u[idx_red_pos_r],
                                       data.boundary.v_out_pos, #TODO: Use local velocity field
                                       temp,
                                       node.region, data)

        ϵL0_inv_squared = 1.0/(ϵL0*ϵL0)
        f[idx_ox_pos_r] = -ϵL0_inv_squared * r.ν_ox/r.ν_el * iᵥ
        f[idx_red_pos_r] = -ϵL0_inv_squared * r.ν_red/r.ν_el * iᵥ

        f[data.idx[data.dom[:el_pos]].ϕₛ] = ϵL0_inv_squared * iᵥ
        f[data.idx[data.dom[:el_pos]].ϕₗ] = -ϵL0_inv_squared * iᵥ

        # if iᵥ*η_r < -eps(typeof(iᵥ))
        #     @warn "iᵥ*η_r = $(iᵥ*η_r)"
        # end

        if data.study.non_isothermal
            f[temp_r] = -ϵL0_inv_squared * iᵥ * (η_r + u[temp_r]/r.ν_el * r.Δs)
        end
    end
end


function system_storage_1d!(f, u, node, data)
    ϵL0 = data.scaling_params.ϵL0
    PE0 = data.scaling_params.PE0

    ϵ = 1e-10 # small value enforcing steady-state solution
    p_l = data.idx[data.dom[:el_neg]].p
    p_r = data.idx[data.dom[:el_pos]].p
    f[p_l] = u[p_l] * ϵ
    f[p_r] = u[p_r] * ϵ

    ϕₛ_l = data.idx[data.dom[:el_neg]].ϕₛ
    ϕₛ_r = data.idx[data.dom[:el_pos]].ϕₛ
    if node.region == data.dom[:cc_neg]
        f[ϕₛ_l] = u[ϕₛ_l] * ϵ
    elseif node.region == data.dom[:el_neg]
        f[ϕₛ_l] = u[ϕₛ_l] * ϵL0 * PE0 * data.el_neg.Cᵥ
    end

    if node.region == data.dom[:cc_pos]
        f[ϕₛ_r] = u[ϕₛ_r] * ϵ
    elseif node.region == data.dom[:el_pos]
        f[ϕₛ_r] = u[ϕₛ_r] * ϵL0 * PE0 * data.el_pos.Cᵥ
    end

    ϕₗ_l = data.idx[data.dom[:el_neg]].ϕₗ
    ϕₗ_r = data.idx[data.dom[:el_pos]].ϕₗ
    f[ϕₗ_l] = -u[ϕₗ_l] * ϵL0 * PE0 * data.el_neg.Cᵥ
    f[ϕₗ_r] = -u[ϕₗ_r] * ϵL0 * PE0 * data.el_pos.Cᵥ

    for idx in data.idx[data.dom[:el_neg]].c # species in negative half cell
        f[idx] = u[idx] * PE0 * data.el_neg.εₗ
    end

    for idx in data.idx[data.dom[:el_pos]].c # species in positive half cell
        f[idx] = u[idx] * PE0 * data.el_pos.εₗ
    end

    if data.study.non_isothermal
        temp_l = data.idx[data.dom[:el_neg]].temp
        temp_r = data.idx[data.dom[:el_pos]].temp

        f[temp_l] = u[temp_l] * PE0
        f[temp_r] = u[temp_r] * PE0
    end
end



function system_edgereaction_1d!(f, u, edge, data)
    # evaluate joule heating term
    if data.study.non_isothermal
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
        elseif region == data.dom[:domain_id_cc_pos]
            σₑ = data.cc_pos.σₑ
            s_ohmₛ = σₑ * (u[ϕₛ_r,1]-u[ϕₛ_r,2]) * (u[ϕₛ_r,1]-u[ϕₛ_r,2])
            f[temp_r] = -s_ohmₛ
        elseif region == data.dom[:domain_id_el_neg]
            σₗ = (σl_eff_1d(u, region, data, deff_factor, 1) +
                 σl_eff_1d(u, region, data, deff_factor, 2)) / 2
            s_ohmₗ = σₗ * (u[ϕₗ_l,1]-u[ϕₗ_l,2])*(u[ϕₗ_l,1]-u[ϕₗ_l,2])
            σₑ = data.el_neg.σₑ
            s_ohmₛ = σₑ * (u[ϕₛ_l,1]-u[ϕₛ_l,2]) * (u[ϕₛ_l,1]-u[ϕₛ_l,2])
            f[temp_l] = -(s_ohmₗ +  s_ohmₛ)
        elseif edge.region == data.dom[:domain_id_el_pos]
            σₗ = (σl_eff_1d(u, region, data, deff_factor, 1) +
                 σl_eff_1d(u, region, data, deff_factor, 2)) / 2
            s_ohmₗ = σₗ * (u[ϕₗ_r,1]-u[ϕₗ_r,2]) * (u[ϕₗ_r,1]-u[ϕₗ_r,2])
            σₑ = data.el_pos.σₑ
            s_ohmₛ = σₑ * (u[ϕₛ_r,1]-u[ϕₛ_r,2])*(u[ϕₛ_r,1]-u[ϕₛ_r,2])
            f[temp_r] = -(s_ohmₗ + s_ohmₛ)
        end
    end
end

function system_breaction_1d!(f, u, bnode, data)
    #alloc = @allocated begin
    @inline system_bcondition_1d!(f, u, bnode, data)
    #end
    #alloc > 0 && @show alloc

    if bnode.region == data.bnd[:sep]
        # joule heating over separator domain integrated along the through-plane direction
        if data.study.non_isothermal
            h = data.geom.lx_sep
            # evaluate Donnan potentials
            z_c = data.electrolyte.counter[col=Key("charge")]
            c_counter_neg = concentration_counter_species(u, data.dom[:el_neg], data, 0)
            c_counter_neg += eps(eltype(u))
            @assert c_counter_neg > 0
            c_counter_pos = concentration_counter_species(u, data.dom[:el_pos], data, 0)
            c_counter_pos += eps(eltype(u))
            @assert c_counter_pos > 0
            c_counter_sep = -data.sep.cf * data.sep.zf / z_c
            @assert c_counter_sep > 0
            ϕₗ_sep_neg = u[data.idx[data.dom[:el_neg]].ϕₗ] #+ 1/z_c * log(c_counter_neg / c_counter_sep) # FIXME
            ϕₗ_sep_pos = u[data.idx[data.dom[:el_pos]].ϕₗ] #+ 1/z_c * log(c_counter_pos / c_counter_sep) # FIXME

            σₑ_sep = data.sep.σₑ
            temp_i = data.idx[data.dom[:sep]].temp
            σₑ_sep += data.sep.∂σₑ∂temp * (u[temp_i] - data.sep.temp_ref)

            f[temp_i] -= σₑ_sep * ((ϕₗ_sep_pos - ϕₗ_sep_neg) / h)^2
        end
    end

end

function num_variables_1d(var)
    nvar = 0
    for domain in var
        nvar = max(nvar, domain.p, domain.ϕₛ, domain.ϕₗ, maximum(domain.c), domain.temp)
    end
    return nvar
end

function physics_1d(prms)
    num_eqs = num_variables_1d(prms.idx)
    VoronoiFVM.Physics(num_species  = num_eqs,
                       storage      = system_storage_1d!,
                       bstorage     = system_bstorage!,
                       flux         = system_flux_1d!,
                       bflux        = system_bflux_1d!,
                       reaction     = system_reaction_1d!,
                       edgereaction = system_edgereaction_1d!,
                       breaction    = system_breaction_1d!,
                       data         = prms)
end
