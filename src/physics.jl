# function velocity_edge_projection(u, edge, data) # returns velocity projected on edge
#     if edge.region == Int(domain_id_el_neg)
#         kₕ = data.el_neg.kₕ
#         μ = data.el_neg.μ
#         vh_edge = -kₕ/μ * (u[Int(p_l),2]-u[Int(p_l),1])
#     elseif edge.region == Int(domain_id_el_pos)
#         kₕ = data.el_pos.kₕ
#         μ = data.el_pos.μ
#         vh_edge = -kₕ/μ * (u[Int(p_r),2]-u[Int(p_r),1])
#     else
#         vh_edge = zero(eltype(u))
#     end
#     return vh_edge
# end

function darcy_outflow!(f, u, node, data)
    @infiltrate


    if node.region==Int(boundary_id_el_neg_outflow)
        @infiltrate
        f[data.var[:el_neg].p] = data.boundary.v_out_neg
    elseif node.region==Int(boundary_id_el_pos_outflow)
        f[data.var[:el_pos].p] = data.boundary.v_out_pos
    end
end

function species_outflow!(f, u, node, data)
    if node.region == boundary_id(:el_neg_outflow) ||
       node.region == boundary_id(:el_pos_outflow)

        PE0 = data.scaling_params.PE0
        if node.region == boundary_id(:el_neg_outflow)
            v_out = PE0*data.boundary.v_out_neg
            domain_id = :el_neg
        elseif node.region == boundary_id(:el_pos_outflow)
            v_out = PE0*data.boundary.v_out_pos
            domain_id = :el_pos
        end
        indices = data.var[domain_id].c
        for idx in indices
            f[idx] = v_out*u[idx]
        end
    end
end

function heat_outflow!(f, u, bnode, data)
    if bnode.region == boundary_id(:el_neg_outflow) ||
       bnode.region == boundary_id(:el_pos_outflow)

        PE0 = data.scaling_params.PE0
        if bnode.region == boundary_id(:el_neg_outflow)
            v_out = PE0 * data.boundary.v_out_neg
            cpᵥ = data.electrolyte.cpᵥ
            domain_sym = :el_neg
        elseif bnode.region == boundary_id(:el_pos_outflow)
            v_out = PE0 * data.boundary.v_out_pos
            cpᵥ = data.electrolyte.cpᵥ
            domain_sym = :el_pos
        end
        idx_temp = data.var[domain_sym].temp
        f[idx_temp] = v_out * cpᵥ * u[idx_temp]
    end
end


function heat_exchange_boundary!(f, u, bnode, data)
    LE0 = data.scaling_params.LE0
    if bnode.region == boundary_id(:cc_neg_left)
        idx_temp = data.var[:cc_neg].temp
        f[idx_temp] = LE0 * data.cc_neg.hₜ * (u[idx_temp]-data.boundary.temp_amb)
    end
    if bnode.region == Int(boundary_id_cc_pos)
        idx_temp = data.var[:cc_neg].temp
        f[idx_temp] = LE0 * data.cc_pos.hₜ * (u[idx_temp]-data.boundary.temp_amb)
    end
end


function concentration_counter_species(u, region, data, node)
    if region == domain_id(:el_neg) || region == domain_id(:el_pos)
        indices = data.var[Int(region)].c
        species = data.electrolyte.species
        c_counter = 0.0
        for idx_species in eachindex(indices)
            z = species[row=idx_species, col="charge"]
            if node == 0
                c = u[indices[idx_species]]
            else
                c = u[indices[idx_species], node]
            end
            c_counter -= z * c
        end
        c_counter /= species[row="counter", col="charge"]
        if c_counter <= 0
            @warn c_counter
        end
        return c_counter
    end
    return 0.0
end

function σl_eff(u, region, data, deff_factor, node=0)
    if region == domain_id(:el_neg) || region == domain_id(:el_pos)
        indices = data.var[Int(region)].c
        temp = data.study.non_isothermal ? u[data.var[Int(region)].temp] : 1.0
        species = data.electrolyte.species
        σl_eff_value = 0.0
        for idx_species in eachindex(indices)
            d = species[row=idx_species, col="diffusivity"]
            z = species[row=idx_species, col="charge"]
            if node == 0
                c = u[indices[idx_species]]
            else
                c = u[indices[idx_species], node]
            end
            σl_eff_value += d * z^2 * c
        end
        # evaluate contribution of counter species
        d = species[row="counter",col="diffusivity"]
        z = species[row="counter",col="charge"]
        c = concentration_counter_species(u, region, data, node)
        σl_eff_value += d * z^2 * c
        return σl_eff_value * deff_factor / temp
    elseif region == boundary_id(:sep)
        σₑ_sep = data.sep.σₑ
        if data.study.non_isothermal
            σₑ_sep += data.sep.∂σₑ∂temp * (u[data.var[region].temp] - data.sep.temp_ref)
        end
        return σₑ_sep
    end
    return 0.0
end

function system_flux!(f, u, edge, data)
    @infiltrate

    # Flux of Darcy's law
    if edge.region == domain_id(:el_neg)
        @infiltrate
        kₕ = data.el_neg.kₕ
        μ = data.el_neg.μ
        idx_p = data.var[Int(edge.region)].p
        vh = kₕ/μ * (u[idx_p, 1] - u[idx_p, 2])
        f[idx_p] = vh
    elseif edge.region == domain_id(:el_pos)
        kₕ = data.el_pos.kₕ
        μ = data.el_pos.μ
        idx_p = data.var[Int(edge.region)].p
        vh = kₕ/μ * (u[idx_p, 1] - u[idx_p, 2])
        f[idx_p] = vh
    else
        vh = 0
    end

    # Flux of the electrostatic potential ϕₛ
    if edge.region == domain_id(:cc_neg)
        σs_eff = data.cc_neg.σₑ
    elseif edge.region == domain_id(:cc_pos)
        σs_eff = data.cc_pos.σₑ
    elseif edge.region == domain_id(:el_neg)
        σs_eff = data.el_neg.σₑ
    elseif edge.region == domain_id(:el_pos)
        σs_eff = data.el_pos.σₑ
    end
    if edge.region == domain_id(:cc_neg) || edge.region == domain_id(:el_neg)
        idx_ϕₛ = data.var[Int(edge.region)].ϕₛ
        f[idx_ϕₛ] = σs_eff * (u[idx_ϕₛ, 1] - u[idx_ϕₛ, 2])
    elseif edge.region == domain_id(:cc_pos) || edge.region == domain_id(:el_pos)
        idx_ϕₛ = data.var[Int(edge.region)].ϕₛ
        f[idx_ϕₛ] = σs_eff * (u[idx_ϕₛ, 1] - u[idx_ϕₛ, 2])
    end
    PE0 = data.scaling_params.PE0

    # Flux of the electrostatic potential ϕₗ
    if edge.region == domain_id(:el_neg) || edge.region == domain_id(:el_pos)
        deff_factor = effective_diffusivity_factor(0.0, 0.0, edge.region, data)
        v_out = edge.region == domain_id(:el_neg) ? data.boundary.v_out_neg : data.boundary.v_out_pos

        dstar_factor = effective_diffusivity_factor(abs(v_out), 0.0, edge.region, data)
        σl_eff_avg = (σl_eff(u, edge.region, data, deff_factor, 1) +
                      σl_eff(u, edge.region, data, deff_factor, 2))/2
        idx_ϕₗ = data.var[Int(edge.region)].ϕₗ
        f[idx_ϕₗ] = σl_eff_avg * (u[idx_ϕₗ, 1] - u[idx_ϕₗ, 2])
        indices = data.var[Int(edge.region)].c
        species = data.electrolyte.species
        for idx_species in eachindex(indices)
            d = species[row=idx_species, col="diffusivity"]
            deff = d * dstar_factor
            z = species[row=idx_species, col="charge"]
            Δc = u[indices[idx_species], 1] - u[indices[idx_species], 2]
            f[idx_ϕₗ] += deff * z * Δc
        end
        d = species[row="counter", col="diffusivity"]
        deff = d * dstar_factor
        z = species[row="counter", col="charge"]
        Δc = concentration_counter_species(u, edge.region, data, 1) -
                concentration_counter_species(u, edge.region, data, 2)
        f[idx_ϕₗ] += deff * z * Δc

        for idx in eachindex(indices)
            idx_c = indices[idx]
            d = data.electrolyte.species[row=idx, col="diffusivity"]
            deff_factor = effective_diffusivity_factor(0.0, 0.0, edge.region, data)
            v_out = edge.region == domain_id(:el_neg) ? data.boundary.v_out_neg : data.boundary.v_out_pos
            dstar_factor = effective_diffusivity_factor(abs(v_out), 0.0, edge.region, data)

            # Migration term
            migration = 0.0
            if data.study.migration
                if data.study.non_isothermal
                    idx_temp = data.var[Int(edge.region)].temp
                    temp = (u[idx_temp, 1] + u[idx_temp, 2]) / 2
                else
                    temp = 1.0
                end

                idx_ϕₗ = data.var[Int(edge.region)].ϕₗ
                ∂ϕ∂n =  (u[idx_ϕₗ, 1] - u[idx_ϕₗ, 2])

                z = data.electrolyte.species[row=idx, col="charge"]
                migration = -z * (d * deff_factor) * ∂ϕ∂n / temp
            end

            f[idx_c] = fvc_flux_Δx(u[idx_c, 1], u[idx_c, 2],
                                PE0*vh - migration, d * dstar_factor)
        end
    end

    if data.study.non_isothermal
        # Heat flux
        LE0 = data.scaling_params.LE0
        if edge.region == domain_id(:cc_neg)
            cpᵥ = data.cc_neg.cpᵥ
            λₜ = data.cc_neg.λₜ
        elseif edge.region == domain_id(:cc_pos)
            cpᵥ = data.cc_pos.cpᵥ
            λₜ = data.cc_pos.λₜ
        elseif edge.region == domain_id(:el_neg)
            εₗ = data.el_neg.εₗ
            cpᵥ = data.electrolyte.cpᵥ # εₗ * data.electrolyte.cpᵥ + (1-εₗ) * data.el_neg.cpᵥ
            λₜ = data.electrolyte.λₜ # εₗ * data.electrolyte.λₜ + (1-εₗ) * data.el_neg.λₜ
        elseif edge.region == domain_id(:el_pos)
            εₗ = data.el_pos.εₗ
            cpᵥ = data.electrolyte.cpᵥ #εₗ * data.electrolyte.cpᵥ + (1-εₗ) * data.el_pos.cpᵥ
            λₜ = data.electrolyte.λₜ #εₗ * data.electrolyte.λₜ + (1-εₗ) * data.el_pos.λₜ
        else
            @assert false
        end

        idx_temp = data.var[Int(edge.region)].temp
        f[idx_temp] = fvc_flux_Δx(u[idx_temp, 1], u[idx_temp, 2],
                                  PE0 * cpᵥ * vh, LE0 * λₜ)
    end
end



# function system_internal_interface_flux!(f, u, bnode, data)
#     if bnode.region == boundary_id(:sep)
# 		h = data.geom.lx_sep #thickness of the separator

#         # flow velocity due to a hydrostatic pressure difference across the separator
#         idx_p_l = data.var[:el_neg].p
#         idx_p_r = data.var[:el_pos].p
# 		vₕ = -data.sep.kh / data.sep.μ * (u[idx_p_r] - u[idx_p_l]) / h
#         # electrokinetic effect
#         idx_ϕₗ_l = data.var[:el_neg].ϕₗ
#         idx_ϕₗ_r = data.var[:el_pos].ϕₗ
#         vₑ = data.sep.kϕ * data.sep.cf * data.sep.zf * (u[idx_ϕₗ_r] - u[idx_ϕₗ_l])
#         vₑ /= (data.sep.μ * h * data.scaling_params.PE0)
#         v = vₕ + vₑ

#         f[idx_p_l] = v # left volumetric flux
#         f[idx_p_r] = -v # right volumetric flux

#         # evaluate Donnan potentials
#         z_c = data.electrolyte.species[row="counter", col="charge"]
#         ϵ = eps(eltype(u))
#         c_c_el_l = concentration_counter_species(u, domain_id(:el_neg), data, 0) + ϵ
#         @infiltrate
#         @assert c_c_el_l > 0
#         c_c_el_r = concentration_counter_species(u, domain_id(:el_pos), data, 0) + ϵ
#         @assert c_c_el_r > 0
#         c_c_sep = -data.sep.cf * data.sep.zf / z_c
#         @assert c_c_sep > 0
#         ϕₗ_sep_l = u[idx_ϕₗ_l] + 1/z_c * log(c_c_el_l / c_c_sep)
#         ϕₗ_sep_r = u[idx_ϕₗ_r] + 1/z_c * log(c_c_el_r / c_c_sep)

#         σₑ_sep = data.sep.σₑ
#         if data.study.non_isothermal
#             @show data.var
#             idx_temp = data.var[:sep].temp

#             σₑ_sep += data.sep.∂σₑ∂temp * (u[idx_temp] - data.sep.temp_ref)
#         end
#         current_flux_sep = -σₑ_sep * (ϕₗ_sep_r - ϕₗ_sep_l) / h

#         f[idx_ϕₗ_l] = current_flux_sep
#         f[idx_ϕₗ_r] = -current_flux_sep

#         if data.study.non_isothermal
#             PE0 = data.scaling_params.PE0
#             LE0 = data.scaling_params.LE0
#             cpᵥ_sep = data.sep.cpᵥ
#             λₜ_sep = data.sep.λₜ

#             idx_temp_l = data.var[:el_neg].temp
#             idx_temp_sep = data.var[:sep].temp
#             idx_temp_r = data.var[:el_pos].temp
#             heat_flux_left = fvc_flux_Δx(u[idx_temp_l], u[idx_temp_sep],
#                                            PE0*cpᵥ_sep*v,
#                                            LE0*λₜ_sep,
#                                            h/2)

#             heat_flux_right = fvc_flux_Δx(u[idx_temp_sep], u[idx_temp_r],
#                                             PE0*cpᵥ_sep*v,
#                                             LE0*λₜ_sep,
#                                             h/2)

#             f[idx_temp_l] = heat_flux_left / (h/2)
#             f[idx_temp_r] = -heat_flux_right / (h/2)
#             f[idx_temp_sep] = -(heat_flux_left - heat_flux_right) / (h/2)
#         end
#     end
# end

function system_internal_interface_flux!(f, u, bnode, data)
    if bnode.region == Int(boundary_id_el_sep_neg)
        @infiltrate
		h = data.geom.lx_sep #thickness of the separator

        # flow velocity due to a hydrostatic pressure difference across the separator
		vₕ = -data.sep.kh / data.sep.μ * (u[Int(p_r)] - u[Int(p_l)]) / h
        # electrokinetic effect
        vₑ = data.sep.kϕ / data.sep.μ * data.sep.cf * data.sep.zf * (u[Int(ϕₗ_r)] - u[Int(ϕₗ_l)]) / h # TODO: CHECK SIGN!!!!
        vₑ /= (data.scaling_params.PE0)
        v = vₕ + vₑ

        f[Int(p_l)] = v # left volumetric flux #FIXME
        f[Int(p_r)] = -v # right volumetric flux # FIXME

        # evaluate Donnan potentials
        z_c = data.electrolyte.species[row="counter", col="charge"]
        c_c_el_l = concentration_counter_species(u, Int(domain_id_el_neg), data, 0) + eps(eltype(u))
        @assert c_c_el_l > 0
        c_c_el_r = concentration_counter_species(u, Int(domain_id_el_pos), data, 0) + eps(eltype(u))
        @assert c_c_el_r > 0
        c_c_sep = -data.sep.cf * data.sep.zf / z_c
        @assert c_c_sep > 0
        ϕₗ_sep_l = u[Int(ϕₗ_l)] + 1/z_c * log(c_c_el_l / c_c_sep)
        ϕₗ_sep_r = u[Int(ϕₗ_r)] + 1/z_c * log(c_c_el_r / c_c_sep)

        σₑ_sep = data.sep.σₑ
        if data.study.non_isothermal
            σₑ_sep += data.sep.∂σₑ∂temp * (u[Int(temp_i)] - data.sep.temp_ref)
        end
        current_flux_sep = -σₑ_sep * (ϕₗ_sep_r - ϕₗ_sep_l) / h

        f[Int(ϕₗ_l)] = current_flux_sep
        f[Int(ϕₗ_r)] = -current_flux_sep

        if data.study.non_isothermal

            PE0 = data.scaling_params.PE0
            LE0 = data.scaling_params.LE0
            cpᵥ_sep = data.sep.cpᵥ
            λₜ_sep = data.sep.λₜ

            heat_flux_left = fvc_flux_Δx(u[Int(temp_l)], u[Int(temp_i)],
                                           PE0*cpᵥ_sep*v,
                                           LE0*λₜ_sep,
                                           h/2)

            heat_flux_right = fvc_flux_Δx(u[Int(temp_i)], u[Int(temp_r)],
                                            PE0*cpᵥ_sep*v,
                                            LE0*λₜ_sep,
                                            h/2)

            f[Int(temp_l)] = heat_flux_left / (h/2)
            f[Int(temp_r)] = -heat_flux_right / (h/2)
            f[Int(temp_i)] = -(heat_flux_left - heat_flux_right) / (h/2)
        end
    end
end



function system_breaction!(f, u, bnode, data)
    system_bcondition!(f, u, bnode, data)

    if bnode.region == Int(boundary_id_el_sep_neg)
        # joule heating over separator domain integrated along the through-plane direction
        if data.study.non_isothermal
            h = data.geom.lx_sep

            # evaluate Donnan potentials
            z_c = data.electrolyte.species[row="counter", col="charge"]
            c_c_el_l = concentration_counter_species(u, Int(domain_id_el_neg), data, 0) + eps(eltype(u))
            @assert c_c_el_l > 0
            c_c_el_r = concentration_counter_species(u, Int(domain_id_el_pos), data, 0) + eps(eltype(u))
            @assert c_c_el_r > 0
            c_c_sep = -data.sep.cf * data.sep.zf / z_c
            @assert c_c_sep > 0
            ϕₗ_sep_l = u[Int(ϕₗ_l)] #+ 1/z_c * log(c_c_el_l / c_c_sep)
            ϕₗ_sep_r = u[Int(ϕₗ_r)] #+ 1/z_c * log(c_c_el_r / c_c_sep)

            σₑ_sep = data.sep.σₑ
            if data.study.non_isothermal
                σₑ_sep += data.sep.∂σₑ∂temp * (u[Int(temp_i)] - data.sep.temp_ref)
            end
            f[Int(temp_i)] -= σₑ_sep * ((ϕₗ_sep_r - ϕₗ_sep_l) / h)^2
        end
    end

end


function system_bcondition!(f, u, bnode, data)
    @infiltrate

    if data.pressure_boundary_type == :p_in_v_out
        darcy_outflow!(f, u, bnode, data)
    end
    species_outflow!(f, u, bnode, data)

    system_internal_interface_flux!(f, u, bnode, data)
    if data.study.non_isothermal
        heat_outflow!(f, u, bnode, data)
    end

    # # p_l
    boundary_dirichlet!(f, u, bnode, 1, Int(boundary_id_el_neg_inflow),
                                        data.boundary.p_in_neg) # pressure inlet
    if data.pressure_boundary_type == :p_in_p_out
        @assert false "Not yet implemented"
        # boundary_dirichlet!(f, u, bnode, 1, Int(boundary_id_el_neg_outflow),
        #                                 data.boundary.p_out_neg) # pressure outlet
    end

    # p_r
    boundary_dirichlet!(f, u, bnode, 2, Int(boundary_id_el_pos_inflow),
                                        data.boundary.p_in_pos) # pressure inlet
    if data.pressure_boundary_type == :p_in_p_out
        @assert false "Not yet implemented"
        # boundary_dirichlet!(f, u, bnode, 2, Int(boundary_id_el_pos_outflow),
        #                                 data.boundary.p_r_out) # pressure outlet
    end

    # ϕₛ_l
    boundary_dirichlet!(f, u, bnode, Int(ϕₛ_l), Int(boundary_id_cc_neg),
                                        data.boundary.ϕₛ_neg)
    # ϕₛ_r
    boundary_dirichlet!(f, u, bnode, Int(ϕₛ_r), Int(boundary_id_cc_pos),
                                        data.boundary.ϕₛ_pos)

    # c_ox_neg_l
    boundary_dirichlet!(f, u, bnode, Int(c_ox_neg_l), Int(boundary_id_el_neg_inflow),
                        data.boundary.species_neg[row="ox_neg"][col="concentration"])
    # c_ox_neg_r
    boundary_dirichlet!(f, u, bnode, Int(c_ox_neg_r), Int(boundary_id_el_pos_inflow),
                        data.boundary.species_pos[row="ox_neg"][col="concentration"])

    # c_red_neg_l
    boundary_dirichlet!(f, u, bnode, Int(c_red_neg_l), Int(boundary_id_el_neg_inflow),
                        data.boundary.species_neg[row="red_neg"][col="concentration"])
    # c_red_neg_r
    boundary_dirichlet!(f, u, bnode, Int(c_red_neg_r), Int(boundary_id_el_pos_inflow),
                        data.boundary.species_pos[row="red_neg"][col="concentration"])

    # c_ox_pos_l
    boundary_dirichlet!(f, u, bnode, Int(c_ox_pos_l), Int(boundary_id_el_neg_inflow),
                        data.boundary.species_neg[row="ox_pos"][col="concentration"])
    # c_ox_pos_r
    boundary_dirichlet!(f, u, bnode, Int(c_ox_pos_r), Int(boundary_id_el_pos_inflow),
                        data.boundary.species_pos[row="ox_pos"][col="concentration"])

    # c_red_pos_l
    boundary_dirichlet!(f, u, bnode, Int(c_red_pos_l), Int(boundary_id_el_neg_inflow),
                        data.boundary.species_neg[row="red_pos"][col="concentration"])
    # c_red_pos_r
    boundary_dirichlet!(f, u, bnode, Int(c_red_pos_r), Int(boundary_id_el_pos_inflow),
                        data.boundary.species_pos[row="red_pos"][col="concentration"])

    if data.study.non_isothermal
        # temp_l
        boundary_dirichlet!(f, u, bnode, Int(temp_l), Int(boundary_id_el_neg_inflow),
                                            data.boundary.temp_amb)
        # temp_r
        boundary_dirichlet!(f, u, bnode, Int(temp_r), Int(boundary_id_el_pos_inflow),
                                            data.boundary.temp_amb)
        heat_exchange_boundary!(f, u, bnode, data)
    end
end

# Boundary flux (i.e. transport processes along boundaries)
function system_bflux!(f, u, edge, data)
    #f .= 0.0

    if edge.region == Int(boundary_id_el_sep_neg) && data.study.non_isothermal
        # Heat flux along separator interface
        LE0 = data.scaling_params.LE0
        λₜ = data.sep.λₜ
        # the y-component of the velocity is assumed to be negligible within the seprator

        f[Int(temp_i)] = fvc_flux_Δx(u[Int(temp_i), 1], u[Int(temp_i), 2],
                                    0.0, LE0 * λₜ)
    end
end

function regularize_non_negative_quantity(u)
    return u <= zero(u) ? zero(u) : u
end

# Note: The reaction term r appears on the left-hand-side of the transport eq., i.e.
#       ∂ₜs(u) + ∇⋅j + r(u) - f = 0
function system_reaction!(f, u, node, data)
    f .= 0
    ϵL0 = data.scaling_params.ϵL0
    if node.region == Int(domain_id_el_neg)
        r = data.el_neg.reactions
        temp = data.study.non_isothermal ? u[Int(temp_l)] : 1.0

        Δϕ_l = u[Int(ϕₛ_l)] - u[Int(ϕₗ_l)]
        η_l = Δϕ_l - ϕ_eq(r, u[Int(c_ox_neg_l)], u[Int(c_red_neg_l)], temp)
        iᵥ = volumetric_current_density(r,
                                       η_l,
                                       u[Int(c_ox_neg_l)], u[Int(c_red_neg_l)],
                                       data.boundary.v_out_neg, #TODO: Use local velocity field
                                       temp,
                                       node.region, data)

        ϵL0_inv_squared = 1.0/(ϵL0*ϵL0)
        f[Int(c_ox_neg_l)] = -ϵL0_inv_squared * r.ν_ox/r.ν_el * iᵥ
        f[Int(c_red_neg_l)] = -ϵL0_inv_squared * r.ν_red/r.ν_el * iᵥ

        f[Int(ϕₛ_l)] = ϵL0_inv_squared * iᵥ
        f[Int(ϕₗ_l)] = -ϵL0_inv_squared * iᵥ

        # if iᵥ*η_l < -eps(typeof(iᵥ))
        #     @warn "iᵥ*η_r = $(iᵥ*η_l)"
        # end

        if data.study.non_isothermal
            f[Int(temp_l)] = -ϵL0_inv_squared * iᵥ * (η_l + u[Int(temp_l)]/r.ν_el * r.Δs)
        end

    elseif node.region == Int(domain_id_el_pos)
        r = data.el_pos.reactions
        temp = data.study.non_isothermal ? u[Int(temp_r)] : 1.0

        Δϕ_r = u[Int(ϕₛ_r)] - u[Int(ϕₗ_r)]
        η_r = Δϕ_r - ϕ_eq(r, u[Int(c_ox_pos_r)], u[Int(c_red_pos_r)], temp)
        iᵥ = volumetric_current_density(r,
                                       η_r,
                                       u[Int(c_ox_pos_r)], u[Int(c_red_pos_r)],
                                       data.boundary.v_out_pos, #TODO: Use local velocity field
                                       temp,
                                       node.region, data)

        ϵL0_inv_squared = 1.0/(ϵL0*ϵL0)
        f[Int(c_ox_pos_r)] = -ϵL0_inv_squared * r.ν_ox/r.ν_el * iᵥ
        f[Int(c_red_pos_r)] = -ϵL0_inv_squared * r.ν_red/r.ν_el * iᵥ

        f[Int(ϕₛ_r)] = ϵL0_inv_squared * iᵥ
        f[Int(ϕₗ_r)] = -ϵL0_inv_squared * iᵥ

        # if iᵥ*η_r < -eps(typeof(iᵥ))
        #     @warn "iᵥ*η_r = $(iᵥ*η_r)"
        # end

        if data.study.non_isothermal
            f[Int(temp_r)] = -ϵL0_inv_squared * iᵥ * (η_r + u[Int(temp_r)]/r.ν_el * r.Δs)
        end
    end
end


function system_storage!(f, u, node, data)
    ϵL0 = data.scaling_params.ϵL0
    PE0 = data.scaling_params.PE0

    ϵ = 1e-10 # small value enforcing steady-state solution
    f[Int(p_l)] = u[Int(p_l)] * ϵ
    f[Int(p_r)] = u[Int(p_r)] * ϵ

    if node.region == Int(domain_id_cc_neg)
        f[Int(ϕₛ_l)] = u[Int(ϕₛ_l)] * ϵ
    elseif node.region == Int(domain_id_el_neg)
        f[Int(ϕₛ_l)] = u[Int(ϕₛ_l)] * ϵL0 * PE0 * data.el_neg.Cᵥ
    end

    if node.region == Int(domain_id_cc_pos)
        f[Int(ϕₛ_r)] = u[Int(ϕₛ_r)] * ϵ
    elseif node.region == Int(domain_id_el_neg)
        f[Int(ϕₛ_r)] = u[Int(ϕₛ_r)] * ϵL0 * PE0 * data.el_pos.Cᵥ
    end

    f[Int(ϕₗ_l)] = -u[Int(ϕₗ_l)] * ϵL0 * PE0 * data.el_neg.Cᵥ
    f[Int(ϕₗ_r)] = -u[Int(ϕₗ_r)] * ϵL0 * PE0 * data.el_pos.Cᵥ

    for idx in 7:2:13 # species in negative half cell
        f[idx] = u[idx] * PE0 * data.el_neg.εₗ
    end

    for idx in 8:2:14 # species in positive half cell
        f[idx] = u[idx] * PE0 * data.el_pos.εₗ
    end

    if data.study.non_isothermal
        f[15] = u[15] * PE0
        f[16] = u[16] * PE0
    end
end

function system_bstorage!(f, u, node, data)
    if node.region == Int(boundary_id_el_sep_neg) && data.study.non_isothermal
        PE0 = data.scaling_params.PE0
        f[Int(temp_i)] = u[Int(temp_i)] * PE0
    end
end

function system_edgereaction!(f, u, edge, data)
    # evaluate joule heating term
    if data.study.non_isothermal
        region = edge.region
        deff_factor = 1.0 # TODO
        if region == Int(domain_id_cc_neg)
            σₑ = data.cc_neg.σₑ
            s_ohmₛ = σₑ * (u[Int(ϕₛ_l),1]-u[Int(ϕₛ_l),2]) * (u[Int(ϕₛ_l),1]-u[Int(ϕₛ_l),2])
            f[Int(temp_l)] = -s_ohmₛ
        elseif region == Int(domain_id_cc_pos)
            σₑ = data.cc_pos.σₑ
            s_ohmₛ = σₑ * (u[Int(ϕₛ_r),1]-u[Int(ϕₛ_r),2]) * (u[Int(ϕₛ_r),1]-u[Int(ϕₛ_r),2])
            f[Int(temp_r)] = -s_ohmₛ
        elseif region == Int(domain_id_el_neg)
            σₗ = (σl_eff(u, region, data, deff_factor, 1) +
                 σl_eff(u, region, data, deff_factor, 2)) / 2
            s_ohmₗ = σₗ * (u[Int(ϕₗ_l),1]-u[Int(ϕₗ_l),2])*(u[Int(ϕₗ_l),1]-u[Int(ϕₗ_l),2])
            σₑ = data.el_neg.σₑ
            s_ohmₛ = σₑ * (u[Int(ϕₛ_l),1]-u[Int(ϕₛ_l),2]) * (u[Int(ϕₛ_l),1]-u[Int(ϕₛ_l),2])
            f[Int(temp_l)] = -(s_ohmₗ +  s_ohmₛ)
        elseif edge.region == Int(domain_id_el_pos)
            σₗ = (σl_eff(u, region, data, deff_factor, 1) +
                 σl_eff(u, region, data, deff_factor, 2)) / 2
            s_ohmₗ = σₗ * (u[Int(ϕₗ_r),1]-u[Int(ϕₗ_r),2]) * (u[Int(ϕₗ_r),1]-u[Int(ϕₗ_r),2])
            σₑ = data.el_pos.σₑ
            s_ohmₛ = σₑ * (u[Int(ϕₛ_r),1]-u[Int(ϕₛ_r),2])*(u[Int(ϕₛ_r),1]-u[Int(ϕₛ_r),2])
            f[Int(temp_r)] = -(s_ohmₗ + s_ohmₛ)
        end
    end
end

function physics_2d(physics_data)
    num_eqs = length(instances(Variables))
    VoronoiFVM.Physics(num_species  = num_eqs,
                       storage      = system_storage!,
                       bstorage     = system_bstorage!,
                       flux         = system_flux!,
                       bflux        = system_bflux!,
                       reaction     = system_reaction!,
                       edgereaction = system_edgereaction!,
                       breaction    = system_breaction!,
                       data         = physics_data)
end
