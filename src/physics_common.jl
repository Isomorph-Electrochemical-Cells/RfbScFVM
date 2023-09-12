function heat_exchange_boundary!(f, u, bnode, data)
    LE0 = data.scaling_params.LE0
    if bnode.region == data.bnd[:cc_neg_left]
        idx_temp = data.idx[data.dom[:cc_neg]].temp
        f[idx_temp] = LE0 * data.cc_neg.hₜ * (u[idx_temp]-data.boundary.temp_amb)
    end
    if bnode.region == data.bnd[:cc_pos_right]
        idx_temp = data.idx[data.dom[:cc_pos]].temp
        f[idx_temp] = LE0 * data.cc_pos.hₜ * (u[idx_temp]-data.boundary.temp_amb)
    end
end

function σl_eff(u, region, data, deff_factor, node=0)
    @assert region == data.dom[:el_neg] || region == data.dom[:el_pos]

    indices = data.idx[region].c
    temp = data.study.non_isothermal ? u[data.idx[region].temp] :
                                       one(eltype(u)) * data.boundary.temp_amb
    species = data.electrolyte.species
    num_species = count_all_species(species)
    σl_eff_value = zero(eltype(u))
    for idx_species in 1:num_species
        diff = effective_diffusion(species[row=idx_species], temp) #FIXME: Consider impact of porous geometry on diffusion!
        z = species[row=idx_species, col=Key("charge")]

        if idx_species == num_species
            c = concentration_eliminated_species(u, region, data, node)
        elseif node == 0
            c = u[indices[idx_species]]
        else
            c = u[indices[idx_species], node]
        end
        σl_eff_value += diff * z^2 * c
    end
    return σl_eff_value * deff_factor / temp
end

function σl_eff(c_all_species::T, deff_factor, temp, data, region) where {T<:AbstractVector}
    species = data.electrolyte.species
    if region == data.dom[:sep]
        species_diffusion = data.sep.species
    else
        species_diffusion = species
    end
    num_species = count_all_species(species)
    σl_eff_value = zero(eltype(c_all_species))
    for idx_species in 1:num_species
        diff = effective_diffusion(species_diffusion[row=idx_species], temp) #FIXME: Consider impact of porous geometry on diffusion!
        z = species[row=idx_species, col=Key("charge")]
        σl_eff_value += diff * z^2 * c_all_species[idx_species]
    end
    return σl_eff_value * deff_factor / temp
end


function concentration_eliminated_species(u, region, data, node)
    if region == data.dom[:el_neg] || region == data.dom[:el_pos]
        indices = data.idx[region].c
        c_counter = zero(eltype(u))
        for idx_species in eachindex(indices)
            z = data.electrolyte.species[row=idx_species, col=Key("charge")]
            if node == 0
                c = u[indices[idx_species]]
            else
                c = u[indices[idx_species], node]
            end
            if c < -1e-5#-SQRT_EPS
                @warn "negative species concentration"
                @show region == data.dom[:el_neg]
                @show region == data.dom[:el_pos]
                @show indices[idx_species]
                @show Float64(c)
            end
            c_counter -= z * c
        end
        c_counter /= data.electrolyte.species[row=end, col=Key("charge")]
        if c_counter < -SQRT_EPS
            @warn "negative species concentration"
        end
        return c_counter
    end
    return zero(eltype(u))
end


function concentration_eliminated_species(c_free_species, data)
    c_counter = zero(eltype(c_free_species))
    for idx_species in eachindex(c_free_species)
        z = data.electrolyte.species[row=idx_species, col=Key("charge")]
        c_counter -= z * c_free_species[idx_species]
    end
    c_counter /= data.electrolyte.species[row=end, col=Key("charge")]
    if c_counter < -SQRT_EPS
        @warn "negative species concentration"
    end
    return max(c_counter, zero(c_counter))
end



function system_internal_interface_flux!(f, u, bnode, data)
    if bnode.region == data.dom[:sep]
		h = data.geom.lx_sep #thickness of the separator

        # flow velocity due to a hydrostatic pressure difference across the separator
        p_l = data.idx[data.dom[:el_neg]].p
        p_r = data.idx[data.dom[:el_pos]].p
		vₕ = -data.sep.kh / data.sep.μ * (u[p_r] - u[p_l]) / h
        # electrokinetic effect
        ϕₗ_l = data.idx[data.dom[:el_neg]].ϕₗ
        ϕₗ_r = data.idx[data.dom[:el_pos]].ϕₗ
        vₑ = data.sep.kϕ / data.sep.μ * data.sep.cf * data.sep.zf * (u[ϕₗ_r] - u[ϕₗ_l]) / h # TODO: CHECK SIGN!!!
        vₑ /= (data.scaling_params.PE0)
        v = vₕ + vₑ # FIXME

        f[p_l] = v  # left volumetric flux
        f[p_r] = -v # right volumetric flux

        z = data.electrolyte.species[col=Key("charge")]
        cf = data.sep.cf
        zf = data.sep.zf

        if data.study.non_isothermal
            temp_neg = u[data.idx[data.dom[:el_neg]].temp]
            temp_pos = u[data.idx[data.dom[:el_pos]].temp]
            temp_i = u[data.idx[data.dom[:sep]].temp]
        else
            temp_neg = one(eltype(u)) * data.boundary.temp_amb
            temp_pos = one(eltype(u)) * data.boundary.temp_amb
            temp_i = one(eltype(u)) * data.boundary.temp_amb
        end

        cc_neg = concentration_eliminated_species(u, data.dom[:el_neg], data, 0)
        c_neg = vcat(u[data.idx[data.dom[:el_neg]].c], cc_neg)
        ϕₗ_sep_neg = donnan_equlibrium_potential(u[ϕₗ_l], c_neg, z, cf, zf, temp_neg)

        cc_pos = concentration_eliminated_species(u, data.dom[:el_pos], data, 0)
        c_pos = vcat(u[data.idx[data.dom[:el_pos]].c], cc_pos)
        ϕₗ_sep_pos = donnan_equlibrium_potential(u[ϕₗ_r], c_pos, z, cf, zf, temp_pos)

        # diffusion-driven fluxes through separator
        num_species = count_all_species(data.electrolyte.species)
        c_sep_l_all = @MVector zeros(num_species)
        donnan_equilibrium_concentrations!(c_sep_l_all, u[ϕₗ_l], ϕₗ_sep_neg, c_neg, z, temp_i)
        c_sep_r_all = @MVector zeros(num_species)
        donnan_equilibrium_concentrations!(c_sep_r_all, u[ϕₗ_r], ϕₗ_sep_pos, c_pos, z, temp_i)

        PE0 = data.scaling_params.PE0
        σₑ_sep_eval = zero(eltype(u))
        diffusion_current = zero(eltype(u))

        species_sep = data.sep.species
        species_el = data.electrolyte.species
        for idx_species in 1:num_species
            diff = effective_diffusion(species_sep[row=idx_species], temp_i)
            z = species_el[row=idx_species, col=Key("charge")]
            c_upwind = v >= 0 ? c_neg[idx_species] : c_pos[idx_species]

            diffusion_flux = -diff * (c_sep_r_all[idx_species] - c_sep_l_all[idx_species]) / h

            ∇ϕₗ_sep = (ϕₗ_sep_pos - ϕₗ_sep_neg) / h
            c_sep_upwind = ∇ϕₗ_sep <= 0 ? c_sep_l_all[idx_species] : c_sep_r_all[idx_species]
            migration_flux = -z * c_sep_upwind * diff / temp_i * ∇ϕₗ_sep

            convective_flux = v * c_upwind * PE0

            species_flux = diffusion_flux #+ convective_flux + migration_flux  # FIXME!!!

            σₑ_sep_eval += z^2 * c_upwind * diff
            diffusion_current += z * diffusion_flux

            if idx_species < num_species # species is a variable of the PDE system
                idx_species_neg = data.idx[data.dom[:el_neg]].c[idx_species]
                idx_species_pos = data.idx[data.dom[:el_pos]].c[idx_species]
                f[idx_species_neg] += species_flux
                f[idx_species_pos] -= species_flux
            end
        end

        current_flux_sep = -(σₑ_sep_eval / temp_i) * (ϕₗ_sep_pos - ϕₗ_sep_neg) / h +
                            diffusion_current #- PE0 * v * zf * cf # FIXME!!!

        f[ϕₗ_l] = current_flux_sep
        f[ϕₗ_r] = -current_flux_sep

        if data.study.non_isothermal
            PE0 = data.scaling_params.PE0
            LE0 = data.scaling_params.LE0
            cpᵥ_sep = data.sep.cpᵥ
            λₜ_sep = data.sep.λₜ

            temp_l = data.idx[data.dom[:el_neg]].temp
            temp_i = data.idx[data.dom[:sep]].temp
            temp_r = data.idx[data.dom[:el_pos]].temp
            heat_flux_left = fvc_flux_Δx(u[temp_l], u[temp_i],
                                           PE0*cpᵥ_sep*v,
                                           LE0*λₜ_sep,
                                           h/2)

            heat_flux_right = fvc_flux_Δx(u[temp_i], u[temp_r],
                                            PE0*cpᵥ_sep*v,
                                            LE0*λₜ_sep,
                                            h/2)

            f[temp_l] = heat_flux_left / (h/2)
            f[temp_r] = -heat_flux_right / (h/2)
            f[temp_i] = -(heat_flux_left - heat_flux_right) / (h/2)
        end
    end
end


# Note: The reaction term r appears on the left-hand-side of the transport eq., i.e.
#       ∂ₜs(u) + ∇⋅j + r(u) - f = 0
function system_reaction!(f, u, node, data)
    f .= 0
    ϵL0 = data.scaling_params.ϵL0
    region = node.region
    if region == data.dom[:el_neg]
        ϵL0_inv_squared = 1.0/(ϵL0*ϵL0)
        idx_c2u = data.idx[region].c
        temp_l = data.idx[region].temp
        temp = data.study.non_isothermal ? u[temp_l] : 1.0

        vec_r = data.el_neg.reactions
        for r in vec_r
            if occursin("1d", data.discr.spatial_discr)
                iᵥ = integrated_volumetric_current_density(r, u, temp, region, data) # 1D
            else
                v = data.boundary.v_out_neg  # FIXME
                c_ox =  u[idx_c2u[r.idx_ox]]
                c_red = u[idx_c2u[r.idx_red]]
                Δϕ = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
                η = Δϕ - ϕ_eq(r, idx -> u[idx_c2u[idx]], temp) # overpotential
                iᵥ = volumetric_current_density(r, η, c_ox, c_red, v, temp, region, data) # 2D
            end

            idx_ox_neg_l = data.idx[region].c[r.idx_ox]
            idx_red_neg_l = data.idx[region].c[r.idx_red]
            f[idx_ox_neg_l] += -ϵL0_inv_squared * r.ν_coeff[r.idx_ox]/r.ν_el * iᵥ
            f[idx_red_neg_l] += -ϵL0_inv_squared * r.ν_coeff[r.idx_red]/r.ν_el * iᵥ

            f[data.idx[region].ϕₛ] += ϵL0_inv_squared * iᵥ
            f[data.idx[region].ϕₗ] += -ϵL0_inv_squared * iᵥ

            if data.study.non_isothermal
                Δϕ_l = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
                idx_c2u = data.idx[region].c
                η_l = Δϕ_l - ϕ_eq(r, idx -> u[idx_c2u[idx]], temp)
                f[temp_l] -= ϵL0_inv_squared * iᵥ * (η_l + u[temp_l]/r.ν_el * r.Δs)
            end
        end

        # homogeneous reactions
        idx_c2u = data.idx[region].c
        reactions = data.electrolyte.hom_reactions
        mat_a = reactions.mat_a
        mat_c = reactions.mat_c
        num_species, num_reactions = size(mat_c)
        @assert num_species == length(idx_c2u)

        for idx_reaction in 1:num_reactions
            reactant_product = one(eltype(u))
            for idx_species in 1:num_species
                a::Int64 = mat_a[idx_reaction, idx_species]
                if a > 0
                    reactant_product *= u[idx_c2u[idx_species]]^a
                end
            end
            for idx_species in 1:num_species
                f[idx_c2u[idx_species]] -= mat_c[idx_species, idx_reaction] * reactant_product
            end
        end

    elseif region == data.dom[:el_pos]
        ϵL0_inv_squared = 1.0/(ϵL0*ϵL0)
        temp_r = data.idx[region].temp
        temp = data.study.non_isothermal ? u[temp_r] : 1.0
        idx_c2u = data.idx[region].c

        vec_r = data.el_pos.reactions
        for r in vec_r
            if occursin("1d", data.discr.spatial_discr)
                iᵥ = integrated_volumetric_current_density(r, u, temp, region, data) # 1D
            else
                v = data.boundary.v_out_pos  # FIXME
                c_ox =  u[idx_c2u[r.idx_ox]]
                c_red = u[idx_c2u[r.idx_red]]
                Δϕ = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
                η = Δϕ - ϕ_eq(r, idx -> u[idx_c2u[idx]], temp) # overpotential
                iᵥ = volumetric_current_density(r, η, c_ox, c_red, v, temp, region, data) # 2D
            end

            idx_ox_pos_r = data.idx[region].c[r.idx_ox]
            idx_red_pos_r = data.idx[region].c[r.idx_red]
            f[idx_ox_pos_r] = -ϵL0_inv_squared * r.ν_coeff[r.idx_ox]/r.ν_el * iᵥ
            f[idx_red_pos_r] = -ϵL0_inv_squared * r.ν_coeff[r.idx_red]/r.ν_el * iᵥ

            f[data.idx[region].ϕₛ] += ϵL0_inv_squared * iᵥ
            f[data.idx[region].ϕₗ] += -ϵL0_inv_squared * iᵥ

            if data.study.non_isothermal
                Δϕ_r = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
                idx_c2u = data.idx[region].c
                η_r = Δϕ_r - ϕ_eq(r, idx -> u[idx_c2u[idx]], temp)
                f[temp_r] -= ϵL0_inv_squared * iᵥ * (η_r + u[temp_r]/r.ν_el * r.Δs)
            end
        end

        # homogeneous reactions
        idx_c2u = data.idx[region].c
        reactions = data.electrolyte.hom_reactions
        mat_a = reactions.mat_a
        mat_c = reactions.mat_c
        num_species, num_reactions = size(mat_c)
        @assert num_species == length(idx_c2u)

        for idx_reaction in 1:num_reactions
            reactant_product = one(eltype(u))
            for idx_species in 1:num_species
                a::Int64 = mat_a[idx_reaction, idx_species]
                if a > 0
                    reactant_product *= u[idx_c2u[idx_species]]^a
                end
            end
            for idx_species in 1:num_species
                f[idx_c2u[idx_species]] -= mat_c[idx_species, idx_reaction] * reactant_product
            end
        end
    end
end


function system_flux!(f, u, edge, data)
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
    idx_ϕₛ = data.idx[edge.region].ϕₛ
    f[idx_ϕₛ] = σs_eff * (u[idx_ϕₛ, 1] - u[idx_ϕₛ, 2])

    PE0 = data.scaling_params.PE0

    if data.study.non_isothermal
        idx_temp = data.idx[edge.region].temp
        temp = (u[idx_temp, 1] + u[idx_temp, 2]) / 2
    else
        temp = data.boundary.temp_amb
    end

    # Flux of the electrostatic potential ϕₗ
    if edge.region == data.dom[:el_neg] || edge.region == data.dom[:el_pos]
        deff_factor = effective_diffusivity_factor(0.0, 0.0, edge.region, data)
        v_out = edge.region == data.dom[:el_neg] ? data.boundary.v_out_neg : data.boundary.v_out_pos

        dstar_factor = effective_diffusivity_factor(abs(v_out), 0.0, edge.region, data)
        σl_eff_avg = (σl_eff(u, edge.region, data, deff_factor, 1) +
                      σl_eff(u, edge.region, data, deff_factor, 2))/2

        idx_ϕₗ = data.idx[edge.region].ϕₗ
        f[idx_ϕₗ] = σl_eff_avg * (u[idx_ϕₗ, 1] - u[idx_ϕₗ, 2])
        indices = data.idx[edge.region].c
        species = data.electrolyte.species
        num_species = count_all_species(species)
        for idx_species in 1:num_species
            diff0 = species[row=idx_species, col=Key("reference_diffusivity")]
            temp0 = species[row=idx_species, col=Key("reference_temperature")]
            α = species[row=idx_species, col=Key("temperature_coefficient")]
            deff = diff0 * (1.0 + α * (temp - temp0)) * dstar_factor
            z = species[row=idx_species, col=Key("charge")]
            if idx_species < num_species
                Δc = u[indices[idx_species], 1] - u[indices[idx_species], 2]
            else
                Δc = concentration_eliminated_species(u, edge.region, data, 1) -
                        concentration_eliminated_species(u, edge.region, data, 2)
            end
            f[idx_ϕₗ] += deff * z * Δc
        end

        for idx in eachindex(indices)
            diff0 = species[row=idx, col=Key("reference_diffusivity")]
            temp0 = species[row=idx, col=Key("reference_temperature")]
            α = species[row=idx, col=Key("temperature_coefficient")]
            d = diff0 * (1.0 + α * (temp - temp0))

            deff_factor = effective_diffusivity_factor(0.0, 0.0, edge.region, data)
            v_out = edge.region == data.dom[:el_neg] ? data.boundary.v_out_neg : data.boundary.v_out_pos
            dstar_factor = effective_diffusivity_factor(abs(v_out), 0.0, edge.region, data)

            # Migration term
            migration = 0.0
            # if data.study.migration
                idx_ϕₗ = data.idx[edge.region].ϕₗ
                ∂ϕ∂n =  (u[idx_ϕₗ, 1] - u[idx_ϕₗ, 2])
                z = data.electrolyte.species[row=idx, col=Key("charge")]
                migration = -z * (d * deff_factor) * ∂ϕ∂n / temp
            #end

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
            cpᵥ = data.electrolyte.cpᵥ
            λₜ = εₗ * data.electrolyte.λₜ + (1-εₗ) * data.el_neg.λₜ
        elseif edge.region == data.dom[:el_pos]
            εₗ = data.el_pos.εₗ
            cpᵥ = data.electrolyte.cpᵥ
            λₜ = εₗ * data.electrolyte.λₜ + (1-εₗ) * data.el_pos.λₜ
        else
            @assert false
        end

        idx_temp = data.idx[edge.region].temp
        f[idx_temp] = fvc_flux_Δx(u[idx_temp, 1], u[idx_temp, 2], PE0 * cpᵥ * vh, LE0 * λₜ)
    end
end

function system_storage!(f, u, node, data)
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

function joule_heating_separator!(f, u, bnode, data)
    # joule heating over separator domain integrated along the through-plane direction
    if bnode.region == data.bnd[:sep] && data.study.non_isothermal
        lx_sep = data.geom.lx_sep

        # flow velocity due to a hydrostatic pressure difference across the separator
        p_l = data.idx[data.dom[:el_neg]].p
        p_r = data.idx[data.dom[:el_pos]].p
		vₕ = -data.sep.kh / data.sep.μ * (u[p_r] - u[p_l]) / lx_sep
        # electrokinetic effect
        ϕₗ_l = data.idx[data.dom[:el_neg]].ϕₗ
        ϕₗ_r = data.idx[data.dom[:el_pos]].ϕₗ
        vₑ = data.sep.kϕ / data.sep.μ * data.sep.cf * data.sep.zf * (u[ϕₗ_r] - u[ϕₗ_l]) / lx_sep
        vₑ /= (data.scaling_params.PE0)
        v = vₕ + vₑ

        # evaluate Donnan potentials
        z = data.electrolyte.species[col=Key("charge")]
        cf = data.sep.cf
        zf = data.sep.zf

        cc_neg = concentration_eliminated_species(u, data.dom[:el_neg], data, 0)
        c_neg = vcat(u[data.idx[data.dom[:el_neg]].c], cc_neg)
        temp_neg = u[data.idx[data.dom[:el_neg]].temp]
        ϕₗ_neg = data.idx[data.dom[:el_neg]].ϕₗ
        ϕₗ_sep_neg = donnan_equlibrium_potential(u[ϕₗ_neg], c_neg, z, cf, zf, temp_neg)

        cc_pos = concentration_eliminated_species(u, data.dom[:el_pos], data, 0)
        c_pos = vcat(u[data.idx[data.dom[:el_pos]].c], cc_pos)
        temp_pos = u[data.idx[data.dom[:el_pos]].temp]
        ϕₗ_pos = data.idx[data.dom[:el_pos]].ϕₗ
        ϕₗ_sep_pos = donnan_equlibrium_potential(u[ϕₗ_pos], c_pos, z, cf, zf, temp_pos)

        idx_temp_i = data.idx[data.dom[:sep]].temp
        temp_i = u[idx_temp_i]

        # diffusion-driven fluxes through separator
        num_species = count_all_species(data.electrolyte.species)
        c_sep_l_all = @MVector zeros(num_species)
        donnan_equilibrium_concentrations!(c_sep_l_all, u[ϕₗ_neg], ϕₗ_sep_neg, c_neg, z, temp_i)
        c_sep_r_all = @MVector zeros(num_species)
        donnan_equilibrium_concentrations!(c_sep_r_all, u[ϕₗ_pos], ϕₗ_sep_pos, c_pos, z, temp_i)

        σₑ_sep_eval = zero(eltype(u))

        species_sep = data.sep.species
        species_el = data.electrolyte.species
        for idx_species in 1:num_species
            diff = effective_diffusion(species_sep[row=idx_species], temp_i)
            z = species_el[row=idx_species, col=Key("charge")]
            c_upwind = v >= 0 ? c_sep_l_all[idx_species] : c_sep_r_all[idx_species]

            σₑ_sep_eval += z^2 * c_upwind * diff
        end

        f[idx_temp_i] -= lx_sep * (σₑ_sep_eval / temp_i) * ((ϕₗ_sep_pos - ϕₗ_sep_neg) / lx_sep)^2 # TODO: CHECK THIS!!
    end
end

function system_bstorage!(f, u, node, data)
    if data.study.non_isothermal
        PE0 = data.scaling_params.PE0
        temp_i = data.idx[data.dom[:sep]].temp
        f[temp_i] = u[temp_i] * PE0
    end
end

function num_variables(var)
    nvar = 0
    for domain in var
        nvar = max(nvar, domain.p, domain.ϕₛ, domain.ϕₗ, maximum(domain.c), domain.temp)
    end
    return nvar
end
