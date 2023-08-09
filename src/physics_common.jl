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


function concentration_counter_species(u, region, data, node)
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
            if c <= -1e-8
                @warn c idx_species indices[idx_species] region
            end
            c_counter -= z * c
        end
        c_counter /= data.electrolyte.counter[col=Key("charge")]
        if c_counter <= -1e-8
            @warn c_counter region
        end
        return c_counter #max(c_counter, zero(c_counter))
    end
    return zero(eltype(u))
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
        vₑ = data.sep.kϕ / data.sep.μ * data.sep.cf * data.sep.zf * (u[ϕₗ_l] - u[ϕₗ_l]) / h
        vₑ /= (data.scaling_params.PE0)
        v = vₕ + vₑ

        f[p_l] = v # left volumetric flux
        f[p_r] = -v # right volumetric flux

        # evaluate Donnan potentials
        z_c = data.electrolyte.counter[col=Key("charge")]
        @inline c_c_el_l = concentration_counter_species(u, data.dom[:el_neg], data, 0) + eps(eltype(u))
        @assert c_c_el_l > 0
        @inline c_c_el_r = concentration_counter_species(u, data.dom[:el_pos], data, 0) + eps(eltype(u))

        @assert c_c_el_r > 0
        c_c_sep = -data.sep.cf * data.sep.zf / z_c
        @assert c_c_sep > 0
        ϕₗ_sep_l = u[ϕₗ_l] + one(eltype(u))/z_c * log(c_c_el_l / c_c_sep)
        ϕₗ_sep_r = u[ϕₗ_r] + one(eltype(u))/z_c * log(c_c_el_r / c_c_sep)

        σₑ_sep = data.sep.σₑ
        if data.study.non_isothermal
            temp_i = data.idx[data.dom[:sep]].temp
            σₑ_sep += data.sep.∂σₑ∂temp * (u[temp_i] - data.sep.temp_ref)
        end
        current_flux_sep = -σₑ_sep * (ϕₗ_sep_r - ϕₗ_sep_l) / h

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


function regularize_non_negative_quantity(u)
    return u <= zero(u) ? zero(u) : u
end

function system_bstorage!(f, u, node, data)
    if node.region == data.dom[:sep] && data.study.non_isothermal
        PE0 = data.scaling_params.PE0
        temp_i = data.idx[data.dom[:sep]].temp
        f[temp_i] = u[temp_i] * PE0
    end
end
