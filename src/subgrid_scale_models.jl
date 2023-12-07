function effective_diffusivity_factor(vel_long, vel_trans, region, data)
    PE0 = data.scaling_params.PE0
    if region == data.dom[:el_neg]
        total_dispersion_factor = data.el_neg.Deff_const
        total_dispersion_factor += data.el_neg.Deff_linear * abs(PE0 * vel_long)
        total_dispersion_factor += data.el_neg.Deff_quadratic * (PE0 * vel_long)^2
        return total_dispersion_factor
    elseif region == data.dom[:el_pos]
        total_dispersion_factor = data.el_pos.Deff_const
        total_dispersion_factor += data.el_pos.Deff_linear * abs(PE0 * vel_long)
        total_dispersion_factor += data.el_pos.Deff_quadratic * (PE0 * vel_long)^2
        return total_dispersion_factor
    end
    return 0.0
end

function effective_diffusion(species_params, temp)
    diff0 = species_params("reference_diffusivity")
    temp0 = species_params("reference_temperature")
    α = species_params("temperature_coefficient")
    return diff0 * (1.0 + α * (temp - temp0))
end


function _donnan_equilibrium_function(ϕ_sep, p)
    ϕ, c, z, cf, zf, temp = p
    f = zero(eltype(c))
    for idx in eachindex(c)
        f += z[idx] * c[idx] * exp(z[idx]/temp * (ϕ - ϕ_sep))
    end
    f += cf * zf
end

# compute the species concentrations inside an ion-exchange membrane
function donnan_equilibrium_concentrations!(c_sep, ϕ, ϕ_sep, c, z, temp)
    for idx in eachindex(c_sep)
        c_sep[idx] = c[idx] * exp(z[idx]/temp * (ϕ - ϕ_sep))
    end
end

# compute the Donnan equilibrium potential inside an ion-exchange membrane
function donnan_equlibrium_potential(ϕ, c, z, cf, zf, temp)
    p = (ϕ=ϕ, c=c, z=z, cf=cf, zf=zf, temp=temp)
    ϕ_init = ϕ
    prob = NonlinearProblem(_donnan_equilibrium_function, ϕ_init, p)
    sol = solve(prob, SimpleTrustRegion(), reltol = 1e-8)

    return sol.u
end
