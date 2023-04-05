
function regularize_non_negative_variable(var)
    (var < 0) ? eps(typeof(var))*exp(var) : eps(typeof(var))+var
end


function volumetric_current_density(r::ReactionParameters{T},
                                    ηb, c_ox_b, c_red_b, v, temp, region, data) where T
    if region == Int(domain_id_el_neg) || region == Int(domain_id_el_pos)
        if region == Int(domain_id_el_neg)
            aᵥ = data.el_neg.aᵥ
            sh_factor = data.el_neg.sh_factor
            sh_exponent = data.el_neg.sh_exponent
        else
            aᵥ = data.el_pos.aᵥ
            sh_factor = data.el_pos.sh_factor
            sh_exponent = data.el_pos.sh_exponent
        end
        # regularize c_ox and c_red
        c_ox_b = regularize_non_negative_variable(c_ox_b)
        c_red_b = regularize_non_negative_variable(c_red_b)
        cr = c_ox_b/c_red_b
        n = -r.ν_el
        KI0 = data.scaling_params.KI0
        sh = sh_factor * v^sh_exponent

        Δϕ₀ = r.Δϕ₀ + r.∂Δϕ₀_∂temp * (temp - r.temp_ref)
        ki = r.ki * exp(n * Δϕ₀ * (1/r.temp_ref - 1/temp))

        # volumetric current density
        iv = aᵥ*n*c_ox_b*KI0*ki*(exp(ηb/temp)-1) / ((1+cr*exp(ηb/temp))*ki/sh +
                                                      (cr^r.α)*exp(r.α*ηb/temp))

        return iv
    end
    return 0.0
end

function ϕ_eq(r::ReactionParameters{T}, c_ox, c_red, temp) where T
    # regularize c_ox and c_red
    c_ox_prime = regularize_non_negative_variable(c_ox)
    c_red_prime = regularize_non_negative_variable(c_red)
    # evaluate reduction potential
    return r.Δϕ₀ + temp * (one(T)/r.ν_el)*log((c_red_prime)^r.ν_red * (c_ox_prime)^r.ν_ox)
end
