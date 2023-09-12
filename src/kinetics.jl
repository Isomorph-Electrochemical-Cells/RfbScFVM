
function regularize_non_negative_variable(var::T) where T
    threshold::T = 100*eps(T)
    min_value::T = 10*eps(T)
    if var < threshold
        return var < zero(var) ? min_value :
                                 min_value + (var/threshold)*(threshold-min_value)
    end
    return var
end

function volumetric_current_density(r::ReactionParameters{T, NSPEC},
                                    ηb, c_ox_b, c_red_b, v, temp,
                                    region, data) where {T, NSPEC}
    iv = zero(c_ox_b)
    if region == data.dom[:el_neg] || region == data.dom[:el_pos]
        if region == data.dom[:el_neg]
            aᵥ = data.el_neg.aᵥ
            sh_factor = data.el_neg.sh_factor
            sh_exponent = data.el_neg.sh_exponent
        else
            aᵥ = data.el_pos.aᵥ
            sh_factor = data.el_pos.sh_factor
            sh_exponent = data.el_pos.sh_exponent
        end

        @inline c_ox_b = regularize_non_negative_variable(c_ox_b)
        @inline c_red_b = regularize_non_negative_variable(c_red_b)

        cr = c_ox_b/c_red_b
        n::T = -r.ν_el

        KI0 = data.scaling_params.KI0
        sh = sh_factor * v^sh_exponent

        Δϕ₀ = r.Δϕ₀ + r.∂Δϕ₀_∂temp * (temp - r.temp_ref)
        ki = r.ki * exp(n * Δϕ₀ * (1/r.temp_ref - 1/temp))

        # volumetric current density
        iv = aᵥ*n*c_ox_b*KI0*ki*(exp(n*ηb/temp)-1) / ((1+cr*exp(n*ηb/temp))*ki/sh +
                                                        (cr^r.α)*exp(n*r.α*ηb/temp))
    end
    return iv
end


function volumetric_current_density(r::ReactionParameters{T,NSPEC},
    ηb, c_ox_b, c_red_b, temp, sh, aᵥ::T, KI0::T) where {T,NSPEC}
    iv = zero(c_ox_b)
    #if region == data.dom[:el_neg] || region == data.dom[:el_pos]
        # if region == data.dom[:el_neg]
        #     aᵥ = data.el_neg.aᵥ
        #     sh_factor = data.el_neg.sh_factor
        #     sh_exponent = data.el_neg.sh_exponent
        # else
        #     aᵥ = data.el_pos.aᵥ
        #     sh_factor = data.el_pos.sh_factor
        #     sh_exponent = data.el_pos.sh_exponent
        # end

        @inline c_ox_b = regularize_non_negative_variable(c_ox_b)
        @inline c_red_b = regularize_non_negative_variable(c_red_b)

        cr = c_ox_b / c_red_b
        n = -r.ν_el
        #KI0 = data.scaling_params.KI0
        #sh = sh_factor * v^sh_exponent

        Δϕ₀ = r.Δϕ₀ + r.∂Δϕ₀_∂temp * (temp - r.temp_ref)
        ki = r.ki * exp(n * Δϕ₀ * (1 / r.temp_ref - 1 / temp))

        # volumetric current density
        iv = aᵥ * n * c_ox_b * KI0 * ki * (exp(n * ηb / temp) - 1) / ((1 + cr * exp(n * ηb / temp)) * ki / sh +
                                                                      (cr^r.α) * exp(n * r.α * ηb / temp))
   #end
    return iv
end


# function ϕ_eq(r::ReactionParameters{T, NSPEC}, u, idx_c2u, temp) where {T, NSPEC}
#     reaction_quotient = one(T)
#     for idx in eachindex(r.ν_coeff) # TODO: Handle case, where counter species is involved in reaction!
#         if r.ν_coeff[idx] == zero(Int64)
#             continue
#         end
#         reaction_quotient *= u[idx_c2u[idx]]^r.ν_coeff[idx]
#     end
#     r.Δϕ₀ + temp * (one(T)/r.ν_el) * log(reaction_quotient)
# end

function ϕ_eq(r::ReactionParameters{T, NSPEC}, c_fun, temp) where {T, NSPEC}
    reaction_quotient = one(T)
    for idx in 1:NSPEC
        if r.ν_coeff[idx] == zero(Int64)
            continue
        end
        reaction_quotient *= c_fun(idx)^r.ν_coeff[idx]
    end
    if (reaction_quotient < - SQRT_EPS)
        @warn "negative reaction quotient" Float64(reaction_quotient)
    end

    r.Δϕ₀ + temp * (one(T)/r.ν_el) * log(abs(reaction_quotient) + eps(reaction_quotient))
end

function ϕ_eq(vec_r::VECTOR_REACTION_PARAMETERS, c_fun, temp) where {T, NSPEC, VECTOR_REACTION_PARAMETERS <: AbstractVector{ReactionParameters{T, NSPEC}}}
    ϕ_eq_result = zero(T)
    num_electrons = zero(Int64)
    for r in vec_r
        reaction_quotient = one(T)
        for idx in 1:NSPEC
            if r.ν_coeff[idx] == zero(Int64)
                continue
            end
            reaction_quotient *= c_fun(idx)^r.ν_coeff[idx]
        end
        num_electrons += r.ν_el

        if (reaction_quotient < - SQRT_EPS)
            @warn "negative reaction quotient" Float64(reaction_quotient)
        end

        ϕ_eq_result += r.ν_el * r.Δϕ₀ + temp * log(abs(reaction_quotient) + eps(reaction_quotient))
    end
    ϕ_eq_result /= num_electrons
end
