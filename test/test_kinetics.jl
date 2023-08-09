module TestKinetics

using RfbScFVM
using Test

function test()

    α = 0.5
    ki = 1.0
    ν_el = -1
    idx_ox = 1
    idx_red = 2
    temp_ref = 1.0
    r = ReactionParameters{Float64, 2}(α=α, ki=ki, ν_el=ν_el, Δϕ₀=0.5,
                                    idx_ox=idx_ox, idx_red=idx_red, temp_ref=temp_ref)
    ηb = 1.0
    c_ox_b = 1.5
    c_red_b = 1.0
    temp = 1.0
    KI0 = 1.0
    sh = 1.0
    aᵥ = 1.0

    function volumetric_current_density_sh_large(r, ηb, c_ox_b, c_red_b, aᵥ, KI0)
        n = -r.ν_el
        return aᵥ * n * KI0 * c_ox_b * r.ki * (c_red_b/c_ox_b)^r.α *
                (exp(n*(1-r.α)*ηb) - exp(-n*r.α*ηb))
    end

    sh = 1e+12
    iᵥ = volumetric_current_density(r, ηb, c_ox_b, c_red_b, temp, sh, aᵥ, KI0)
    iᵥ_sh_large = volumetric_current_density_sh_large(r, ηb, c_ox_b, c_red_b, aᵥ, KI0)
    @test iᵥ ≈ iᵥ_sh_large

    function volumetric_current_density_sh_small(r, ηb, c_ox_b, c_red_b, aᵥ, KI0, sh)
        n = -r.ν_el
        return aᵥ * n * KI0 * c_ox_b * sh * (exp(n*ηb) - 1) / (1 + c_ox_b/c_red_b * exp(n*ηb))
    end

    sh = 1e-8
    iᵥ = volumetric_current_density(r, ηb, c_ox_b, c_red_b, temp, sh, aᵥ, KI0)
    iᵥ_sh_small = volumetric_current_density_sh_small(r, ηb, c_ox_b, c_red_b, aᵥ, KI0, sh)
    @test iᵥ ≈ iᵥ_sh_small
end

end
