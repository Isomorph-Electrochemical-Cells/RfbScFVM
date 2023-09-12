module TestDonnanPotential

using RfbScFVM
using Test
#using StaticArrays
#using NonlinearSolve

function test()
    ϕ = 0.5
    c = [0.001, 0.001]
    z = [1.0, -1.0]
    cf = 1.25
    temp = 2.5
    c_sep = similar(c)

    # Test Donnan exclusion principle for an anion-exchange membrane
    zf = 1.0
    ϕ_sep = donnan_equlibrium_potential(ϕ, c, z, cf, zf, temp)
    donnan_equilibrium_concentrations!(c_sep, ϕ, ϕ_sep, c, z, temp)

    @test c_sep[1] ≈ zero(c_sep[1]) atol=1e-4
    @test c_sep[2] ≈ cf rtol=1e-4

    # Test Donnan exclusion principle for a cation-exchange membrane
    zf = -1.0
    ϕ_sep = donnan_equlibrium_potential(ϕ, c, z, cf, zf, temp)
    donnan_equilibrium_concentrations!(c_sep, ϕ, ϕ_sep, c, z, temp)

    @test c_sep[1] ≈ cf rtol=1e-4
    @test c_sep[2] ≈ zero(c_sep[2]) atol=1e-4

    # Special case when only counter-ion species are present
    zf = 1.0
    c = [0.0, 1.0]
    ϕ_sep = donnan_equlibrium_potential(ϕ, c, z, cf, zf, temp)

    # check if the numerical solution well-approximates the analytical solution of ϕ
    ϕ_sep_exact = ϕ - log(-cf * zf / (c[2] * z[2])) * temp / z[2]
    @test ϕ_sep ≈ ϕ_sep_exact rtol=1e-6

    # check that only counter-ion species are present in the membrane
    donnan_equilibrium_concentrations!(c_sep, ϕ, ϕ_sep, c, z, temp)
    @test c_sep[1] ≈ zero(c_sep[1]) atol=1e-4
    @test c_sep[2] ≈ cf rtol=1e-4

    # check Donnan potential in case of arbitrary number of species
    n = 4 # number of species
    c = vcat(1e-3 * ones(n-1), 1.0)
    z = vcat(1:(n-1), -1.0)
    c_sep = similar(c)

    cf = 1.5
    temp = 2.5
    zf = 4.0

    ϕ_sep = donnan_equlibrium_potential(ϕ, c, z, cf, zf, temp)
    donnan_equilibrium_concentrations!(c_sep, ϕ, ϕ_sep, c, z, temp)

    # check Donnan exclusion principle for an anion-exchange membrane
    @test all(c_sep[1:(n-1)] .< 1e-3)
    # check electroneutrality condition within membrane
    @test sum(z.*c_sep) + zf*cf ≈ 0.0 atol=1e-8
    # check equality of electrochemical potentials
    μ̃ = temp * log.(c) .+ z .* ϕ
    μ̃_sep = temp * log.(c_sep) .+ z .* ϕ_sep
    @test μ̃ ≈ μ̃_sep rtol=1e-6
end

end
