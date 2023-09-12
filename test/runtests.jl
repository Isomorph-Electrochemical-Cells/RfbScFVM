using Test
using Aqua
using RfbScFVM

@testset "Kinetics" begin
    include("test_kinetics.jl")
    TestKinetics.test()
end

@testset "DonnanPotential" begin
    include("test_donnan_potential.jl")
    TestDonnanPotential.test()
end

@testset "TestModel1Dvs2D" begin
    include("test_model_1d_vs_2d.jl")
    TestModel1Dvs2D.test()
end

@testset "TestHomogeneousReactions" begin
    include("test_homogeneous_reactions.jl")
    TestHomogeneousReactions.test()
end

@testset "Aqua" begin
    Aqua.test_project_toml_formatting(RfbScFVM)
    Aqua.test_deps_compat(RfbScFVM)
    Aqua.test_stale_deps(RfbScFVM; ignore=[:Pluto, :PlutoVista, :PlutoUI])
    Aqua.test_project_extras(RfbScFVM)
    Aqua.test_piracy(RfbScFVM)
    Aqua.test_undefined_exports(RfbScFVM)
    Aqua.test_unbound_args(RfbScFVM)
    Aqua.test_ambiguities(RfbScFVM)
end
