using Test
using Aqua
using RfbScFVM

@testset "RfbScFVM" begin
    @testset "Kinetics" begin
        include("test_kinetics.jl")
        TestKinetics.test()
    end
    @testset "Aqua" begin
        Aqua.test_project_toml_formatting(RfbScFVM)
        Aqua.test_deps_compat(RfbScFVM)
        Aqua.test_stale_deps(RfbScFVM; ignore=[:Pluto, :PlutoVista])
        Aqua.test_project_extras(RfbScFVM)
        Aqua.test_piracy(RfbScFVM)
        Aqua.test_undefined_exports(RfbScFVM)
        Aqua.test_unbound_args(RfbScFVM)
        Aqua.test_ambiguities(RfbScFVM)
    end
end
