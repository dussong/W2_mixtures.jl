using Test

@testset "MW2.jl" begin
    # --------------------------------------------
    @testset "gaussians" begin include("tests_gaussians.jl") end
    @testset "mixtures" begin include("tests_mixtures.jl") end
    @testset "gamma1d" begin include("tests_gamma1d.jl") end
    @testset "semicircle" begin include("tests_semicircle.jl") end
    @testset "slater" begin include("tests_slater.jl") end
    @testset "symm_mixtures" begin include("tests_group_sym_atoms.jl") end
end