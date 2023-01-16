using Test
using QuadGK
using Cubature

include("../src/main.jl")

# define 1d gamma distribution
α = 3.0
G3 = gamma1d(1.0, 1.0, α)
@test(isvalid(G3))

# 1d test on normalisation
f(x) = density(G3, [x])
R, e = quadgk(f, -20.0, 20.0, rtol=1e-8)
@test(abs.(R - 1.0) < 1e-5)


g(x) = x*density(G3, [x])
R, e = quadgk(g, -20.0, 20.0, rtol=1e-8)
@test(abs.(R - 1.0) < 1e-5)
