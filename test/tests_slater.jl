using Test
using QuadGK
using Cubature

include("../src/main.jl")

# define 1d slaterian -> transformed into multi-d slaterian
G3 = slater(3.0, 4.0)
G2 = slater([3.0], [4.0;;])
@test(G2.m == G3.m)
@test(G2.Σ == G3.Σ)

# define several slaterian distribution
G0 = slater([3.0; 0.0], [4.0 0.1; 0.1 3.0])
G1 = slater([3.0; 3.0], [4.0 0.0; 0.0 2.0])
@test(isvalid(G0))
@test(isvalid(G1))

G0b = slater([3.0; 0.0], [4.0 0.1; 0.2 3.0]) #not symmetric
@test(!isvalid(G0b))

# testing function density
G0 = slater([3.0; 0.0], [4.0 0.1; 0.1 3.0])
G1 = slater([3.0; 3.0], [4.0 0.0; 0.0 2.0])
density(G1, [3.0; 4.0])
density.(Ref(G0), [[3.0; 4.0], [3.0; 4.0]])


# 1d test on normalisation
f(x) = density(G3, [x])
R, e = quadgk(f, -20.0, 20.0, rtol=1e-8)
@test(abs.(R - 1.0) < 1e-5)


# 2d test on normalisation
f(x) = density(G0, [x[1]; x[2]])
R, e = hcubature(f, [-20.0, -20.0], [20.0, 20.0], abstol=1e-8)
@test(abs.(R - 1.0) < 1e-6)



# more tests... W2, w_kl, bar
# + mixtures
t = rand()
gs = [G0, G1]
λs = [t, 1 - t]

gbar, c = barmulti(λs, gs::Vector{loc_scat{Slater}}; Niter=25)
gbar2 = bar(1 - t, G0, G1)
@test(gbar.m == gbar2.m)
@test(abs.(gbar.Σ[1] - gbar2.Σ[1]) < 1e-5)
