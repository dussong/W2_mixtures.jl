# The goal here is to compare the defined barycenter based on gamma1d mixtures with the standard W2 barycenter.
using Plots
pyplot()

using PyCall
ot = pyimport("ot")

include("../../src/main.jl")

nx = 200
xmin = -4.0
xmax = 2.0
X = range(xmin, xmax, length=nx) |> collect
XX = [[X[i]] for i = 1:nx]

α = 3.0

# Barycenter between two sym gamma1d functions
m0 = -2.0
Σ0 = 0.1
g0 = gamma1d(m0, Σ0, α)
@assert(isvalid(g0))


m1 = 1.0
Σ1 = 1e-2
g1 = gamma1d(m1, Σ1, α)
@assert(isvalid(g1))


m2 = 1.0
Σ2 = 0.3
g2 = gamma1d(m2, Σ2, α)
@assert(isvalid(g2))

m3 = 0.0
Σ3 = 0.03
g3 = gamma1d(m3, Σ3, α)
@assert(isvalid(g3))


M0 = mixture(2, [0.7; 0.3], [g0; g1])
M1 = mixture(2, [0.4; 0.6], [g2; g3])

# MW2 barycenter
baryvec = []
barydens = []
T = range(0.0, 1.0, length=5)
for t in T
   @time bt = bar(t, M0, M1)
   push!(baryvec, bt)
   push!(barydens, density.(Ref(bt), XX))
end

# W2 barycenter
g0dens = 1e-16 .+ density.(Ref(M0), XX)
g1dens = 1e-16 .+ density.(Ref(M1), XX)
A = zeros(length(g0dens), 2)
A[:, 1] = g0dens
A[:, 2] = g1dens

x = reshape(X, length(X), 1)
M = ot.dist(x, x)
M = M ./ maximum(M)
reg = 1e-4

# true barycenter (W2, computed with a Sinkhorn algorithm)
baryw2 = [g0dens]
for t in T[2:end-1]
   @show t
   @time bt = ot.barycenter(A, M, reg, [1 - t; t], method="sinkhorn_log", numItermax=10000)
   push!(baryw2, bt)
end
push!(baryw2, g1dens)

# Comparison between the two barycenters
for (i, t) in enumerate(T)
   P = plot(X, baryw2[i], label="", ylim=(0.0, 2.0), color=1, linewidth=4, yformatter=_ -> "", xformatter=_ -> "") #, label="W2"
   P = plot!(X, barydens[i], label="", color=2, tickfontsize=18, linewidth=4) #, label="MW2"

   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section4_gamma1d/bary_w2_mw2_1d$i.pdf")
end