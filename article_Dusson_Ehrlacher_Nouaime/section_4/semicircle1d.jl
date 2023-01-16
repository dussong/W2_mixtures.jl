#We compare the MW2 barycenter based on 1d Wigner semicircle mixtures with the standard W2 barycenter.
using Plots
pyplot()

using PyCall
ot = pyimport("ot")

include("../../src/main.jl")

nx = 200
xmin = 0.0
xmax = 0.8
X = range(xmin, xmax, length=nx) |> collect
XX = [[X[i]] for i = 1:nx]


# Barycenter between two sym semicircle functions
m0 = 0.2
Σ0 = 1e-2
@assert(isvalid(semicircle(m0, Σ0)))
g0 = semicircle(m0, Σ0)

m1 = 0.4
Σ1 = 1e-3
@assert(isvalid(semicircle(m1, Σ1)))
g1 = semicircle(m1, Σ1)

m2 = 0.4
Σ2 = 0.002
@assert(isvalid(semicircle(m2, Σ2)))
g2 = semicircle(m2, Σ2)

m3 = 0.6
Σ3 = 0.003
@assert(isvalid(semicircle(m3, Σ3)))
g3 = semicircle(m3, Σ3)

M0 = mixture(2, [0.3; 0.7], [g0; g1])
M1 = mixture(2, [0.6; 0.4], [g2; g3])

# MW2 barycenter
baryvec = []
barydens = []
T = range(0.0, 1.0, length=5)
for t in T
   bt = bar(t, M0, M1)
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

#computed with a Sinkhorn algorithm
baryw2 = [g0dens]
for t in T[2:end-1]
   @show t
   bt = ot.barycenter(A, M, reg, [1 - t; t], method="sinkhorn_log", numItermax=10000)
   push!(baryw2, bt)
end
push!(baryw2, g1dens)

# comparison between the two barycenters
for (i, t) in enumerate(T)
   P = plot(X, baryw2[i], label="", ylim=(0.0, 8.0), color=1, linewidth=4, yformatter=_ -> "", xformatter=_ -> "") #, label="W2"
   P = plot!(X, barydens[i], label="", color=2, tickfontsize=18, linewidth=4) #, label="MW2"

   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section4_semicircle/1d/bary_w2_mw2_1d$i.pdf")
end


