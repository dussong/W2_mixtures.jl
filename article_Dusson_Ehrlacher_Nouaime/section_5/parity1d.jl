# Comparison between W2 and MW2 barycenters of parity symmetrized densities

using Plots
pyplot()

using PyCall
ot = pyimport("ot")

include("../../src/main.jl")

nx = 200
xmin = -1.0
xmax = 1.0
X = range(xmin, xmax, length=nx) |> collect
XX = [[X[i]] for i = 1:nx]

SN = perm_group(2) #define group
ga = parity() #define group action

#Definition of mixtures
m0 = -0.2
Σ0 = 0.0009
@assert(isvalid(slater(m0, Σ0)))
a0 = slater(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = 0.4
Σ1 = 0.0016
@assert(isvalid(slater(m1, Σ1)))
a1 = slater(m1, Σ1)
A1 = gsa(a1, SN, ga)

m2 = 0.3
Σ2 = 0.0036
@assert(isvalid(slater(m2, Σ2)))
a2 = slater(m2, Σ2)
A2 = gsa(a2, SN, ga)

m3 = 0.5
Σ3 = 0.0050
@assert(isvalid(slater(m3, Σ3)))
a3 = slater(m3, Σ3)
A3 = gsa(a3, SN, ga)

M0 = mixture(2, [0.3; 0.7], [A0; A1])
M1 = mixture(2, [0.6; 0.4], [A2; A3])

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
a0dens = 1e-16 .+ density.(Ref(M0), XX)
a1dens = 1e-16 .+ density.(Ref(M1), XX)
A = zeros(length(a0dens), 2)
A[:, 1] = a0dens
A[:, 2] = a1dens

x = reshape(X, length(X), 1)
M = ot.dist(x, x)
M = M ./ maximum(M)
reg = 1e-4 #regularization parameter

baryw2 = [a0dens]
for t in T[2:end-1]
   @show t
   bt = ot.barycenter(A, M, reg, [1 - t; t], method="sinkhorn_log", numItermax=10000)
   push!(baryw2, bt)
end
push!(baryw2, a1dens)

# Comparison between the two barycenters
for (i, t) in enumerate(T)
   P = plot(X, baryw2[i], label="", ylim=(0.0, 6.0), color=1, linewidth=4, yformatter=_ -> "", xformatter=_ -> "") #, label="W2"
   P = plot!(X, barydens[i], label="", color=2, tickfontsize=18, linewidth=4) #, label="MW2"

   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_parity/bary_w2_mw2_1d$i.pdf")
end

# Comparison between the optimal transport plans
nx = 100
xmin = -1.0
xmax = 1.0
X = range(xmin, xmax, length=nx) |> collect
XX = [[X[i]] for i = 1:nx]

SN = perm_group(2)
ga = parity()

# Barycenter between two slater functions
m0 = -0.2
Σ0 = .05
@assert(isvalid(slater(m0, Σ0)))
a0 = slater(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = 0.4
Σ1 = .05
@assert(isvalid(slater(m1, Σ1)))
a1 = slater(m1, Σ1)
A1 = gsa(a1, SN, ga)

a0dens_sym = 1e-16 .+ density.(Ref(A0), XX)
a1dens_sym = 1e-16 .+ density.(Ref(A1), XX)

plot(X, a0dens_sym, label="", color=1, linewidth=4, yformatter=_ -> "", xformatter=_ -> "")
plot!(X, a1dens_sym, label="", color=2, linewidth=4, yformatter=_ -> "", xformatter=_ -> "")
savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_parity/otplan_densities.pdf")

x = reshape(X, length(X), 1)
M = ot.dist(x, x)
M = M ./ maximum(M)

d, g = _W2(A0, A1) #extract optimal solution (group element)
a1_adapted = group_action(ga, a1, g)

a0dens = 1e-16 .+ density.(Ref(a0), XX)
a1dens = 1e-16 .+ density.(Ref(a1_adapted), XX)

G0_W2 = ot.emd(a0dens_sym ./ sum(a0dens_sym), a1dens_sym ./ sum(a1dens_sym), M)
heatmap(X, X, G0_W2, c=:viridis, grid=:false, colorbar=:false, zlim=(0.0, 0.005), aspect_ratio=:equal, tickfontsize=14, yformatter=_ -> "", xformatter=_ -> "")
savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_parity/otplan_w2.pdf")

G0_MW2 = ot.emd(a0dens ./ sum(a0dens), a1dens ./ sum(a1dens), M)
heatmap(X, X, G0_MW2, c=:viridis, grid=:false, colorbar=:false, zlim=(0.0, 0.005), aspect_ratio=:equal, tickfontsize=14)

G0_MW2sym = G0_MW2[end:-1:1, end:-1:1]

heatmap(X, X, 1 / 2 * (G0_MW2 + G0_MW2sym), c=:viridis, grid=:false, colorbar=:false, zlim=(0.0, 0.005), aspect_ratio=:equal, tickfontsize=14, yformatter=_ -> "", xformatter=_ -> "")
savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_parity/otplan_w2sym.pdf")