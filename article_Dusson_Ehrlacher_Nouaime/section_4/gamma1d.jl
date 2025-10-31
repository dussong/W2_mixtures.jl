# The goal here is to compare the defined barycenter based on gamma1d mixtures with the standard W2 barycenter.
using Plots
pyplot()

using Dierckx

# using PyCall
# ot = pyimport("ot")

include("../../src/main.jl")

nx = 200
xmin = -4.0
xmax = 2.0
X = range(xmin, xmax, length=nx) |> collect
XX = [[X[i]] for i = 1:nx]


nf = 10000
Xfin = range(xmin, xmax, length=nf) |> collect
XXfin = [[Xfin[i]] for i = 1:nf]
X01 = range(0.,1.,length = nf) |> collect

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

ρ0 = density.(Ref(M0), XXfin)
ρ1 = density.(Ref(M1), XXfin)

# cdf for ρ
sρ0 = Spline1D(Xfin, ρ0, bc="nearest", k=3) #spline representation of density
Tf0 = [integrate(sρ0, Xfin[1], y) for y in Xfin] #T: cumulative distr. function
Tf0 = [max(Tf,0.) for Tf in Tf0]
Tf0 = Tf0 ./ maximum(Tf0) #transport map
indf0 = [findall(x -> (x > 1e-6), Tf0[2:end] - Tf0[1:end-1]); length(Tf0)] #remove singularity
icdf0 = Spline1D(Tf0[indf0], Xfin[indf0], bc="nearest", k=3) #spline representation of F (inverse of cdf)
evalicdf0 = evaluate(icdf0, X01)
evalicdf0[1] = X[1]
evalicdf0[end] = X[end]

# cdft = Spline1D(evalicdf0, X01, bc="nearest", k=3)
# bt = Dierckx.derivative(cdft, Xfin)


# plot(evaluate(icdf0, collect(0:.01:1)))
# cdft = Spline1D(X, icdfbt, bc="nearest", k=1)

sρ1 = Spline1D(Xfin, ρ1, bc="nearest", k=3) #spline representation of density
Tf1 = [integrate(sρ1, Xfin[1], y) for y in Xfin] #T: cumulative distr. function
Tf1 = [max(Tf,0.) for Tf in Tf1]
Tf1 = Tf1 ./ maximum(Tf1) #transport map
indf1 = [findall(x -> (x > 1e-6), Tf1[2:end] - Tf1[1:end-1]); length(Tf1)] #remove singularity
icdf1 = Spline1D(Tf1[indf1], Xfin[indf1],bc="nearest", k=3) #spline representation of F (inverse of cdf)
evalicdf1 = evaluate(icdf1, X01)
evalicdf1[1] = X[1]
evalicdf1[end] = X[end]


# true barycenter (W2, computed with icdf0)
baryw2 = [ρ0]
for t in T[2:end-1]
   @show t
   icdfbt = (1-t)*evalicdf0 + t*evalicdf1
   cdft = Spline1D(icdfbt,X01, bc="nearest", k=1)
   bt = Dierckx.derivative(cdft, Xfin)
   push!(baryw2, bt)
end
push!(baryw2, ρ1)


# # W2 barycenter
# g0dens = 1e-16 .+ density.(Ref(M0), XX)
# g1dens = 1e-16 .+ density.(Ref(M1), XX)
# A = zeros(length(g0dens), 2)
# A[:, 1] = g0dens
# A[:, 2] = g1dens

# x = reshape(X, length(X), 1)
# M = ot.dist(x, x)
# M = M ./ maximum(M)
# reg = 1e-4

# # true barycenter (W2, computed with a Sinkhorn algorithm)
# baryw2 = [g0dens]
# for t in T[2:end-1]
#    @show t
#    @time bt = ot.barycenter(A, M, reg, [1 - t; t], method="sinkhorn_log", numItermax=10000)
#    push!(baryw2, bt)
# end
# push!(baryw2, g1dens)


# Comparison between the two barycenters
for (i, t) in enumerate(T)
   P = plot(Xfin, baryw2[i], label="W2", ylim=(0.0, 2.0), color=1, linewidth=4, yformatter=_ -> "", xformatter=_ -> "") #, 
   P = plot!(X, barydens[i], label="", color=2, tickfontsize=18, linewidth=4) #, label="MW2"

   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section4_gamma1d/bary_w2_mw2_1d$i.pdf")
end