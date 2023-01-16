using Plots
pyplot()

using PyCall
ot = pyimport("ot")

include("../../src/main.jl")
# ---------------------------------------
# 
# Example with 2d Wigner semicircle mixtures
# 
# ---------------------------------------
nx = 50 #grid parameter
xmin = 0.0
xmax = 1.0
X = range(xmin, xmax, length=nx) |> collect;
XX = [[X[i]; X[j]] for i = 1:nx, j = 1:nx]

m0 = [0.3; 0.6]
Σ0 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(semicircle(m0, Σ0)))
g0 = semicircle(m0, Σ0)

m1 = [0.7; 0.7]
Σ1 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(semicircle(m1, Σ1)))
g1 = semicircle(m1, Σ1)

m2 = [0.5; 0.6]
Σ2 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(semicircle(m2, Σ2)))
g2 = semicircle(m2, Σ2)

m3 = [0.4; 0.25]
Σ3 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(semicircle(m3, Σ3)))
g3 = semicircle(m3, Σ3)

M0 = mixture(2, [0.3; 0.7], [g0; g1])
M1 = mixture(2, [0.4; 0.6], [g2; g3])

# MW2 barycenter
baryvec = []
barydens = []
T = range(0.0, 1.0, length=5)
for t in T
   bt = bar(t, M0, M1)
   push!(baryvec, bt)
   push!(barydens, density.(Ref(bt), XX))
end

for (i, t) in enumerate(T)
   P = contour(X, X, barydens[i], c=:viridis, aspect_ratio=:equal, clim=(0.0, 8.0), fill=:false, yformatter=_->"", xformatter=_->"", colorbar =:false)
   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section4_semicircle/2d/bary_mw2_2d_contour$i.pdf")
end

# W2 barycenter computed with Sinkhorn algorithm
XY = [[X[i]; X[j]] for i = 1:nx for j = 1:nx]
g0dens = 1e-16 .+ density.(Ref(M0), XY)
g1dens = 1e-16 .+ density.(Ref(M1), XY)
A = zeros(length(g0dens), 2)
A[:, 1] = g0dens
A[:, 2] = g1dens

M = zeros(nx^2,nx^2)
for i in 1:nx^2 
   for j in 1:nx^2
      M[i,j] = norm(XY[i]-XY[j])^2
   end
end
M = M ./ maximum(M)
reg = 1e-4 #regularization parameter

baryw2 = [g0dens]
for t in T[2:end-1]
   @show t
   bt = ot.barycenter(A, M, reg, [1 - t; t], method="sinkhorn_log", numItermax=10000)
   push!(baryw2, bt)
end
push!(baryw2,g1dens)

for (i, t) in enumerate(T)
   P = contour(X, X, (reshape(baryw2[i], nx, nx))', c=:viridis, aspect_ratio=:equal, clim=(0.0, 8.0), fill=:false, yformatter=_->"", xformatter=_->"", colorbar =:false)
   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section4_semicircle/2d/bary_w2_2d_contour$i.pdf")
end