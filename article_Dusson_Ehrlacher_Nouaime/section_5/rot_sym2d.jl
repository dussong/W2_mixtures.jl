using Plots
pyplot()

using PyCall
ot = pyimport("ot")

include("../../src/main.jl")

# ---------------------------------------
# 
# Second example 2d
# 
# ---------------------------------------
nx = 50
xmin = -3.0
xmax = 3.0
X = range(xmin, xmax, length=nx) |> collect;
XX = [[X[i]; X[j]] for i = 1:nx, j = 1:nx]

SN = rot_group()
ga = rot_sym()

m0 = [0.0; 0.5]
Σ0 = [0.001 0.0; 0.0 0.001]
@assert(isvalid(slater(m0, Σ0)))
a0 = slater(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = [0.0; 1.5]
Σ1 = [0.02 0.0; 0.0 0.02]
@assert(isvalid(slater(m1, Σ1)))
a1 = slater(m1, Σ1)
A1 = gsa(a1, SN, ga)

m2 = [0.0; 0.5]
Σ2 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(slater(m2, Σ2)))
a2 = slater(m2, Σ2)
A2 = gsa(a2, SN, ga)

m3 = [0.0; 2.5]
Σ3 = [0.01 0.00; 0.00 0.01]
@assert(isvalid(slater(m3, Σ3)))
a3 = slater(m3, Σ3)
A3 = gsa(a3, SN, ga)

M0 = mixture(2, [0.5; 0.5], [A0; A1])
M1 = mixture(2, [0.2; 0.8], [A2; A3])

baryvec = []
barydens = []
T = range(0.0, 1.0, length=5)
for t in T
   bt = bar(t, M0, M1)
   push!(baryvec, bt)
   push!(barydens, density.(Ref(bt), XX))
end

for (i, t) in enumerate(T)
   P = contour(X, X, barydens[i], c=:viridis, aspect_ratio=:equal, clim=(0.0, 0.22), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false)
   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_rotsym/bary_mw2_2d_contour$i.pdf")
end
# 


# for (i, t) in enumerate(T)
#    P = surface(X, X, barydens[i], c=:viridis, grid=:false, zlim=(0.0, 3.0), camera=(10, 50), colorbar=:false, yformatter=_ -> "", xformatter=_ -> "", zformatter=_ -> "")
#    display(P)
#    sleep(0.5)
#    savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_rotsym/bary_mw2_2d_surface$i.pdf")
# end

XY = [[X[i]; X[j]] for i = 1:nx for j = 1:nx]
a0dens = 1e-16 .+ density.(Ref(M0), XY)
a1dens = 1e-16 .+ density.(Ref(M1), XY)
A = zeros(length(a0dens), 2)
A[:, 1] = a0dens
A[:, 2] = a1dens

M = zeros(nx^2, nx^2)
for i in 1:nx^2
   for j in 1:nx^2
      M[i, j] = norm(XY[i] - XY[j])^2
   end
end
M = M ./ maximum(M)
reg = 1e-4

# true barycenter (W2, computed with a Sinkhorn algorithm)
baryw2 = [a0dens]
for t in T[2:end-1]
   @show t
   @time bt = ot.barycenter(A, M, reg, [1 - t; t], method="sinkhorn_log", numItermax=10000)
   push!(baryw2, bt)
end
push!(baryw2, a1dens)

for (i, t) in enumerate(T)
   P = contour(X, X, (reshape(baryw2[i], nx, nx))', c=:viridis, aspect_ratio=:equal, clim=(0.0, 0.22), colorbar=:false, yformatter=_ -> "", xformatter=_ -> "", zformatter=_ -> "")
   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_rotsym/bary_w2_2d_contour$i.pdf")
end

# for (i, t) in enumerate(T)
#    P = surface(X, X, (reshape(baryw2[i], nx, nx))', c=:viridis, grid=:false, zlim=(0.0, 3.0), camera=(10, 50), colorbar=:false, yformatter=_ -> "", xformatter=_ -> "", zformatter=_ -> "")
#    display(P)
#    sleep(0.5)
#    savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_rotsym/bary_w2_2d_surface$i.pdf")
# end