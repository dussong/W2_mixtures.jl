using Plots
pyplot()

using PyCall
ot = pyimport("ot")

include("../../src/main.jl")

# ---------------------------------------
# 
# Permutation symmetry arising from antisymmetry in 2D
# 
# ---------------------------------------
nx = 50
xmin = 0.0
xmax = 1.0
X = range(xmin, xmax, length=nx) |> collect;
XX = [[X[i]; X[j]] for i = 1:nx, j = 1:nx]

SN = perm_group(2) #define group
ga = perm_antisym() #define group action

# Define four a-sym mixtures based on gaussian distributions
m0 = [0.3; 0.7]
Σ0 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(gauss(m0, Σ0)))
a0 = gauss(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = [0.69; 0.7]
Σ1 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(gauss(m1, Σ1)))
a1 = gauss(m1, Σ1)
A1 = gsa(a1, SN, ga)

m2 = [0.55; 0.6]
Σ2 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(gauss(m2, Σ2)))
a2 = gauss(m2, Σ2)
A2 = gsa(a2, SN, ga)

m3 = [0.4; 0.2]
Σ3 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(gauss(m3, Σ3)))
a3 = gauss(m3, Σ3)
A3 = gsa(a3, SN, ga)

m4 = [0.5; 0.3]
Σ4 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(gauss(m4, Σ4)))
a4 = gauss(m4, Σ4)
A4 = gsa(a4, SN, ga)

m5 = [0.7; 0.75]
Σ5 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(gauss(m5, Σ5)))
a5 = gauss(m5, Σ5)
A5 = gsa(a5, SN, ga)

m6 = [0.6; 0.2]
Σ6 = [0.01 0.0; 0.0 0.01]
@assert(isvalid(gauss(m6, Σ6)))
a6 = gauss(m6, Σ6)
A6 = gsa(a6, SN, ga)

M0 = mixture(2, [0.5; 0.5], [A0; A1])
M1 = mixture(2, [0.4; 0.6], [A2; A3])
M2 = mixture(3, [0.3; 0.3; 0.4], [A0; A5; A6])
M3 = mixture(3, [0.3; 0.3; 0.4], [A4; A5; A6])

v1 = [ 1.; 0.; 0.; 0. ]
v2 = [ 0.; 1.; 0.; 0. ]
v3 = [ 0.; 0.; 1.; 0. ]
v4 = [ 0.; 0.; 0.; 1. ]
ms = [ M0; M1; M2; M3 ]
nb_images = 5


contour(X, X, density.(Ref(M0), XX), c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false)

contour(X, X, density.(Ref(M1), XX), c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false)

contour(X, X, density.(Ref(M2), XX), c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false)

contour(X, X, density.(Ref(M3), XX), c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false)

originals = []
for i in 1:nb_images
   for j in 1:nb_images
      tx = (i-1) / (nb_images - 1)
      ty = (j-1) / (nb_images - 1)
      tmp1 = (1 - tx) * v1 + tx * v2
      tmp2 = (1 - tx) * v3 + tx * v4
      λ = (1 - ty) * tmp1 + ty * tmp2
      @show λ
      @show typeof(λ)
      @show typeof(ms)
      bt = bar(λ, ms)
      P = contour(X, X, density.(Ref(bt), XX), c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false)
      display(P)

      push!(originals, contour(X, X, density.(Ref(bt), XX), c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false))

      println("Contour ($i, $j) done.")
   end
end

poriginals = plot(originals..., size = nb_images .* (512, 512), layout = grid(nb_images, nb_images))
savefig(poriginals, save_string * "originals.pdf")
save_object(save_string * "originals.jld2", save_originals)

pmodifieds = plot(modifieds..., size = nb_images .* (512, 512), layout = grid(nb_images, nb_images))
savefig(pmodifieds, save_string * "modifieds.pdf")
save_object(save_string * "modifieds.jld2", save_modifieds)

save_object(save_string * "distances.jld2", distances)




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
   P = contour(X, X, barydens[i], c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.), fill=:false, yformatter=_ -> "", xformatter=_ -> "", colorbar=:false)
   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_permantisym/bary_mw2_2d_contour$i.pdf")
end

# W2 barycenter
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
reg = 1e-4 #regularization parameter

baryw2 = [a0dens]
for t in T[2:end-1]
   @show t
   bt = ot.barycenter(A, M, reg, [1 - t; t], method="sinkhorn_log", numItermax=10000)
   push!(baryw2, bt)
end
push!(baryw2, a1dens)

for (i, t) in enumerate(T)
   P = contour(X, X, (reshape(baryw2[i], nx, nx))', c=:viridis, aspect_ratio=:equal, clim=(0.0, 10.0), colorbar=:false, yformatter=_ -> "", xformatter=_ -> "", zformatter=_ -> "")
   display(P)
   sleep(0.5)
   savefig("data_article_Dusson_Ehrlacher_Nouaime/section5_permantisym/bary_w2_2d_contour$i.pdf")
end