using LinearAlgebra
using SpecialFunctions
# ------------------------------------------------------------------------------
# Location-scatter distributions
# ------------------------------------------------------------------------------
struct loc_scat{T}
   m::Vector
   Σ::Matrix
   d::T
end

struct Slater end
struct Gauss end
struct Semicircle end
struct Gamma1d
   α::Float64
end

loc_scat(m::Float64, Σ::Float64, d::T) where {T} = loc_scat([m], [Σ;;], d)

slater(m::Vector, Σ::Matrix) = loc_scat(m, Σ, Slater())
slater(m::Float64, Σ::Float64) = slater([m], [Σ;;])
gauss(m::Vector, Σ::Matrix) = loc_scat(m, Σ, Gauss())
gauss(m::Float64, Σ::Float64) = gauss([m], [Σ;;])
semicircle(m::Vector, Σ::Matrix) = loc_scat(m, Σ, Semicircle())
semicircle(m::Float64, Σ::Float64) = semicircle([m], [Σ;;])
gamma1d(m::Float64, Σ::Float64, α) = loc_scat([m], [Σ;;], Gamma1d(α))


isvalid(g::loc_scat) = issymmetric(g.Σ) && isposdef(g.Σ) && (length(g.m) == size(g.Σ, 1))

function W2(g0::loc_scat{T}, g1::loc_scat{T}) where T
   # compute the square of quadratic Wasserstein distance
   Σ00 = sqrt(g0.Σ)
   Σ010 = sqrt(Σ00 * g1.Σ * Σ00)
   return norm(g0.m - g1.m)^2 + tr(g0.Σ + g1.Σ - 2 * Σ010)
end

function density(g::loc_scat{Gauss}, x)
   d = size(g.Σ, 1)
   y = (x - g.m)' * inv(g.Σ) * (x - g.m)
   Z = sqrt(det(2 * pi * g.Σ))
   return 1 / Z * exp(-1 / 2 * y)
end

function density(g::loc_scat{Slater}, x)
   d = size(g.Σ, 1)
   y = (x - g.m)' * inv(g.Σ) * (x - g.m)
   a = sqrt(d + 1)
   Z = 2 * (pi / (d + 1))^(d / 2) * gamma(d) / gamma(d / 2) * sqrt(det(g.Σ))
   return 1 / Z * exp(-a * sqrt(y))
end

function density(g::loc_scat{Semicircle}, x)
   d = size(g.Σ, 1)
   y = (x - g.m)' * inv(g.Σ) * (x - g.m)
   a = 1 / (d + 3)
   Z = (pi^((d + 1) / 2) * (d + 3)^(d / 2)) / (2 * gamma((d + 3) / 2)) * sqrt(det(g.Σ))
   if y <= 1 / a
      return 1 / Z * sqrt(1 - a * y)
   else
      return 0
   end
end

function density(g::loc_scat{Gamma1d}, x)
   d = size(g.Σ, 1)
   @assert (d == 1)
   α = g.d.α
   β = sqrt(α)
   y = 1 / sqrt(g.Σ[1, 1]) * (x[1] - g.m[1] + α / β)
   if y >= 0
      return 1 / sqrt(det(g.Σ)) * (β^α / gamma(α)) * y^(α - 1) * exp(-β * y)
   else
      return 0
   end
end

function bar(t, g0::loc_scat{T}, g1::loc_scat{T}) where T #barycenter between two atoms
   @assert(t <= 1 && t >= 0)
   sΣ1 = sqrt(g1.Σ)
   i, j = size(g0.Σ)
   C = sΣ1 * inv(sqrt(sΣ1 * g0.Σ * sΣ1)) * sΣ1
   Σt = ((1 - t) * Matrix(I, i, j) .+ t .* C) * g0.Σ * ((1 - t) * Matrix(I, i, j) .+ t .* C)
   mt = (1 - t) .* g0.m + t .* g1.m
   return loc_scat(mt, Σt, g0.d)
end

function barmulti(λs, gs::Vector{loc_scat{T}}; Niter=15) where T
   l = length(λs)
   @assert(l == length(gs))
   n = size(gs[1].Σ, 1)
   Σb = Matrix(I, n, n)

   for niter = 1:Niter
      Σb = sum(λs[j] * sqrt(sqrt(Σb) * gs[j].Σ * sqrt(Σb)) for j = 1:l)
   end
   mb = sum(λs[j] * gs[j].m for j in 1:l)
   gb = loc_scat(mb, Σb, gs[1].d)
   cost = sum(λs[j] * W2(gs[j], gb) for j in 1:l)
   return gb, real(cost)
end

