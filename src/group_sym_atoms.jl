using LinearAlgebra
using Combinatorics
using Optim
# ------------------------------------------------------------------------------
# Location-scatter distributions
# ------------------------------------------------------------------------------
struct gsa{T1,T2,T3}
   a::T1 #atom type
   G::T2 #group
   ga::T3 #group action
end
isvalid(ga::gsa) = isvalid(ga.a)

struct perm_group
   n::Int
   G::Vector{Vector{Int}} #vector of group elements (permutations)
end
perm_group(n::Int) = perm_group(n, collect(permutations(1:n)))

struct rot_group end

struct parity end
struct perm_sym end
struct perm_sym_square end
struct perm_antisym end
struct rot_sym end

function group_action(::perm_sym, a::loc_scat, g)
   ms = a.m[g]
   Σs = a.Σ[g, g]
   return loc_scat(ms, Σs, a.d)
end

function group_action(::perm_sym_square, a::loc_scat, g)
   ms = a.m[g]
   Σs = a.Σ[g, g]
   return loc_scat(ms, Σs, a.d)
end

function group_action(::perm_antisym, a::loc_scat, g)
   ms = a.m[g]
   Σs = a.Σ[g, g]
   return loc_scat(ms, Σs, a.d)
end

function group_action(::parity, a::loc_scat, g) #TODO does not work for gamma1d
   if g == [2, 1]
      ms = -a.m
   else
      ms = a.m
   end
   return loc_scat(ms, a.Σ, a.d)
end

function group_action(::rot_sym, a::loc_scat, g)
   ms = g' * a.m
   Σs = g' * a.Σ * g
   return loc_scat(ms, Σs, a.d)
end

# function _W2(A0::gsa{T1,rot_group,rot_sym}, A1::gsa{T1,rot_group,rot_sym}) where {T1,T3}
#    @assert length(A0.a.m) == 2
#    # solve optim problem
# end

function _W2(A0::gsa{T1,rot_group,rot_sym}, A1::gsa{T1,rot_group,rot_sym}) where {T1}
   function g(θ)
      Q = [cos(θ) sin(θ); -sin(θ) cos(θ)]
      W2(A0.a, group_action(ga, A1.a, Q))
   end
   res = optimize(g, 0.0, 2 * pi)
   θopt = res.minimizer
   W2value = res.minimum
   Qopt = [cos(θopt) sin(θopt); -sin(θopt) cos(θopt)]
   return W2value, Qopt
end

function W2(A0::gsa{T1,rot_group,rot_sym}, A1::gsa{T1,rot_group,rot_sym}) where {T1}
   return _W2(A0, A1)[1]
end

function _W2(A0::gsa{T1,T2,T3}, A1::gsa{T1,T2,T3}) where {T1,T2,T3}
   d, ind = findmin(x -> W2(A0.a, group_action(A1.ga, A1.a, x)), A0.G.G)
   return d, A0.G.G[ind]
end

W2(A0::gsa{T1,T2,T3}, A1::gsa{T1,T2,T3}) where {T1,T2,T3} = _W2(A0, A1)[1]

function density(A::gsa{T1,T2,parity}, x) where {T1,T2}
   N = length(A.G.G)
   return 1.0 / N * sum(density(group_action(A.ga, A.a, g), x) for g in A.G.G)
end

function density(A::gsa{T1,T2,perm_sym}, x) where {T1,T2}
   N = length(A.G.G)
   return 1.0 / N * sum(density(group_action(A.ga, A.a, g), x) for g in A.G.G)
end

function density(A::gsa{loc_scat{Gauss},perm_group,perm_sym_square}, x)
   N = length(A.G.G)
   @assert A.G.n == 2
   a1 = group_action(A.ga, A.a, A.G.G[1])
   a2 = group_action(A.ga, A.a, A.G.G[2])
   c = (1 / sqrt(det(inv(a1.Σ) / pi))
        + 1 / sqrt(det(inv(a2.Σ) / pi))
        + 2 * (exp(-1 / 2 * (a1.m - a2.m)' * inv(a1.Σ + a2.Σ) * (a1.m - a2.m)) / sqrt(det((inv(a1.Σ) + inv(a2.Σ)) / (2 * pi))))
   )
   return 1 / c * (exp(-1 / 2 * (x - a1.m)' * inv(a1.Σ) * (x - a1.m))
                   +
                   exp(-1 / 2 * (x - a2.m)' * inv(a2.Σ) * (x - a2.m)))^2
end

function density(A::gsa{loc_scat{Gauss},perm_group,perm_antisym}, x)
   N = length(A.G.G)
   @assert A.G.n == 2
   a1 = group_action(A.ga, A.a, A.G.G[1])
   a2 = group_action(A.ga, A.a, A.G.G[2])
   c = (1 / sqrt(det(inv(a1.Σ) / pi))
        +
        1 / sqrt(det(inv(a2.Σ) / pi))
        -
        2 * (exp(-1 / 2 * (a1.m - a2.m)' * inv(a1.Σ + a2.Σ) * (a1.m - a2.m)) / sqrt(det((inv(a1.Σ) + inv(a2.Σ)) / (2 * pi))))
   )
   return 1 / c * (exp(-1 / 2 * (x - a1.m)' * inv(a1.Σ) * (x - a1.m))
                   -
                   exp(-1 / 2 * (x - a2.m)' * inv(a2.Σ) * (x - a2.m)))^2
end

function density(A::gsa{T1,rot_group,rot_sym}, x) where {T1}
   # density evaluated with numerical sum
   N = 200
   θgrid = range(0.0, 2 * pi, length=N + 1)
   return 1.0 / N * sum(density(group_action(A.ga, A.a, [cos(θ) sin(θ); -sin(θ) cos(θ)]), x) for θ in θgrid[1:N])
end

function bar(t, A0::gsa{T1,T2,T3}, A1::gsa{T1,T2,T3}) where {T1,T2,T3} #barycenter between two sym atoms
   @assert(t <= 1 && t >= 0)
   _, g = _W2(A0, A1)
   return gsa(bar(t, A0.a, group_action(A0.ga, A1.a, g)), A0.G, A0.ga)
end

# hack for 4 functions and antisymmetric group action
function _mW2(A0::gsa{T1,T2,T3}, A1::gsa{T1,T2,T3}, A2::gsa{T1,T2,T3}, A3::gsa{T1,T2,T3}, λs) where {T1,T2,T3}
   gmin = A0
   mw2g = 1e10
   for g1 in A0.G.G
      for g2 in A0.G.G
         for g3 in A0.G.G
            v = [A0;group_action(A1.ga, A1.a, g1);
                                  group_action(A2.ga, A2.a, g2);
                                  group_action(A3.ga, A3.a, g3)]
            gb, c = barmulti(λs, [A0.a;group_action(A1.ga, A1.a, g1);
                                  group_action(A2.ga, A2.a, g2);
                                  group_action(A3.ga, A3.a, g3)])
            if c < mw2g
               mw2g = c 
               gmin = gb
            end
         end
      end
   end
   return gmin, mw2g
end

function barmulti(λs, gs::Vector{gsa{T1, T2, T3}}) where {T1, T2, T3}
   gmin, c = _mW2(gs[1], gs[2], gs[3], gs[4], λs)
   return gsa(gmin, perm_group(2), perm_antisym()), c
end

# function barmulti(λs, gs::Vector{loc_scat{T}}; Niter=15) where {T}
#    l = length(λs)
#    @show gs
#    @assert(l == length(gs))
#    n = size(gs[1].Σ, 1)
#    Σb = Matrix(I, n, n)

#    for niter = 1:Niter
#       Σb = sum(λs[j] * sqrt(sqrt(Σb) * gs[j].Σ * sqrt(Σb)) for j = 1:l)
#    end
#    mb = sum(λs[j] * gs[j].m for j in 1:l)
#    gb = loc_scat(mb, Σb, gs[1].d)
#    cost = sum(λs[j] * W2(gs[j], gb) for j in 1:l)
#    return gb, real(cost)
# end

