using JuMP
using GLPK

struct mixture{T}
   K::Int
   p::Vector
   at::Vector{T}
end

function isvalid(M::mixture)
   out = true
   k = M.K
   if (length(M.p) != k) || (length(M.at) != k) || (abs.(sum(M.p) - 1.0) > 1e-8)
      out = false
   else
      for i in 1:k
         if !(isvalid((M.at)[i]))
            out = false
         end
      end
   end
   return out
end

function density(M::mixture, x)
   return sum(M.p[k] * density((M.at)[k], x) for k in 1:M.K)
end

function w_kl(M0, M1) #solves (3.1) and find the values of the wstar
   function f(M0, M1)
      c = 0
      for k in 1:M0.K
         for l in 1:M1.K
            c += (w[k, l]) * W2(M0.at[k], M1.at[l])
         end
      end
      return c
   end
   model = Model(GLPK.Optimizer)
   @variable(model, w[1:M0.K, 1:M1.K])
   @constraint(model, [k = 1:M0.K], sum(w[k, :]) == M0.p[k])
   @constraint(model, [l = 1:M1.K], sum(w[:, l]) == M1.p[l])
   @constraint(model, [k = 1:M0.K, l = 1:M1.K], 0 <= w[k, l])
   @objective(model, Min, f(M0, M1))
   optimize!(model)
   return value.(w)
end

function bar(t, M0::mixture{T}, M1::mixture{T}) where {T} 
   #barycenter between two mixtures for t in [0,1]
   @assert(t <= 1 && t >= 0)
   K, p, at = 0, Float64[], T[]
   wkl = w_kl(M0, M1)
   for k = 1:M0.K
      for l = 1:M1.K
         if wkl[k, l] != 0.0
            gt = bar(t, M0.at[k], M1.at[l])
            K += 1
            push!(p, wkl[k, l])
            push!(at, gt)
         end
      end
   end
   return mixture(K, p, at)
end

function cost_matrix_mix(gmm::Vector{mixture{T}}) where {T}
   nMarginal = length(gmm) # number of marginals
   KK = [gmm[k].K for k in 1:nMarginal] #number of elements in each mixture
   C = Array{Float64}(undef, KK...)
   λs = 1 / nMarginal * ones(1, nMarginal)
   for k in CartesianIndices(C)
      _, cost = barmulti(λs, [gmm[j].at[k[j]] for j in 1:nMarginal])
      C[k] = cost
   end
   return C
end

function wmulti(gmm::Vector{mixture{T}}) where {T}
   n = length(gmm)
   KK = [gmm[k].K for k in 1:n]
   nb_unknowns = prod(KK)
   Cmatrix = cost_matrix_mix(gmm)
   Cflat = reshape(Cmatrix, nb_unknowns)
   Indices = reshape(CartesianIndices(Cmatrix), nb_unknowns)

   f() = sum(w .* Cflat) #function to optimize
   # Model variables
   model = Model(GLPK.Optimizer)
   @variable(model, w[1:nb_unknowns])
   # Build constraint indices (over which the sum happens)
   cst = [[findall(x -> x == k, [Indices[l][j] for l in 1:nb_unknowns])
           for k in 1:KK[j]]
          for j in 1:n]
   # Constraints
   @constraint(model, [k = 1:nb_unknowns], 0 <= w[k])
   for j = 1:n
      @constraint(model, [k = 1:KK[j]], sum(w[cst[j][k]]) == gmm[j].p[k])
   end
   @objective(model, Min, f())
   optimize!(model)
   return value.(w), Indices
end

function bar(λs::Vector, gmm::Vector{mixture{T}}) where {T}
   # barycenter between n mixtures
   @assert(abs.(sum(λs) - 1.0) < 1e-10)
   @assert(length(λs) == length(gmm))
   n = length(λs)
   w, indices = wmulti(gmm)
   k = 0
   Pi = Float64[]
   gtot = T[]
   for (ni, i) in enumerate(indices)
      if abs.(w[ni]) > 1e-13
         k += 1
         gs = [gmm[j].at[i[j]] for j in 1:n]
         gnew, _ = barmulti(λs, gs)
         push!(gtot, gnew)
         push!(Pi, w[ni])
      end
   end
   return mixture(k, Pi, gtot)
end