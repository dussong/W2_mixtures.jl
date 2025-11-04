using Test

include("../src/main.jl")

# Test validity of mixtures

# define several gaussian distribution
G0 = gauss([3.0; 0.0], [4.0 0.1; 0.1 3.0])
G1 = gauss([3.0; 3.0], [4.0 0.0; 0.0 2.0])
M = mixture(2, [0.3; 0.7], [G0; G1])
@test(isvalid(M))

M1 = mixture(4, [0.3; 0.7], [G0; G1])
@test(!isvalid(M1))
M2 = mixture(2, [0.3; 0.8], [G0; G1])
@test(!isvalid(M2))

# test barycentre pour mixture with one element
M0 = mixture(1, [1.0], [G0])
M1 = mixture(1, [1.0], [G1])
t = rand()
mbar = bar(t, M0, M1)
gbar = bar(t, G0, G1)
@test(mbar.K == 1)
@test(mbar.p[1] == 1.0)
@test(mbar.K == 1)
@test(mbar.at[1].m == gbar.m)
@test(mbar.at[1].Σ == gbar.Σ)

# test barycenter multimarginal 2 mixtures - 1d
# test barycentre multimarginal with two mixtures
G0 = gauss(3.0, 4.0)
G1 = gauss(2.0, 2.5)
G2 = gauss(6.0, 3.0)
G3 = gauss(-4.0, 0.3)

M0 = mixture(1, [1.0], [G0])
M1 = mixture(1, [1.0], [G2])

t = rand()
bar1 = bar(t, M0, M1)
bar2 = bar([1 - t; t], [M0; M1])
x = [[i] for i in range(-15, 15, length=1000)]
d1 = density.(Ref(bar1), x)
d2 = density.(Ref(bar2), x)
@test(norm(d1 - d2, Inf) < 1e-5)

M0 = mixture(2, [0.2; 0.8], [G0; G1])
M1 = mixture(2, [0.5; 0.5], [G2; G3])

t = rand()
bar1 = bar(t, M0, M1)
bar2 = bar([1 - t; t], [M0; M1])
x = [[i] for i in range(-15, 15, length=1000)]
d1 = density.(Ref(bar1), x)
d2 = density.(Ref(bar2), x)
@test(norm(d1 - d2, Inf) < 1e-5)


# test barycentre multimarginal with two mixtures - 2d
G0 = gauss([3.0; 0.0], [4.0 0.1; 0.1 3.0])
G1 = gauss([3.0; 3.0], [4.0 0.0; 0.0 2.0])
G2 = gauss([5.0; -2.0], [4.0 0.1; 0.1 5.5])
G3 = gauss([2.0; 3.0], [6.0 0.2; 0.2 2.0])
M0 = mixture(2, [0.2; 0.8], [G0; G1])
M1 = mixture(2, [0.5; 0.5], [G2; G3])

t = rand()
bar1 = bar(t, M0, M1)
bar2 = bar([1 - t; t], [M0; M1])

nx = 10
xmin, xmax = -15.0, 15.0
X = range(xmin, xmax, length=nx) |> collect;
xy = [[X[i], X[j]] for i in 1:nx for j = 1:nx]
d1 = density.(Ref(bar1), xy)
d2 = density.(Ref(bar2), xy)
@test(norm(d1 - d2, Inf) < 1e-6)