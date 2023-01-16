using Test
using QuadGK
using Cubature

include("../src/main.jl")

# Test parity 1d
SN = perm_group(2)
ga = parity()

m0 = [0.3]
Σ0 = [0.01;;]
a0 = slater(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = [0.0]
Σ1 = [0.2;;]
a1 = slater(m1, Σ1)
A1 = gsa(a1, SN, ga)

W2(A0, A1)
# 1d test on normalisation
f(x) = density(A0, [x])
R, e = quadgk(f, -20.0, 20.0, rtol=1e-8)
@test(abs.(R - 1.0) < 1e-5)

# Test parity 2d
SN = perm_group(2)
ga = parity()

m0 = [0.3; 0.6]
Σ0 = [0.01 0.0; 0.0 0.02]
a0 = slater(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = [0.0; -10.0]
Σ1 = [0.2 0.0; 0.0 0.5]
a1 = slater(m1, Σ1)
A1 = gsa(a1, SN, ga)

W2(A0, A1)
# 2d test on normalisation
f(x) = density(A0, [x[1]; x[2]])
R, e = hcubature(f, [-20.0, -20.0], [20.0, 20.0], abstol=1e-8)
@test(abs.(R - 1.0) < 1e-6)

# ------------------------------------------
# 
#  Test sym group
# 
# ------------------------------------------
SN = perm_group(2)
ga = perm_sym()

m0 = [0.3; 0.6]
Σ0 = [0.01 0.0; 0.0 0.02]
a0 = slater(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = [0.0; 1.0]
Σ1 = [0.2 0.0; 0.0 0.5]
a1 = slater(m1, Σ1)
A1 = gsa(a1, SN, ga)

W2(A0, A1)
# 2d test on normalisation
f(x) = density(A0, [x[1]; x[2]])
R, e = hcubature(f, [-20.0, -20.0], [20.0, 20.0], abstol=1e-8)
@test(abs.(R - 1.0) < 1e-6)


# ------------------------------------------
# 
#  Test sym group square
# 
# ------------------------------------------
SN = perm_group(2)
ga = perm_sym_square()

m0 = [0.3; 0.6]
Σ0 = [0.01 0.0; 0.0 0.02]
a0 = gauss(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = [0.0; 1.0]
Σ1 = [0.2 0.0; 0.0 0.5]
a1 = gauss(m1, Σ1)
A1 = gsa(a1, SN, ga)

W2(A0, A1)
# 2d test on normalisation
f(x) = density(A0, [x[1]; x[2]])
R, e = hcubature(f, [-40.0, -40.0], [40.0, 40.0], abstol=1e-8)
@test(abs.(R - 1.0) < 1e-4)


# ------------------------------------------
# 
#  Test sym group antisym
# 
# ------------------------------------------
SN = perm_group(2)
ga = perm_antisym()

m0 = [0.4; 0.6]
Σ0 = [0.01 0.0; 0.0 0.02]
a0 = gauss(m0, Σ0)
A0 = gsa(a0, SN, ga)

m1 = [0.0; 1.0]
Σ1 = [0.2 0.0; 0.0 0.5]
a1 = gauss(m1, Σ1)
A1 = gsa(a1, SN, ga)

W2(A0, A1)
# 2d test on normalisation
f(x) = density(A0, [x[1]; x[2]])
R, e = hcubature(f, [-21.0, -21.0], [21.0, 21.0], abstol=1e-8)
# @test(abs.(R - 1.0) < 1e-6)


# ------------------------------------------
# 
#  Test rotation group
# 
# ------------------------------------------
# Test rotation 2d
S0 = rot_group()
ga = rot_sym()

m0 = [0.3; 0.6]
Σ0 = [0.01 0.0; 0.0 0.02]
a0 = slater(m0, Σ0)
A0 = gsa(a0, S0, ga)

m1 = [0.0; -10.0]
Σ1 = [0.2 0.0; 0.0 0.5]
a1 = slater(m1, Σ1)
A1 = gsa(a1, S0, ga)

θ = 0.5
Q = [cos(θ) sin(θ); -sin(θ) cos(θ)]

group_action(ga, A1.a, Q)

W2(A0, A1)
# 2d test on normalisation
f(x) = density(A0, [x[1]; x[2]])
R, e = hcubature(f, [-40.0, -40.0], [40.0, 40.0], abstol=1e-8)
@test(abs.(R - 1.0) < 1e-6)