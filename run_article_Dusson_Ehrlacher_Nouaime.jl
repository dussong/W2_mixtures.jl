# Loading packages
#-----------------
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# run test cases section 4
println("Test on 1D Slater distribution")
include("./article_Dusson_Ehrlacher_Nouaime/section_4/slater1d.jl")

println("Test on 1D gamma distribution")
include("./article_Dusson_Ehrlacher_Nouaime/section_4/gamma1d.jl")

println("Test on 1D Wigner semicircle distribution")
include("./article_Dusson_Ehrlacher_Nouaime/section_4/semicircle1d.jl")

println("Test on 2D Slater distribution")
include("./article_Dusson_Ehrlacher_Nouaime/section_4/slater2d.jl")

println("Test on 2D Wigner semicircle distribution")
include("./article_Dusson_Ehrlacher_Nouaime/section_4/semicircle2d.jl")

# run test cases section 5
println("Test on 1D even distributions")
include("./article_Dusson_Ehrlacher_Nouaime/section_5/parity1d.jl")

println("Test on 2D permutation symmetric distributions")
include("./article_Dusson_Ehrlacher_Nouaime/section_5/perm_sym2d.jl")

println("Test on 2D permutation symmetric - from antisymmetric distributions")
include("./article_Dusson_Ehrlacher_Nouaime/section_5/perm_antisym2d.jl")

println("Test on 2D rotation symmetric distributions")
include("./article_Dusson_Ehrlacher_Nouaime/section_5/rot_sym2d.jl")