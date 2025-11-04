# W2_mixtures.jl

This is the code used to generate the figures of the article entitled *A Wasserstein-type interpolation for generic mixture models, including 
location-scatter and 
group invariant measures* by Geneviève Dusson, Virginie Ehrlacher, and Nathalie Nouaime.

## Requirements

To run the code, you need Julia 1.7 with the following packages:

LinearAlgebra, Combinatorics, Optim, Cubature, GLPK, JuMP, Plots, PyCall, QuadGK, SpecialFunctions, Test

You also need the Python Optimal Transport library [POT](https://pythonot.github.io/) to be installed on the Python distribution called with PyCall.

## Run the numerical tests

To run the simulations, simply run 

`run_article_Dusson_Ehrlacher_Nouaime.jl`

All the figures from the article should be generated in the corresponding subfolders of the folder `data_article_Dusson_Ehrlacher_Nouaime`.