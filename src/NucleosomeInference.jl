module NucleosomeInference

# Write your package code here.
using Turing
using Random
using DataFrames
using Plots, StatsPlots, Measures

export randclasses, simpletwoclass, plotconfig, simplemulticlass, sample_mh_jump

include("sample.jl")
include("plotting.jl")
include("model.jl")
include("sampling.jl")



end
