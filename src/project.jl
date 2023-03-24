using Turing
using Random
using DataFrames
using StatsFuns
using Plots, StatsPlots, Measures

using FASTX, BioSequences

#export randclasses, simpletwoclass, plotconfig, simplemulticlass, sample_mh_jump

include("sample.jl")
include("plotting.jl")
include("model.jl")
include("sampling.jl")
include("smf.jl")
include("plots.jl")


"""
    getprojectdir()
    function to determine the where the project directory is relative to current working directory which may be in a sub directory of the project dir
"""
function getprojectdir()
    d = pwd()
    if basename(d) == "notebooks"
        return dirname(d)
    else
        return d
    end
end