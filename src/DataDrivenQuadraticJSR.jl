module DataDrivenQuadraticJSR

using LinearAlgebra
using StaticArrays
using Printf
using JuMP
using MosekTools
using PyPlot
using Distributions

include("proba_bounds/gui.jl")
include("proba_bounds/bounds.jl")

include("compute_jsr/gui.jl")
include("compute_jsr/optimize.jl")

include("spherical_measure.jl")

end # module
