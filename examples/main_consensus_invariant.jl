using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

α = 0.5
S1 = @SMatrix [α (1-α) 0 0; (1-α) α 0 0; 0 0 α (1-α); 0 0 (1-α) α]
A1 = SMatrix{8,8}([S1 zeros(4, 4); zeros(4, 4) S1])

S2 = @SMatrix [α 0 (1-α) 0; 0 α 0 (1-α); (1-α) 0 α 0; 0 (1-α) 0 α]
A2 = SMatrix{8,8}([S2 zeros(4, 4); zeros(4, 4) S2])

print("")
sleep(0.1)