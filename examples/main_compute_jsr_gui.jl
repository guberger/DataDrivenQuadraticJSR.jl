using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

# α = 0.25*2*π
# P = [8.0 0.0; 0.0 1.0]
# A = P*0.5*(@SMatrix [cos(α) -sin(α); sin(α) cos(α)])/P
A = @SMatrix [0 1; 1 0]
A1 = @SMatrix [2 1; 0.1 1]
A2 = @SMatrix [2 1; -0.1 1]

my_gui = DataDrivenQuadraticJSR.new_gui_compute_jsr((A1, A2))

print("")
sleep(0.1)