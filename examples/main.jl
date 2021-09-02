using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

α = 0.25*2*π
P = [8.0 0.0; 0.0 1.0]
A = P*0.5*(@SMatrix [cos(α) -sin(α); sin(α) cos(α)])/P
A = @SMatrix [0 1; 1 0]

GUI_ELEMS = new_gui(A)

slider_γ = GUI_ELEMS[4]
slider_tr = GUI_ELEMS[5]
the_P_opt = GUI_ELEMS[8]

print("")
sleep(0.1)