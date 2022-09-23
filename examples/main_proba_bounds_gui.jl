include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

my_gui = DataDrivenQuadraticJSR.new_gui_proba_bounds()

print("")
sleep(0.1)