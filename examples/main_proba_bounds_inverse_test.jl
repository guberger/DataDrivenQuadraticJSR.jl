include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

β = 0.001
ϵ = 0.01
k = 5

N = DataDrivenQuadraticJSR.varying_bound_2_inv(β, ϵ, k)
display(N)
ϵ_int_1 = DataDrivenQuadraticJSR.varying_bound_2(N - 1, k, β)
ϵ_int_2 = DataDrivenQuadraticJSR.varying_bound_2(N, k, β)

display(N)
display(ϵ_int_1)
display(ϵ_int_2)

print("")
sleep(0.1)