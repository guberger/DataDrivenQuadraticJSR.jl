include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

β_list = 0.01
ϵ_list = 0.5 .^ (1:10)
k = 1

f = open(string(@__DIR__, "/proba_bounds_inverse_results.txt"), "w")

for ϵ in ϵ_list
    N = DataDrivenQuadraticJSR.varying_bound_2_inv(β, ϵ, k)
    println(f, (N, ϵ))
end

close(f)

print("")
sleep(0.1)