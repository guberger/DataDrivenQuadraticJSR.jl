include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

η_list = range(0, 0.1, length=20)
n = 2

f = open(string(@__DIR__, "/spherical_cap_results.txt"), "w")

for η in η_list
    s = DataDrivenQuadraticJSR.area_2sided_cap_inv(2*η, n)
    println(f, (η, sqrt(1/s^2 - 1)))
end

close(f)

print("")
sleep(0.1)