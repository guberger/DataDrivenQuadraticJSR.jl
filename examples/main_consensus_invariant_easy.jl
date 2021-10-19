using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

α = 0.5
S = @SMatrix [α (1-α); (1-α) α]
A1 = SMatrix{4,4}([S zeros(2, 2); zeros(2, 2) zeros(2, 2)])
A2 = SMatrix{4,4}([zeros(2, 2) zeros(2, 2); zeros(2, 2) S])
A3 = SMatrix{4,4}([S zeros(2, 2); zeros(2, 2) S])

A_list = (A1, A2, A3)

nsample = 5000
x_list = [SVector{4}(randn(4)) for i = 1:nsample]
σ_list = rand((1, 2, 3), nsample)
y_list = [A_list[σ_list[i]]*x_list[i] for i = 1:nsample]
ηm_list = zeros(nsample)
ηa_list = zeros(nsample)

γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
    ηm_list, ηa_list, 0.0, 1.2; max_trace=nothing)

print("")
sleep(0.1)