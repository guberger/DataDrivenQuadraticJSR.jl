using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

α = 0.5
S1 = @SMatrix [α (1-α) 0; (1-α) α 0; 0 0 0]
A1 = SMatrix{6,6}([S1 zeros(3, 3); zeros(3, 3) S1])

S2 = S1[SVector(2, 3, 1), SVector(2, 3, 1)]
A2 = SMatrix{6,6}([S2 zeros(3, 3); zeros(3, 3) S2])

S3 = S1[SVector(3, 1, 2), SVector(3, 1, 2)]
A3 = SMatrix{6,6}([S3 zeros(3, 3); zeros(3, 3) S3])

A_list = (A1, A2, A3)

nsample = 5000
x_list = [SVector{6}(randn(6)) for i = 1:nsample]
σ_list = rand((1, 2, 3), nsample)
y_list = [A_list[σ_list[i]]*x_list[i] for i = 1:nsample]
ηm_list = zeros(nsample)
ηa_list = zeros(nsample)

γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
    ηm_list, ηa_list, 0.0, 1.2; max_trace=nothing)

print("")
sleep(0.1)