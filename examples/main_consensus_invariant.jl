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

S3 = @SMatrix [α 0 0 (1-α); 0 α (1-α) 0; 0 (1-α) α 0; (1-α) 0 0 α]
A3 = SMatrix{8,8}([S3 zeros(4, 4); zeros(4, 4) S3])

A = (A1, A2, A3)

nsample = 3000
x_list = [SVector{8}(randn(8)) for i = 1:nsample]
σ_list = rand((1, 2, 3), nsample)
y_list = [A[σ_list[i]]*x_list[i] for i = 1:nsample]
ηm_list = zeros(nsample)
ηa_list = zeros(nsample)

γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
    ηm_list, ηa_list, 0.0, 1e2; max_trace=nothing)

print("")
sleep(0.1)