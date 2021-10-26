using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

S1_edge = [
    1   0   1   0   0   0   0   1
    0   1   1   0   0   0   1   0
    0   0   1   0   0   1   0   0
    0   1   1   1   0   1   0   0
    1   0   0   0   1   1   0   0
    0   1   0   0   0   1   0   0
    0   0   1   0   0   0   1   1
    0   0   0   0   1   1   1   1]
T1_edge = [
    1 1 1
    1 1 0
    1 0 1]
A1_edge = SMatrix{11,11}([S1_edge zeros(8, 3); zeros(3, 8) T1_edge])
A1 = A1_edge ./ sum(A1_edge, dims=2)

S2_edge = [
    1   0   0   1   0   0   1   0
    0   1   0   1   0   0   0   1
    0   0   1   1   1   1   0   0
    1   1   0   1   0   1   0   0
    1   0   0   0   1   0   0   0
    0   1   0   0   1   1   0   0
    1   1   0   0   0   1   1   0
    1   1   0   0   1   1   0   1]
T2_edge = [
    1 1 0
    1 1 1
    0 1 1]
A2_edge = SMatrix{11,11}([S2_edge zeros(8, 3); zeros(3, 8) T2_edge])
A2 = A2_edge ./ sum(A2_edge, dims=2)

S3_edge = [
    1   0   1   0   1   0   0   0
    0   1   1   0   1   0   0   1
    1   0   1   0   1   1   0   1
    0   1   0   1   1   0   0   1
    0   1   0   1   1   0   0   1
    0   1   0   0   1   1   0   0
    0   1   1   0   0   0   1   0
    1   1   0   0   0   1   1   1]
T3_edge = [
    1 0 1
    0 1 1
    1 1 1]
A3_edge = SMatrix{11,11}([S3_edge zeros(8, 3); zeros(3, 8) T3_edge])
A3 = A3_edge ./ sum(A3_edge, dims=2)

A_list = [A1, A2, A3]

nsample = 2000
x_list = [SVector{11}(randn(11)) for i = 1:nsample]
σ_list = rand((1, 2, 3), nsample)
y_list = [A_list[σ_list[i]]*x_list[i] for i = 1:nsample]
ηm_list = zeros(nsample)
ηa_list = zeros(nsample)

# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
#     ηm_list, ηa_list, 0.0, 1.01; max_trace=nothing, ϵ=1e-3)

U1 = [
    1 0
    1 0
    1 0
    1 0
    1 0
    1 0
    1 0
    1 0
    0 1
    0 1
    0 1] ./ [sqrt(8) sqrt(3)]
U2 = nullspace(U1')
P = SMatrix{11,11}(U2*U2')

β = 0.01
ϵ = 0.054
sbar = 1/sqrt(1 + ϵ^2)
A = DataDrivenQuadraticJSR.area_2sided_cap(sbar, 11)
η = A/2
ϵbar = η/3
nsample = DataDrivenQuadraticJSR.single_monotone_bound_inv(β, ϵbar)
γ2_opt = 0.0

# for i = 1:nsample
#     x = SVector{11}(randn(11))
#     σ = rand((1, 2, 3))
#     y = A_list[σ]*x
#     local γ2_opt_temp = (y'*P*y) / (x'*P*x + eps(1.0))
#     global γ2_opt = max(γ2_opt,γ2_opt_temp)
#     # display((γ2_opt_temp, γ2_opt))
# end

γ_opt = sqrt(γ2_opt)


# x_list = [SVector{11}(randn(11)) for i = 1:nsample]
# σ_list = rand((1, 2, 3), nsample)
# y_list = [A_list[σ_list[i]]*x_list[i] for i = 1:nsample]
# ηm_list = zeros(nsample)
# ηa_list = zeros(nsample)
# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
#     ηm_list, ηa_list, 0.0, 0.0; P=P)

print("")
sleep(0.1)