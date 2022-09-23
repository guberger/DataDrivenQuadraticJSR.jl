using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

S1_edge = [
    1 1 0 0 1
    0 1 1 1 0
    1 0 1 1 0
    0 0 1 1 1
    1 1 0 0 1]
T1_edge = [
    1 1 1
    1 1 0
    1 0 1]
A1_edge = SMatrix{8,8}([S1_edge zeros(5, 3); zeros(3, 5) T1_edge])
A1 = A1_edge ./ sum(A1_edge, dims=2)

S2_edge = [
    1 0 0 1 1
    1 1 1 0 0
    1 1 1 0 0
    0 1 1 1 0
    1 0 0 1 1]
T2_edge = [
    1 1 0
    1 1 1
    0 1 1]
A2_edge = SMatrix{8,8}([S2_edge zeros(5, 3); zeros(3, 5) T2_edge])
A2 = A2_edge ./ sum(A2_edge, dims=2)

S3_edge = [
    1 0 1 1 0
    0 1 0 1 1
    1 0 1 0 1
    0 1 0 1 1
    1 0 1 0 1]
T3_edge = [
    1 0 1
    0 1 1
    1 1 1]
A3_edge = SMatrix{8,8}([S3_edge zeros(5, 3); zeros(3, 5) T3_edge])
A3 = A3_edge ./ sum(A3_edge, dims=2)

A_list = [A1, A2, A3]

nsample = 2000
x_list = [SVector{8}(randn(8)) for i = 1:nsample]
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
    0 1
    0 1
    0 1] ./ [sqrt(5) sqrt(3)]
U2 = nullspace(U1')
P = SMatrix{8,8}(U2*U2')

β = 0.01
ϵ = 0.122
sbar = 1/sqrt(1 + ϵ^2)
A = DataDrivenQuadraticJSR.area_2sided_cap(sbar, 8)
η = A/2
ϵbar = η/3
nsample = DataDrivenQuadraticJSR.single_monotone_bound_inv(β, ϵbar)
γ2_opt = 0.0
counter = 0

for i = 1:nsample
    global counter +=1
    x = SVector{8}(randn(8))
    σ = rand((1, 2, 3))
    y = A_list[σ]*x
    local γ2_opt_temp = (y'*P*y) / (x'*P*x + eps(1.0))
    global γ2_opt = max(γ2_opt,γ2_opt_temp)
    if mod(counter - 1, 1_000_000) == 0
        display((counter, counter/nsample, nsample, γ2_opt, sqrt(γ2_opt)))
    end
end

γ_opt = sqrt(γ2_opt)

print("")
sleep(0.1)