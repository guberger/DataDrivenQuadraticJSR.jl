using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

A1_edge = [
    1   -1  0   1
    -1  1   -1  0
    1   0   1   1
    0   -1  1   1]
A1 = SMatrix{4,4}(A1_edge ./ sum(abs, A1_edge, dims=2))

A2_edge = [
    1   -1  1   0
    -1  1   -1  -1
    0   -1  1   1
    1   0   1   1]
A2 = SMatrix{4,4}(A2_edge ./ sum(abs, A2_edge, dims=2))

A3_edge = [
    1   0   1   1
    0   1   -1  -1
    1   -1  1   1
    1   -1  0   1]
A3 = SMatrix{4,4}(A3_edge ./ sum(abs, A3_edge, dims=2))

A_list = [A1, A2, A3]

# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic_white(
#     A_list, 0.0, 1.2; max_trace=nothing)

nsample = 2000
x_list = [SVector{4}(randn(4)) for i = 1:nsample]
σ_list = rand((1, 2, 3), nsample)
y_list = [A_list[σ_list[i]]*x_list[i] for i = 1:nsample]
ηm_list = zeros(nsample)
ηa_list = zeros(nsample)

# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
#     ηm_list, ηa_list, 0.0, 1.01; max_trace=nothing, ϵ=1e-3)

U1 = ([1, -1, 1, 1]/sqrt(4))[:, :]
U2 = nullspace(U1')
P = SMatrix{4,4}(U2*U2')

β = 0.01
ϵ = 0.02
sbar = 1/sqrt(1 + ϵ^2)
A = DataDrivenQuadraticJSR.area_2sided_cap(sbar, 4)
η = A/2
ϵbar = η/3
nsample = DataDrivenQuadraticJSR.single_monotone_bound_inv(β, ϵbar)
γ2_opt = 0.0
counter = 0

for i = 1:nsample
    global counter +=1
    x = SVector{4}(randn(4))
    σ = rand((1, 2, 3))
    y = A_list[σ]*x
    local γ2_opt_temp = (y'*P*y) / (x'*P*x + eps(1.0))
    global γ2_opt = max(γ2_opt,γ2_opt_temp)
    if mod(counter - 1, 1_000_000) == 0
        display((counter, counter/nsample, nsample, γ2_opt, sqrt(γ2_opt)))
    end
end

γ_opt = sqrt(γ2_opt)

# x_list = [SVector{4}(randn(4)) for i = 1:nsample]
# σ_list = rand((1, 2, 3), nsample)
# y_list = [A_list[σ_list[i]]*x_list[i] for i = 1:nsample]
# ηm_list = zeros(nsample)
# ηa_list = zeros(nsample)
# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
#     ηm_list, ηa_list, 0.0, 0.0; P=P)

print("")
sleep(0.1)