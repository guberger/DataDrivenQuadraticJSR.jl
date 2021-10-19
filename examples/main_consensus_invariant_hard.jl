using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

EYE = SMatrix{4,4}(I)
α = 0.5
S1_1 = @SMatrix [0 1 1 1; 1 0 0 1; 1 0 0 1; 1 1 1 0]
S1_2 = S1_1 + EYE
S1 = (S1_2 ./ sum(S1_2, dims=1))'
A1 = SMatrix{8,8}([S1 zeros(4, 4); zeros(4, 4) S1])

S2_1 = @SMatrix [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0]
S2_2 = S2_1 + EYE
S2 = (S2_2 ./ sum(S2_2, dims=1))'
A2 = SMatrix{8,8}([S2 zeros(4, 4); zeros(4, 4) S2])

S3_1 = @SMatrix [0 0 1 1; 0 0 1 1; 1 1 0 1; 1 1 1 0]
S3_2 = S3_1 + EYE
S3 = (S3_2 ./ sum(S3_2, dims=1))'
A3 = SMatrix{8,8}([S3 zeros(4, 4); zeros(4, 4) S3])

A_list = [A1, A2, A3]

# U1 = [1 0; 1 0; 1 0; 1 0; 0 1; 0 1; 0 1; 0 1]
# G = nullspace(U1')
# U = [U1 G]
# Ap_list = [U'*A*U for A in A_list]
# display(Ap_list[1])
# App_list = [SMatrix{6,6}(A[3:end, 3:end]) for A in Ap_list]
# display(App_list[1])

# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic_white(
#     A_list, 0.0, 1.2; max_trace=nothing)
# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic_white(
#     App_list, 0.0, 1.2; max_trace=nothing)

# U = nullspace([1, 1, 1, 1][:, :]')
# Appp_list = [SMatrix{3,3}(U'*S*U) for S in (S1, S2, S3)]
# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic_white(
#     Appp_list, 0.0, 1.2; max_trace=nothing)

nsample = 20000
x_list = [SVector{8}(randn(8)) for i = 1:nsample]
σ_list = rand((1, 2, 3), nsample)
y_list = [A_list[σ_list[i]]*x_list[i] for i = 1:nsample]
ηm_list = zeros(nsample)
ηa_list = zeros(nsample)

# γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
#     ηm_list, ηa_list, 0.5, 0.75; max_trace=nothing, ϵ=1e-3)

U1 = [1 0; 1 0; 1 0; 1 0; 0 1; 0 1; 0 1; 0 1]
U2 = nullspace(U1')
P = SMatrix{8,8}(U2*U2')
γ_opt, P_opt = DataDrivenQuadraticJSR.jsr_quadratic(x_list, y_list,
    ηm_list, ηa_list, 0.0, 0.0; P=P)

print("")
sleep(0.1)