function jsr_quadratic(x_list::Vector{SVector{N,T}}, y_list::Vector{SVector{N,T}},
        ηm_list::Vector{T}, ηa_list::Vector{T},
        γ_min::Float64, γ_max::Float64;
        max_trace::Float64 = 1e5, ϵ::Float64 = 1e-5) where {N,T}
    EYE = SMatrix{N,N}(I)
    γ_up = γ_max
    γ_do = γ_min
    γ_opt = 0.0
    P_opt = zero(SMatrix{N,N,T})

    while γ_up - γ_do > ϵ
        γ = (γ_up + γ_do)/2
        γ2 = γ^2
        @printf("γ = %f: ", γ)

        model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))
        P = @variable(model, [1:N, 1:N], Symmetric)

        @constraint(model, Symmetric(P - EYE) in PSDCone())
        @constraint(model, tr(P) <= max_trace)

        for (x, y, ηm, ηa) in zip(x_list, y_list, ηm_list, ηa_list)
            aff = y'*P*y - γ2*(x'*P*x + ηa)*(1.0 + ηm)
            @constraint(model, aff <= 0.0)
        end

        @objective(model, Min, sum(P[:].^2))

        optimize!(model)

        @printf("%s - %s - %s\n",
            termination_status(model), primal_status(model), dual_status(model))

        if Int(primal_status(model)) == 1
            γ_up = γ
            P_opt = SMatrix{N,N,T}(value.(P))
            γ_opt = γ
        else
            γ_do = γ
        end
    end

    @printf("Interval for γ: %f - %f\n", γ_do, γ_up)
    return γ_opt, P_opt
end
