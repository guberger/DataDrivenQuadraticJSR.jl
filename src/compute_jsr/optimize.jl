function jsr_quadratic(x_list::Vector{SVector{N,T}}, y_list::Vector{SVector{N,T}},
        ηm_list::Vector{T}, ηa_list::Vector{T},
        γ_min::Float64, γ_max::Float64;
        max_trace::Union{Nothing,Float64}=1e5, ϵ::Float64=1e-5) where {N,T}
    γ_up = γ_max
    γ_do = γ_min
    γ_opt = 0.0
    P_opt = zero(SMatrix{N,N,T})

    while γ_up - γ_do > ϵ
        γ = (γ_up + γ_do)/2
        @printf("γ = %f: ", γ)

        if !isnothing(max_trace)
            model, P = _optim_model_regular(Val(N),
                x_list, y_list, ηm_list, ηa_list, γ, max_trace)
        else
            model, P = _optim_model_singular(Val(N),
                x_list, y_list, ηm_list, ηa_list, γ)
        end

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

function _optim_model_regular(::Val{N},
        x_list, y_list, ηm_list, ηa_list, γ, max_trace) where {N}
    γ2 = γ^2
    EYE = SMatrix{N,N}(I)
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))
    P = @variable(model, [1:N, 1:N], Symmetric)

    @constraint(model, Symmetric(P - EYE) in PSDCone())
    @constraint(model, tr(P) <= max_trace)

    for (x, y, ηm, ηa) in zip(x_list, y_list, ηm_list, ηa_list)
        aff = y'*P*y - γ2*(x'*P*x + ηa)*(1.0 + ηm)
        @constraint(model, aff <= 0.0)
    end

    @objective(model, Min, sum(P[:].^2))

    return model, P
end

function _optim_model_singular(::Val{N},
        x_list, y_list, ηm_list, ηa_list, γ) where {N}
    γ2 = γ^2
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))
    P = @variable(model, [1:N, 1:N], PSD)

    @constraint(model, tr(P) >= 1)

    for (x, y, ηm, ηa) in zip(x_list, y_list, ηm_list, ηa_list)
        aff = y'*P*y - γ2*(x'*P*x + ηa)*(1.0 + ηm)
        @constraint(model, aff <= 0.0)
    end

    @objective(model, Min, sum(P[:].^2))

    return model, P
end