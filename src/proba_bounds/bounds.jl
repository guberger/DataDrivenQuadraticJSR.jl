function bisection(f::Function, a::Number, b::Number;
        tol::AbstractFloat=1e-5, maxiter::Integer=100)
    fa = f(a)
    fa*f(b) ≤ 0 || error("Invalid starting interval [a,b]")
    i = 0
    while b - a > tol
        i += 1
        i != maxiter || error("Max iteration exceeded")
        c = (a + b)/2
        fc = f(c)
        if iszero(fc)
            return c, c
        elseif fa*fc > 0
            a = c  # Root is in the right half of [a,b].
            fa = fc
        else
            b = c  # Root is in the left half of [a,b].
        end
    end
    return a, b
end

# 0 ≤ i < N && 0 < ϵ < 1
function monomial_uniform(N, i, ϵ)
    if i == 0
        return (1 - ϵ)^N
    end
    α = max(i, N - i) + 1
    β = N - α + 1
    val = sum(log.(α:N)) + i*log(ϵ) + (N - i)*log(1 - ϵ) - sum(log.(1:β))
    return exp(val)
end

# Bound in Campi & Garatti 2018, Corollary 1
# 0 ≤ ζ ≤ N
# ϵ(N, ζ, β)
function uniform_bound(N::Integer, ζ::Integer, β::Real; tol=1e-5)
    if iszero(β)
        return 1.0, 1.0
    elseif β ≥ 1 || iszero(ζ)
        return 0.0, 0.0
    end

    f(ϵ) = if iszero(ϵ)
        1 - β
    elseif ϵ ≥ 1
        -β
    else
        sum(i -> monomial_uniform(N, i - 1, ϵ), 1:ζ) - β
    end

    return bisection(f, 0.0, 1.0, tol=tol, maxiter=10000)
end

# 0 ≤ k ≤ i ≤ N && 0 ≤ ϵ < 1
function monomial_varying(N, i, k, ϵ)
    val = sum(log.((i - k + 1):(N - k))) -
        (N - i)*log(1 - ϵ) - sum(log.((i + 1):N))
    return exp(val)
end

# Bound in Campi & Garatti 2018, Theorem 4
# ϵ(N, k, β)
function varying_bound_1(N::Integer, k::Integer, β::Real; tol=1e-5)
    if k == N
        return 1.0, 1.0
    elseif iszero(β)
        return 1.0, 1.0
    elseif β ≥ 1
        return 0.0, 0.0
    end

    f(ϵ) = if ϵ ≥ 1
        1.0
    else
        β*sum(i -> monomial_varying(N, i, k, ϵ), k:N)/(N + 1) - 1
    end

    return bisection(f, 0.0, 1.0, tol=tol, maxiter=10000)
end

# Bound in Garatti & Campi 2021, Theorem 1
# ϵ(N, k, β)
function varying_bound_2(N::Integer, k::Integer, β::Real; tol=1e-5)
    if k == N
        return 1.0, 1.0
    elseif iszero(β)
        return 1.0, 1.0
    elseif β ≥ 1
        return 0.0, 0.0
    end

    f(ϵ) = if ϵ ≥ 1
        1.0
    else
        β*sum(i -> monomial_varying(N, i, k, ϵ), k:(N - 1))/N - 1
    end

    return bisection(f, 0.0, 1.0, tol=tol, maxiter=10000)
end

# Bound in Garatti & Campi 2021, Theorem 1
# N(β, ϵ, k)
function varying_bound_2_inv(β::Real, ϵ::Real, k::Integer)
    if β ≥ 1
        return 0
    elseif ϵ ≤ 0 || β ≤ 0
        return -1
    end

    N = k
    F = 0

    while F < N/β
        F += 1
        F *= ((N - k + 1)/(N + 1))/(1 - ϵ)
        N += 1
    end

    return N
end

# Bound β = (1 - ϵ)^N, valid only for the single-variable monotonic case.
# N(β, ϵ, k=1)
function single_monotone_bound_inv(β::Real, ϵ::Real)
    if β ≥ 1 || ϵ ≥ 1
        return 0
    elseif ϵ ≤ 0 || β ≤ 0
        return -1
    end
    return ceil(Int, log(β)/log(1 - ϵ))
end