# Area of the two-sided spherical cap, i.e., going from abscissa -1 to s and
# from s to 1; assuming the total area (i.e., when s=0) is 1.
function area_2sided_cap(s::Real, n::Integer)
    Dist = Beta((n - 1)/2, 1/2)
    return cdf(Dist, 1 - s^2)
end

# Abscissa of the two-sided spherical cap, i.e., going from abscissa -1 to s and
# from s to 1 such that the area is A; assuming the total area is 1.
function area_2sided_cap_inv(A::Real, n::Integer)
    A = max(0, min(1, A))
    Dist = Beta((n - 1)/2, 1/2)
    return sqrt(1 - quantile(Dist, A))
end