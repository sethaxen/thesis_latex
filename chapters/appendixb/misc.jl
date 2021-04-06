using StatsFuns, Distributions, QuadGK, Plots

@inline legendre_recur(l, x, Plx, Plm1x) = ((2l + 1) * x * Plx - l * Plm1x) / (l + 1)

function legendre_series(lmax, x)
    T = eltype(x / one(lmax))
    vals = Vector{T}(undef, lmax + 1)
    @inbounds Plm1x = vals[1] = one(T)
    lmax > 0 || return vals
    @inbounds Plx = vals[2] = T(x)
    l = 1
    while l < lmax
        (Plx, Plm1x) = (legendre_recur(l, x, Plx, Plm1x), Plx)
        l += 1
        @inbounds vals[l + 1] = Plx
    end
    return vals
end

using StatsFuns, Distributions, QuadGK, Plots

# struct ModifiedWrappedNormal{T} <: Distributions.ContinuousUnivariateDistribution
#     μ::T
#     σ::T
#     lognorm::T
# end

# function ModifiedWrappedNormal(μ, σ; nmax = 100)
#     (inorm, _) = quadgk(0, π) do x
#         return exp(_modified_wrapped_normal_logpdf_unnorm(x, μ, σ; nmax = nmax)) * sin(x)
#     end
#     return ModifiedWrappedNormal(μ, σ, -log(inorm))
# end

# function _modified_wrapped_normal_logpdf_unnorm(x, μ, σ; nmax = 100)
#     dnorm = Normal(μ, σ)
#     lpold = -Inf
#     lp = logaddexp(logpdf(dnorm, x), logpdf(dnorm, -x))
#     twopi = 2π
#     n = 0
#     while n < nmax
#         n += 1
#         lpold = lp
#         lp = logaddexp(lp, logpdf(dnorm, n * twopi + x))
#         lp = logaddexp(lp, logpdf(dnorm, n * twopi - x))
#     end
#     return lp
# end

# function Distributions.logpdf(d::ModifiedWrappedNormal, x::Real; nmax = 100)
#     lp = _modified_wrapped_normal_logpdf_unnorm(x, d.μ, d.σ; nmax = nmax)
#     return lp + d.lognorm
# end

# function Distributions.pdf(d::ModifiedWrappedNormal, x::Real; nmax = 100)
#     return exp(logpdf(d, x; nmax = nmax))
# end

struct PolarNormal{T} <: Distributions.ContinuousUnivariateDistribution
    μ::T
    σ::T
end

PolarNormal(μ, σ) = PolarNormal{Base.promote_typeof(μ, σ)}(μ, σ)

function Distributions.pdf(d::PolarNormal, x::Real; nmax = 100)
    Pnx = legendre_series(nmax, cos(x))
    Pnμ = legendre_series(nmax, cos(d.μ))
    n = 0:nmax
    gn = (2 .* n .+ 1) .* exp.((-n .* (n .+ 1)) .* d.σ^2 / 2) .* Pnμ
    p = (gn'Pnx) / 2
    return p
end

function Distributions.pdf(d::PolarNormal, xs::AbstractArray; nmax = 100)
    Pnμ = legendre_series(nmax, cos(d.μ))
    n = 0:nmax
    gn = (2 .* n .+ 1) .* exp.((-n .* (n .+ 1)) .* d.σ^2 / 2) .* Pnμ
    p = [(gn'legendre_series(nmax, cos(x))) / 2 for x in xs]
    return p
end

function Distributions.logpdf(d::PolarNormal, x::Real; nmax = 100)
    return log(abs(pdf(d, x; nmax = nmax)))
end

function Distributions.logpdf(d::PolarNormal, xs::AbstractArray; nmax = 100)
    return log.(abs.(pdf(d, x; nmax = nmax)))
end

struct PolarApprox{T} <: Distributions.ContinuousUnivariateDistribution
    a::Vector{T}
end

function Distributions.pdf(d::PolarApprox, x::Real)
    nmax = length(d.a) - 1
    n = 0:nmax
    Pnx = legendre_series(nmax, cos(x))
    c = (2 .* n .+ 1) .* Pnx
    return d.a'c / 2 / d.a[1] # ensure normalized
end

function Distributions.logpdf(d::PolarApprox, x::Real)
    return log(abs(pdf(d, x)))
end

function PolarApprox(d::PolarNormal, nmax = 100)
    Pnμ = legendre_series(nmax, cos(d.μ))
    n = 0:nmax
    an = exp.((-n .* (n .+ 1)) .* d.σ^2 / 2) .* Pnμ
    return PolarApprox(an)
end

function expect(d, f)
    return quadgk(0, π) do θ
        f(θ) * pdf(d, θ) * sin(θ)
    end[1]
end
