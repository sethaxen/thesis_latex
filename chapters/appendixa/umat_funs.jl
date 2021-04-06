using BlockArrays, Distributions, LinearAlgebra, Test
using PGFPlotsX, Colors, LaTeXStrings, Random

struct E{D} end
const E1 = E{1}()
const E2 = E{2}()
const E3 = E{3}()
const E4 = E{4}()
const E5 = E{5}()
const E6 = E{6}()

Base.Broadcast.broadcastable(b::E) = Ref(b)

@inline δ(i, j) = i == j
@inline c(l, n) = sqrt((l - n) * (l + n + 1))
@inline function γ(s, l, m)
    return sqrt(
        (l^2 - s^2) * (l - m) * (l - m - 1) /
        ifelse(iszero(l), one(l), l^2 * (4l^2 - 1)),
    )
end
@inline function λ(s, l, m)
    return s * sqrt((l - m) * (l + m + 1)) / ifelse(iszero(l), one(l), l^2 + l)
end
@inline function κ(s, l, m)
    return sqrt(
        (l^2 - m^2) * (l^2 - s^2) / ifelse(iszero(l), one(l), l^2 * (4l^2 - 1)),
    )
end

function u(::typeof(E1), l′, m′, l, m, p, s)
    return im *
           (-(c(l, -m) * δ(l, l′) * δ(m′ + 1, m) + c(l, m) * δ(l, l′) * δ(m′ - 1, m)) / 2)
end

function u(::typeof(E2), l′, m′, l, m, p, s)
    return (c(l, -m) * δ(l, l′) * δ(m′ + 1, m) - c(l, m) * δ(l, l′) * δ(m′ - 1, m)) / 2
end

u(::typeof(E3), l′, m′, l, m, p, s) = im * (-m * δ(l, l′) * δ(m′, m))

function u(::typeof(E4), l′, m′, l, m, p, s)
    return im * (
        (
            -γ(s, l′, -m′) * δ(m′, m + 1) * δ(l′ - 1, l) +
            γ(s, l′, m′) * δ(m′, m - 1) * δ(l′ - 1, l) +
            γ(s, l, m) * δ(m′, m + 1) * δ(l′ + 1, l) -
            γ(s, l, -m) * δ(m′, m - 1) * δ(l′ + 1, l) +
            λ(s, l, m) * δ(m′, m + 1) * δ(l′, l) +
            λ(s, l, -m) * δ(m′, m - 1) * δ(l′, l)
        ) * p / 2
    )
end

function u(::typeof(E5), l′, m′, l, m, p, s)
    return (
        -γ(s, l′, -m′) * δ(m′, m + 1) * δ(l′ - 1, l) -
        γ(s, l′, m′) * δ(m′, m - 1) * δ(l′ - 1, l) +
        γ(s, l, m) * δ(m′, m + 1) * δ(l′ + 1, l) +
        γ(s, l, -m) * δ(m′, m - 1) * δ(l′ + 1, l) +
        λ(s, l, m) * δ(m′, m + 1) * δ(l′, l) - λ(s, l, -m) * δ(m′, m - 1) * δ(l′, l)
    ) * p / 2
end

function u(::typeof(E6), l′, m′, l, m, p, s)
    return im * (
        p *
        δ(m′, m) *
        (
            κ(s, l′, m′) * δ(l′ - 1, l) +
            ((s * m) / ifelse(iszero(l), one(l), l * (l + 1))) * δ(l′, l) +
            κ(s, l, m) * δ(l′ + 1, l)
        )
    )
end

u(b::E, l′, l, p, s) = u.(b, l′, (-l′):l′, l, transpose((-l):l), p, s)
function u(b::E, L, p, s)
    T = typeof(sqrt(s * L) * p)
    S = E ∈ (E2, E5) ? T : complex(T)
    dims = 1:2:(2L + 1)
    umat = BlockArray{S}(undef, dims, dims)
    for l′ = 0:L, l = 0:L
        ublock = u(b, l′, l, p, s)
        setblock!(umat, ublock, l′ + 1, l + 1)
    end
    return umat
end

function Bmatrix(D, L, p, s)
    us = u.((E1, E2, E3, E4, E5, E6), L, p, s)
    T = complex(typeof(zero(eltype(D)) * zero(eltype(first(us)))))
    bmat = fill!(Matrix{T}(undef, size(first(us))), zero(T))
    for i in 1:6, j in 1:6
        bmat .+= D[i, j] .* (us[i] * us[j]) ./ 2
    end
    return bmat
end

@inline blockind(l′, l) = ((l′^2 + 1):((l′ + 1)^2), (l^2 + 1):((l + 1)^2))

function diffusion_mat(; σr = 1, σt = 1, digits = 2)
    d = Diagonal([σr, σr, σr, σt, σt, σt])
    D = Symmetric(round.(d * round.(rand(LKJ(6, 1)); digits=digits) * d; digits=digits))
    return D
end

function diffusion_modify(S, mode = :dense)
    D = Symmetric(copyto!(similar(S), S))
    if mode === :indep
        D.data[1:3, 4:6] .= 0
        return D
    elseif mode === :diag
        return copyto!(D.data, diagm(diag(D)))
    elseif mode === :rot
        D.data[1:6, 4:6] .= 0
        return D
    elseif mode === :trans
        D.data[1:3, 1:6] .= 0
        return D
    elseif mode === :scalarindep
        r = sum(diag(S)[1:3]) / 3
        t = sum(diag(S)[4:6]) / 3
        copyto!(D.data, diagm([r, r, r, t, t, t]))
    else
        return D
    end
end

function checkproperties(L, p, s)
    @testset "$b" for b in (E1, E2, E3, E4, E5, E6)
        mat = u(b, L, p, s)
        @testset "is skew-hermitian" begin
            @test mat ≈ -mat'
        end
        b in (E1, E2, E3) && @testset "is tridiagonal" begin
            @test Tridiagonal(mat) ≈ mat
        end
        b in (E4, E5, E6) && @testset "is block-tridiagonal" begin
            for l′ in 0:L, l in 0:L
                abs(l - l′) > 1 && @test iszero(getblock(mat, l′ + 1, l + 1))
            end
        end
        b in (E4, E5, E6) && @testset "blocks are tridiagonal" begin
            @testset "block ($l′, $l)" for l′ in 0:L, l in max(0, l′ - 1):min(l′ + 1, L)
                matblock = getblock(mat, l′ + 1, l + 1)
                for k in (1 - size(matblock, 1)):(size(matblock, 2) - 1)
                    abs(k) > 2 && @test iszero(diag(matblock, k))
                end
            end
        end
    end
end

function colorize(z)
    r = abs(z)
    h = rad2deg((2π + angle(z)) % 2π)
    l = 1 / (1 + r^0.3)
    s = 0.5
    return HSL(h, s, l)
end

function make_plot(A, L = Int(sqrt(size(A, 1))) - 1)
    l = 0:L
    block_centers = @. l^2 + l + 1
    block_edges = @. block_centers - l - 0.5

    x, y = axes(A)
    xy = Pair.(x, y')
    meta = colorize.(vec(A))
    coords = Coordinates(vec(first.(xy)), vec(last.(xy)), meta = meta)
    @pgf Axis(
        {
            enlargelimits = false,
            "axis equal image",
            "y dir" = "reverse",
            "axis on top" = true,
            xtick = block_centers[2:end],
            xticklabels = l[2:end],
            xlabel = L"$\ell$",
            "axis x line*" = "top",
            ytick = block_centers[2:end],
            yticklabels = l[2:end],
            ylabel = L"$\ell'$",
            "xticklabel style" = {
                font = raw"\tiny"
            },
            "yticklabel style" = {
                font = raw"\tiny"
            },
            tickwidth = 0,
            clip = false,
        },
        PlotInc(
            {
                "matrix plot",
                no_marks,
                "mesh/color input" = "explicit",
                "mesh/cols" = length(x),
            },
            coords,
        ),
        VLine.(Ref({gray, "line width" = "0.1mm"}), block_edges[2:end])...,
        HLine.(Ref({gray, "line width" = "0.1mm"}), block_edges[2:end])...,
    )
end

Random.seed!(57)
D = diffusion_mat(σr = 0.5, σt = 0.5)
L, p, s = 7, 1, 0
ext = "pdf"
for modify in (:dense, :rot, :trans, :indep, :scalarindep)
    B = Bmatrix(diffusion_modify(D, modify), L, p, s)
    pgfsave("Bmat_se3_$modify.$ext", make_plot(B))
    ϕ = exp(B)
    pgfsave("ftrans_se3_$modify.$ext", make_plot(ϕ))
end
