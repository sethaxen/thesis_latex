using ChainRules, LaTeXStrings, PGFPlotsX

diffquotmat(f, λi, Δλ) = (f(λi + Δλ) - f(λi)) / Δλ

function diffquotmatapprox(f, λi, Δλ)
    df(x) = last(frule((Zero(), One()), f, x))
    return (df(λi) + df(λi + Δλ)) / 2
end

function compute_errors(f, λi, n = 1000)
    Δλ = exp.(range(log(eps()), 0; length = n))
    dq = diffquotmat.(f, λi, Δλ)
    dqapprox = diffquotmatapprox.(f, λi, Δλ)
    abserr = abs.(dq .- dqapprox)
    totalerr_exp = abs.(Δλ .^ 2 + eps(typeof(λi)) ./ Δλ)
    return Δλ, (abserr, totalerr_exp)
end

function make_plot(f, λi, n = 1000)
    Δλ, (err, err_pred) = compute_errors(f, λi, n)
    return @pgf LogLogAxis(
        {
            xlabel = L"\Delta\lambda",
            ylabel = "Absolute Error",
            xmode = "log",
            ymode = "log",
            height = "8cm",
            width = "8cm",
        },
        PlotInc({no_marks, smooth, color = "black"}, Coordinates(Δλ, err_pred)),
        LegendEntry("predicted"),
        PlotInc({no_marks, smooth, color = "blue!50!white"}, Coordinates(Δλ, err)),
        LegendEntry("actual"),
    )
end

pgfsave("diffquot_error_exp.tikz", make_plot(exp, 0.0))
pgfsave("diffquot_error_log.tikz", make_plot(log, 1.0))
