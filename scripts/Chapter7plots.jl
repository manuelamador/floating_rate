# -*- coding: utf-8 -*-
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using LTBonds 
using Plots
using BenchmarkTools
using LaTeXStrings 


modelLB, modelSB, modelLB2 = map((0.05, 1.0, 0.025)) do λ
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 200, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = (λ < 0.05 ? 0.006 : 0.003), span = 2.0)
    bond = BondCE2012(n = 350, min = 0.0, max = 1.5, κ = 0.03, λ = λ)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    generate_workspace(LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R
    ))
end;

for m in (modelLB, modelSB, modelLB2)
    @time solve!(m; max_iters = 1500, g = 1.0, err = 1e-10, print_every = 50)
end 


# ## Figures

#to set spacing of markers and other attributes
@recipe function f(::Type{Val{:samplemarkers}}, x, y, z; step = 10)
    n = length(y)
    sx, sy = x[1:step:n], y[1:step:n]
    # add an empty series with the correct type for legend markers
    @series begin
        seriestype := :path
        markershape --> :auto
        x := []
        y := []
    end
    # add a series for the line
    @series begin
        primary := false # no legend entry
        markershape := :none # ensure no markers
        seriestype := :path
        seriescolor := get(plotattributes, :seriescolor, :auto)
        x := x
        y := y
    end
    # return  a series for the sampled markers
    primary := false
    seriestype := :scatter
    markershape --> :auto
    x := sx
    y := sy
end


function plot_q(idx,a1,a2,a3)
    #set linewidth for plots:
    lw = 2
    ms = 3
    default(size = (600,400), xtickfontsize = 12, ytickfontsize = 12, yguidefontsize = 14, xguidefontsize = 14)
    msdiamond = 5

    qstar1, qstar2, qstar3 = map(risk_free_price ∘ get_base_pars, (a1, a2, a3))
    b_grid1, b_grid2, b_grid3 = map(get_b_grid, (a1, a2, a3))
    q1, q2, q3 = map(get_q, (a1, a2, a3))

    f = plot(b_grid1, q1[:, idx]./qstar1, line = (lw, :black), legend = false, xlabel = (L"$b'$"), ylabel = (L"$q(b',y=y_{min})/q^{RF}$"))
    plot!(f, b_grid2, q2[:, idx]./qstar2, line = (lw,:dash, :black),st = :samplemarkers, step = 20, markercolor=:black, shape = :circle, markersize = ms)
    plot!(f, b_grid3, q3[:, idx]./qstar3, line = (lw,:dashdot, :black),st = :samplemarkers, step = 20, markercolor=:black, shape = :diamond, markersize = msdiamond)

    plot!(f, xlims = (0, b_grid1[end]))
    return f
end


#plot prices at lowest Y (Figure 7-6 (a))
plot_q(1, modelLB, modelSB, modelLB2)

#plot prices at mean Y (Figure 7-6 (b))
midY = length(get_y_grid(modelLB)) ÷ 2
plot_q(midY, modelLB, modelSB, modelLB2)

#plot prices at max Y (Figure 7-6 (c))
plot_q(length(get_y_grid(modelLB)), modelLB, modelSB, modelLB2)



