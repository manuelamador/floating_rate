# -*- coding: utf-8 -*-
# # Floating rate bond plots
#
# Plots several specification of the model without runs, with runs, with short short and long bonds and with the associated floating rate bond. 

import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using LTBonds
using Random 
using Plots
using PrettyTables
using LaTeXStrings 

SAVE_FIGS = false; # set to true to save the figures to files. 

# the floating rate model:
mFR = let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 50, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    bond = FloatingRateBond(;n = 350, min = 0.0, max = 1.5, λ = 0.05, κbar = 1.0)  
    fr = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R,
        η = 0.1
    )   
    generate_workspace(fr)
end;

mEGLT, mCKLT = let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 50, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    bond = Bond(n = 350, min = 0.0, max = 1.5, κ = R - 1, λ = 0.05)  

    eg = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R,
    )

    ck = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R,
        η = 0.1
    )

        
    (generate_workspace(eg), generate_workspace(ck))
end;


mEGST, mCKST = let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 50, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    bond = Bond(n = 350, min = 0.0, max = 1.5, κ = R - 1, λ = 1.0)  

    eg = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R,
    )

    ck = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R,
        η = 0.1
    )

        
    (generate_workspace(eg), generate_workspace(ck))
end;

for m ∈ (mFR, mEGST, mCKST, mEGLT, mCKLT)
    @time solve!(m; print_every = 200, max_iters = 5000)
end


# ## Plots 

#set linewidth for plots:
lw = 2
ms = 3
msdiamond = 5
default(size = (600, 400), xtickfontsize = 12, ytickfontsize = 12, yguidefontsize = 14, xguidefontsize = 14)

###Value at zero debt plots###
#Global variables rock!
ygrid = get_y_grid(mFR)
nY = length(ygrid)
midY = nY ÷ 2
β = mFR.model.preferences.β;


f = plot(ygrid, get_v(mFR)[1,:,midY], line = (lw, :black), legend = false, xlabel = (L"$y$"),
        ylabel = (L"$V(y,b=0)$"))
plot!(f, ygrid, get_v(mCKST)[1,:] , line = (lw, :red), markershape = :circle)
plot!(f, ygrid, get_v(mEGST)[1,:] , line = (lw, :blue), markershape = :diamond)
SAVE_FIGS && savefig(f, joinpath(@__DIR__, "..", "output", "V0ST.pdf"))
f

f = plot(ygrid, get_v(mFR)[1,:,midY], line = (lw, :black), legend=false,
        xlabel = (L"$y$"), ylabel = (L"$V(y,b=0)$"))
plot!(f, ygrid, get_v(mCKLT)[1,:] , line = (lw, :red), markershape = :circle)
plot!(f, ygrid, get_v(mEGLT)[1,:] , line = (lw, :blue), markershape = :diamond)
SAVE_FIGS && savefig(f, joinpath(@__DIR__, "..", "output", "V0LT.pdf"))
f

# ### Welfare


welfareST = [
    inv_u.(Ref(get_u(mFR)), get_v(mFR)[1, i, 1]) ./ 
    inv_u.(Ref(get_u(mCKST)), get_v(mCKST)[1,i]) for i = 1:length(-get_v(mFR)[1, :, 1])];

f = plot(ygrid, welfareST, line=(lw, :black), legend = false, xlabel = L"$y$",
        ylabel = "Welfare Gain")
SAVE_FIGS && savefig(f, joinpath(@__DIR__, "..", "output", "WST.pdf"))
f

welfareLT = [
    inv_u.(Ref(get_u(mFR)), get_v(mFR)[1, i, 1]) ./ 
    inv_u.(Ref(get_u(mCKLT)), get_v(mCKLT)[1,i]) for i = 1:length(-get_v(mFR)[1, :, 1])];

f = plot(ygrid, welfareLT, line = (lw, :black), legend = false, xlabel = L"$y$",
    ylabel = "Welfare Gain")
SAVE_FIGS && savefig(f, joinpath(@__DIR__, "..", "output", "WLT.pdf"))
f

# ## Pareto Frontiers####

_get_frontier_y_state(::AbstractFixedRateBond, yi, lag_y) = (yi, )
_get_frontier_y_state(::AbstractFloatingRateBond, yi, lag_y) = (yi, lag_y)
_get_frontier_κ(bond::AbstractFixedRateBond, _, _, _) = get_κ(bond)
_get_frontier_κ(::AbstractFloatingRateBond, m, bi, lag_y) = m.current.κ[bi, lag_y]
_get_frontier_mv(bond::AbstractFixedRateBond, q, b, _) = find_bond_return(bond; q) * b
_get_frontier_mv(bond::AbstractFloatingRateBond, q, b, κ) = find_bond_return(bond; q, κ) * b

get_frontier(m, yi; kwargs...) = get_frontier(get_bond(m), m, yi; kwargs...)

function get_frontier(bond, m::LTBonds.WorkSpace, yi; lag_y = midY)
    bgrid = get_b_grid(m)
    λ = get_λ(bond)
    mv = []
    v = []
    y_state = _get_frontier_y_state(bond, yi, lag_y)
    for bi=1:length(bgrid)
        default = m.policies.d[bi, y_state...] == m.model.m.m_min
        biprime = get_b_pol(m)[bi, y_state...][end, 1].idx
        q = get_q(m)[biprime, yi]
        κ = _get_frontier_κ(bond, m, bi, lag_y)
        b = get_b_grid(m)[bi]
        push!(mv, (1 - default) * _get_frontier_mv(bond, q, b, κ))
        push!(v, (1 - default) * get_v(m)[bi, y_state...] + default * get_vD(m)[yi])
    end
    return v, mv
end

frontiers = [(m, get_frontier(m,midY)...) for m in (mFR, mEGST, mCKST, mEGLT, mCKLT)];

markerevery(series; n = 1) = collect(view(series, 1:n:length(series)))

let 
    vD = get_vD(mFR)[midY]
    mi = 1
    xx = [x  for x in zip(frontiers[mi][2], frontiers[mi][3]) if x[1] > vD && x[2] > 0]

    f = plot([x[1] for x in xx], [x[2] for x in xx], 
        line=(lw, :black), legend=false, xlabel=(L"$V$"), ylabel=(L"MV")
    )
    plot!(f, [vD], [0], marker = 4, markercolor=:black)
    plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash,line=(lw, :black))

    vD = get_vD(mEGST)[midY]
    mi = 2
    xx = [x  for x in zip(frontiers[mi][2], frontiers[mi][3]) if x[1] > vD && x[2] > 0]
    plot!(f, [x[1] for x in xx], [x[2] for x in xx], 
        line=(lw, :blue))
    xx_f = markerevery(xx; n = 10)
    scatter!(f, [x[1] for x in xx_f], [x[2] for x in xx_f], 
        markershape=:diamond,markercolor=:blue)
    plot!(f, [vD], [0], markershape=:diamond, marker = 4,markercolor=:blue)
    plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash,line=(lw, :blue))

    vD = get_vD(mCKST)[midY]
    mi = 3
    xx = [x  for x in zip(frontiers[mi][2], frontiers[mi][3]) if x[1] > vD && x[2] > 0]
    plot!(f, [x[1] for x in xx], [x[2] for x in xx], 
        line=(lw, :red))
    xx_f = markerevery(xx; n = 15)
    scatter!(f, [x[1] for x in xx_f], [x[2] for x in xx_f], 
        markershape=:circle,markercolor=:red)
    plot!(f, [vD], [0], markershape=:circle, marker = 4,markercolor=:red)
    plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash,line=(lw, :red))

    SAVE_FIGS && savefig(f, joinpath(@__DIR__,"..","output","FrontierST.pdf"))
    f
end 

let 
    vD = get_vD(mFR)[midY]
    mi = 1
    xx = [x  for x in zip(frontiers[mi][2], frontiers[mi][3]) if x[1] > vD && x[2] > 0]

    f = plot([x[1] for x in xx], [x[2] for x in xx], 
        line=(lw, :black), legend=false, xlabel=(L"$V$"), ylabel=(L"MV")
    )
    plot!(f, [vD], [0], marker = 4, markercolor=:black)
    plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash,line=(lw, :black))

    vD = get_vD(mEGLT)[midY]
    mi = 4
    xx = [x  for x in zip(frontiers[mi][2], frontiers[mi][3]) if x[1] > vD && x[2] > 0]
    plot!(f, [x[1] for x in xx], [x[2] for x in xx], 
        line=(lw, :blue))
    xx_f = markerevery(xx; n = 10)
    scatter!(f, [x[1] for x in xx_f], [x[2] for x in xx_f], 
        markershape=:diamond,markercolor=:blue)
    plot!(f, [vD], [0], markershape=:diamond, marker = 4,markercolor=:blue)
    plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash,line=(lw, :blue))

    vD = get_vD(mCKLT)[midY]
    mi = 5
    xx = [x  for x in zip(frontiers[mi][2], frontiers[mi][3]) if x[1] > vD && x[2] > 0]
    plot!(f, [x[1] for x in xx], [x[2] for x in xx], 
        line=(lw, :red))
    xx_f = markerevery(xx; n = 15)
    scatter!(f, [x[1] for x in xx_f], [x[2] for x in xx_f], 
        markershape=:circle,markercolor=:red)
    plot!(f, [vD], [0], markershape=:circle, marker = 4,markercolor=:red)
    plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash,line=(lw, :red))

    SAVE_FIGS && savefig(f, joinpath(@__DIR__,"..","output","FrontierLT.pdf"))
    f
end 



