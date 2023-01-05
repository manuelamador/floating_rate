# -*- coding: utf-8 -*-
_LOAD_MODELS = false
_SAVE_MODELS = true;

import Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()

using Revise
using LTBonds
using Random 
using PrettyTables
using Plots
using Serialization

benchmark_parameters =  let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 50, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    η=0.1

    (R=R, pref=pref,  y=y, m=m, penalty=penalty, η=η )
end;


if _LOAD_MODELS
    models = deserialize(joinpath(@__DIR__, "..", "output", "tmpcal1.dat"))
else 
    models = let
        R, pref, y, m, penalty, η = benchmark_parameters
        N = 350

        bondLT = Bond(n = N, min = 0.0, max = 1.5, κ = R - 1, λ = 0.05)  
        bondST = Bond(n = N, min = 0.0, max = 1.5, κ = R - 1, λ = 1.0)  
        bondperp = Bond(n = N, min = 0.0, max = 1.5, κ = R - 1, λ = 0.0) 
        
        bondFR = FloatingRateBond(;n = N, min = 0.0, max = 1.5, λ = 0.05, κbar = 0.05) 
        bondFRlowκ = FloatingRateBond(;n = N, min = 0.0, max = 1.5, λ = 0.05, κbar = 0.015) 
        bondFRperp = FloatingRateBond(; n = N, min = 0.0, max = 1.5, λ = 0.0, κbar = 0.05) 
        bondFRST = FloatingRateBond(n = N, min = 0.0, max = 1.5, λ = 1.0, κbar = 0.05)
   
        # EG models

        egperp = LTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondperp, 
            def_costs = penalty, 
            R = R,
        )

        egST = LTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondST, 
            def_costs = penalty, 
            R = R,
        )

        egLT = LTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondLT, 
            def_costs = penalty, 
            R = R,
        )

        frEG = LTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondFR, 
            def_costs = penalty, 
            R = R
        )

        frEGST = LTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondFRST, 
            def_costs = penalty, 
            R = R
        )   

        frEGperp = LTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondFRperp, 
            def_costs = penalty, 
            R = R
        )

        # CK models 

        ckST = CKLTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondST, 
            def_costs = penalty, 
            R = R,
            η = η
        )

        ckLT = CKLTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondLT, 
            def_costs = penalty, 
            R = R,
            η = η
        )

        fr = CKLTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondFR, 
            def_costs = penalty, 
            R = R,
            η = η
        )   

        frlowκ = CKLTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondFRlowκ, 
            def_costs = penalty, 
            R = R,
            η = η
        )   

        frperp = CKLTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bondFRperp, 
            def_costs = penalty, 
            R = R,
            η = η
        )   

        models = (; egLT, egST,  egperp, ckLT, ckST, fr, frlowκ, frperp, frEG, frEGperp, frEGST)
        map(generate_workspace, models)
    end
    
    @time foreach(models) do m 
        @time solve!(m, print_every = 100, max_iters = 5_000)
    end 
end;

_SAVE_MODELS && serialize(joinpath(@__DIR__, "..", "output", "tmpcal1.dat"), models);

big_T = 20_000 
big_N = 1_000
rng = Random.seed!(1234)

shocks, paths = create_shocks_paths(models.fr, big_T, big_N; rng);  # use the same draws. Make sure we use a CK model.

computed_moments = map(models) do (m) 
    simulation!(paths, shocks, m; n = big_T, trim = 1000, trim_def = 20)
    moments(paths, m)
end;

pretty_table(
    collect(map(m -> pairs(m), computed_moments)),
    row_names = collect(keys(computed_moments)), 
    backend = Val(:html)
)

# # Plots 

#set linewidth for plots:
lw = 2
ms = 3
msdiamond = 5
default(size = (600, 400), xtickfontsize = 12, ytickfontsize = 12, yguidefontsize = 14, xguidefontsize = 14)

###Value at zero debt plots###
#Global variables rock!
ygrid = get_y_grid(first(models))
nY = length(ygrid)
midY = nY ÷ 2
β = first(models).model.preferences.β;


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


markerevery(series; n = 1) = collect(view(series, 1:n:length(series)))

let 
    frontiers = Dict( [key => (m, get_frontier(m, midY)...) for (key, m) in pairs(models)])
    
    f = plot()
    for (m, marker, color, s) in zip(
            [:fr, :frEG, :egST, :frEGperp, :frEGST, :frperp], 
            [:circle, :diamond, :square, :ltriangle, :utriangle, :star], 
            [:black, :blue, :red, :gray, :purple, :yellow],
            [1, 2, 3, 4, 5, 6])
        vD = get_vD(models[m])[midY]
        xx = [x  for x in zip(frontiers[m][2], frontiers[m][3]) if x[1] >= vD && x[2] > 0]

        plot!(f, [x[1] for x in xx], [x[2] for x in xx], 
            line = (lw, color), xlabel = ("V"), ylabel = ("MV"), 
            label = nothing)
        plot!(f, [vD], [0], marker = marker, markercolor = color, markersize = 3/2 * (7 - s) + 1, label = nothing)
        plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash, markersize = 3/2 * (7 - s) + 1, 
            line = (lw, color), label = nothing)
        xx_f = markerevery(xx; n = 10)
        scatter!(f, [x[1] for x in xx_f], [x[2] for x in xx_f], 
            markershape = marker, markercolor = color, markersize = 3/2 * (7 - s) + 1,
            label = String(m))
    end 

    f
end 


let f = plot()
    plot!(f, ygrid, get_v(models.egST)[1, :] , line = (lw, :red), markershape = :square)
    plot!(f, ygrid, get_v(models.fr)[1, :, midY], line = (lw, :black), legend = false, xlabel = ("y"), ylabel = ("V(y, b = 0)"), markershape = :circle)
    plot!(f, ygrid, get_v(models.frEG)[1, :, midY] , line = (lw, :blue), markershape = :diamond)
    plot!(f, ygrid, get_v(models.frEGperp)[1, :, midY] , line = (lw, :gray), markershape = :ltriangle)
    f
end 

# # Multiplicity ?

# +
frEG = let 
    R, pref, y, m, penalty, η = benchmark_parameters
    N = 350

    bondFR = FloatingRateBond(;n = N, min = 0.0, max = 1.5, λ = 0.05, κbar = 1.0) 
    
    extra_model = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bondFR, 
        def_costs = penalty, 
        R = R
    )

    ws = generate_workspace(extra_model)
    
    @time solve!(ws, print_every = 100, max_iters = 5_000)
    ws
end;

simulation!(paths, shocks, frEG; n = big_T, trim = 1000, trim_def = 20)
moments(paths, frEG)

# +
frEG_2 = let 
    R, pref, y, m, penalty, η = benchmark_parameters
    N = 350

    bondFR = FloatingRateBond(;n = N, min = 0.0, max = 1.5, λ = 0.05, κbar = 1.0) 
    
    extra_model = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bondFR, 
        def_costs = penalty, 
        R = R
    )

    ws = copy_workspace(extra_model, models.frEGST)
    
    @time solve!(ws, print_every = 100, max_iters = 5_000)
    ws
end;

simulation!(paths, shocks, frEG_2; n = big_T, trim = 1000, trim_def = 20)
moments(paths, frEG_2)
# -

let 
    models = (; frEG, frEG_2)
    
    frontiers = Dict( [key => (m, get_frontier(m, midY)...) for (key, m) in pairs(models)])
    
    f = plot()
    for (m, marker, color) in zip(
            keys(models), 
            [:circle, :diamond, :square, :ltriangle, :utriangle, :star], 
            [:black, :blue, :red, :gray, :purple, :yellow])
        vD = get_vD(models[m])[midY]
        xx = [x  for x in zip(frontiers[m][2], frontiers[m][3]) if x[1] >= vD && x[2] > 0]

        plot!(f, [x[1] for x in xx], [x[2] for x in xx], 
            line = (lw, color), xlabel = ("V"), ylabel = ("MV"), 
            label = nothing)
        plot!(f, [vD], [0], marker = marker, markercolor = color, label = nothing)
        plot!(f, [vD, xx[end][1]], [0, xx[end][2]], ls = :dash,line = (lw, color), label = nothing)
        xx_f = markerevery(xx; n = 10)
        scatter!(f, [x[1] for x in xx_f], [x[2] for x in xx_f], 
            markershape = marker, markercolor = color, label = String(m))
    end 

    f
end 


