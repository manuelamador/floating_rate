# # An example with a permanent interest rate increase

import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using LTBonds 
using KernelDensity
using Random
using PrettyTables
using Plots

# The following are the parameters of Chatterjee and Eyigungor (2012)
# But we start with a one period bond
mLB = let
    R = 1.01   # annual interest rate of ~4% 
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 200, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.003, span = 2.0)
    bond = BondCE2012(n = 350, min = 0.0, max = 1.5, κ = 0.03, λ = 0.05)
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

@time solve!(mLB; max_iters = 1000, g = 1.0, err = 1e-10, print_every = 50)


mLB_2 = let
    R = 1.03  # permanent interest rate increase to ~12.6%
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 200, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.003, span = 2.0)
    bond = BondCE2012(n = 350, min = 0.0, max = 1.5, κ = 0.03, λ = 0.05)
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

@time solve!(mLB_2; max_iters = 1000, g = 1.0, err = 1e-10, print_every = 50)

# ### Simulations and moments

big_T = 20_000 
big_N = 1_000
rng = Random.seed!(1234)


shocks, paths = create_shocks_paths(mLB, big_T, big_N; rng);
simulation!(paths, shocks, mLB; n = big_T, trim = 1000, trim_def = 20)
moments_LB = moments(paths, mLB);
pretty_table(pairs(moments_LB))
dens = kde((vcat([p.y for p in paths]...), vcat([p.b for p in paths]...)))


shocks_2, paths_2 = create_shocks_paths(mLB_2, big_T, big_N; rng);
simulation!(paths_2, shocks_2, mLB_2; n = big_T, trim = 1000, trim_def = 20)
moments_LB_2 = moments(paths_2, mLB_2);
pretty_table(pairs(moments_LB_2))
dens_2 = kde((vcat([p.y for p in paths_2]...), vcat([p.b for p in paths_2]...)))


# #### Plots 

f = plot()
contour!(f, get_y_grid(mLB), get_b_grid(mLB), get_d_pol(mLB), xlabel = "y", ylabel = "b", c = palette([:black, :white]),  levels = [-0.06, 0.0, 0.06])
contour!(f, get_y_grid(mLB_2), get_b_grid(mLB_2), get_d_pol(mLB_2), color = palette([:blue, :red]), levels = [-0.06, 0.0, 0.06])
contour!(f, collect(dens.x), collect(dens.y), collect(dens.density'), colorbar=false, levels = 25)
contour!(f, collect(dens_2.x), collect(dens_2.y), collect(dens_2.density'), colorbar=false, levels = 10, color = :blue)


f2 = let b1 = get_b_grid(mLB), 
        q1 = get_q(mLB) ./ LTBonds.risk_free_price(mLB.model), 
        b2 = get_b_grid(mLB_2)
        q2 = get_q(mLB_2) ./ LTBonds.risk_free_price(mLB_2.model)

    plot(b1, q1[:, 1], lw = 3, c = :black, title = "q(b)/qRF", xlabel = "b") 
    plot!(b1, q1[:, 100], lw = 3, c = :black) 
    plot!(b1, q1[:, end], lw = 3, c = :black) 

    plot!(b2, q2[:, 1], c = :blue, lw = 3) 
    plot!(b2, q2[:, 100], c = :blue, lw = 3) 
    plot!(b2, q2[:, end], c = :blue, lw = 3) 
end


f3 = let b1 = get_b_grid(mLB), 
    q1 = get_q(mLB), 
    b2 = get_b_grid(mLB_2)
    q2 = get_q(mLB_2)

plot(b1, q1[:, 1], lw = 3, c = :black, title = "q(b)", xlabel = "b") 
plot!(b1, q1[:, 100], lw = 3, c = :black) 
plot!(b1, q1[:, end], lw = 3, c = :black) 

plot!(b2, q2[:, 1], c = :blue, lw = 3) 
plot!(b2, q2[:, 100], c = :blue, lw = 3) 
plot!(b2, q2[:, end], c = :blue, lw = 3) 
end