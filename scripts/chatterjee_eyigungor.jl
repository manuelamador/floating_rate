# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Julia 20 threads 1.6.3
#     language: julia
#     name: julia-20-threads-1.6
# ---

# # Chatterjee Eyigungor AER 2012
#
# Replicates the Chatterjee and Eyigungor, American Economic Review 2012 paper.
#
# The link to the paper is here: 
#
# https://www.aeaweb.org/articles?id=10.1257/aer.102.6.2674

# +
using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
using LTBonds 
using Plots
using PrettyTables 
using BenchmarkTools
using Random

default(labels = "")
# -

# # Long bond model

# The following are exactly the parameters of Chatterjee and Eyigungor (2012)
# But we start with a one period bond
mLB = let
    R = 1.01
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

# ### Some policy plots 

let a = mLB
    q = get_q(a)
    plot(q[:, 1], marker = :auto, title = "Price functions", xlabel = "b") 
    plot!(q[:, end], marker = 1) 
    plot!(q[:, 100], marker = 2) 
end

let yi = 100, ylow = 1, a = mLB
    yhigh = length(get_y_grid(a))
    b_pol = get_b_pol(a)
    b_grid = get_b_grid(a)
    out_low = [ 
        all_default_at(a, i, ylow) ? missing : b_grid[b_pol[i, ylow][end].idx]  
        for i in eachindex(b_grid)
    ]
    out_high = [ 
        all_default_at(a, i, yhigh) ? missing : b_grid[b_pol[i, yhigh][1].idx]  
        for i in eachindex(b_grid)
    ]
    out_mid = [ 
        all_default_at(a, i, yi) ? missing : b_grid[b_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]
    plot(b_grid, out_mid, lw = 3)
    plot!(b_grid, out_low, lw = 3)
    plot!(b_grid, out_high, lw = 3)

    plot!([0, 1.0], [0, 1.0], ls = :dash)
    plot!(title = "Policy functions")
end

# ### Simulations and moments

big_T = 20_000 
big_N = 1_000
rng = Random.seed!(1234)


@time shocks, paths = create_shocks_paths(mLB, big_T, big_N; rng);

@time simulation!(paths, shocks, mLB; n = big_T, trim = 1000, trim_def = 20)
@time moments_LB = moments(paths, mLB);

pretty_table(pairs(moments_LB))


