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

# # Eaton-Gersovitz and Cole-Kehoe model
#
# Compares the model without self-fulfilling shocks with the same version of the model with them. It does this for several specifications.

# + tags=[]
import Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.resolve(); Pkg.instantiate()
# -

using Revise 
using LTBonds
using Plots 
using BenchmarkTools

default(label = "")

# ## One period bond

mEG, mCK = let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 200, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    
    # one period period bond:
    bond = Bond(n = 350, min = 0.0, max = 1.5, κ = R - 1, λ = 1.0)  # one period  

    eg = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R
    )

    ck = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R,
        η = 0.05
    )
        
    (generate_workspace(eg),  generate_workspace(ck))
end;

# + tags=[]
@time solve!(mEG; max_iters = 2000, g = 1.0, err = 1e-10, print_every = 50)
@time solve!(mCK; max_iters = 2000, g = 1.0, err = 1e-10, print_every = 50)

# +
let yi = 100, ylow = 1
    yhigh = length(get_y_grid(mEG))
    b_grid = get_b_grid(mEG)
    mCKb_pol = get_b_pol(mCK)
    mEGb_pol = get_b_pol(mEG)
    
    out_low = [ 
        all_default_at(mEG, i, ylow) ? missing : b_grid[mEGb_pol[i, ylow][end].idx]  
        for i in eachindex(b_grid)
    ]
    out_high = [ 
        all_default_at(mEG, i, yhigh) ? missing : b_grid[mEGb_pol[i, yhigh][1].idx]  
        for i in eachindex(b_grid)
    ]
    out_mid = [ 
        all_default_at(mEG, i, yi) ? missing : b_grid[mEGb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]
    
    out_low_CK = [ 
        all_default_at(mCK, i, ylow) ? missing : b_grid[mCKb_pol[i, ylow][end].idx]  
        for i in eachindex(b_grid)
    ]
    out_high_CK = [ 
        all_default_at(mCK, i, yhigh) ? missing : b_grid[mCKb_pol[i, yhigh][1].idx]  
        for i in eachindex(b_grid)
    ]
    out_mid_CK = [ 
        all_default_at(mCK, i, yi) ? missing : b_grid[mCKb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]

    plot(b_grid, out_mid, lw = 2)
    plot!(b_grid, out_low, lw = 2)
    plot!(b_grid, out_high, lw = 2)

    plot!([0, 1.5], [0, 1.5], lw = 1, ls = :dot)
    plot!(title = "Policy functions: EG (solid) CK (dashed)")

    plot!(b_grid, out_mid_CK, lw = 2, ls = :dash)
    plot!(b_grid, out_low_CK, lw = 2, ls = :dash)
    plot!(b_grid, out_high_CK, lw = 2, ls = :dash)

end

# -

let 
    b_grid = get_b_grid(mEG)
    qEG = get_q(mEG)
    qCK = get_q(mCK)
    
    plot(b_grid, qEG[:, 100], lw = 3) 
    plot!(b_grid, qEG[:, 1], lw = 3, title = "Prices: EG (solid) CK (dashed)", xlabel = "b") 
    plot!(b_grid, qEG[:, end], lw = 3) 
    
    plot!(b_grid, qCK[:, 100], lw = 2, ls = :dash) 
    plot!(b_grid, qCK[:, 1], lw =  2, ls = :dash) 
    plot!(b_grid, qCK[:, end], lw = 2, ls = :dash) 
end 


# ## Long Bonds

mEGLB, mCKLB = let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 200, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    
    bond = Bond(n = 350, min = 0.0, max = 1.5, κ = 0.03, λ = 0.2) 
    
    eg = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R
    )

    ck = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R,
        η = 0.05
    )
        
    (generate_workspace(eg),  generate_workspace(ck))
end;

# + tags=[]
@time solve!(mEGLB; max_iters = 2000, g = 1.0, err = 1e-7, print_every = 50)
@time solve!(mCKLB; max_iters = 2000, g = 1.0, err = 1e-7, print_every = 50)

# + tags=[]
let yi = 100, ylow = 1, mEG = mEGLB, mCK = mCKLB, 
    yhigh = length(get_y_grid(mEG))
    b_grid = get_b_grid(mEG)
    mCKb_pol = get_b_pol(mCK)
    mEGb_pol = get_b_pol(mEG)
    
    out_low = [ 
        all_default_at(mEG, i, ylow) ? missing : b_grid[mEGb_pol[i, ylow][end].idx]  
        for i in eachindex(b_grid)
    ]
    out_high = [ 
        all_default_at(mEG, i, yhigh) ? missing : b_grid[mEGb_pol[i, yhigh][1].idx]  
        for i in eachindex(b_grid)
    ]
    out_mid = [ 
        all_default_at(mEG, i, yi) ? missing : b_grid[mEGb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]
    
    out_low_CK = [ 
        all_default_at(mCK, i, ylow) ? missing : b_grid[mCKb_pol[i, ylow][end].idx]  
        for i in eachindex(b_grid)
    ]
    out_high_CK = [ 
        all_default_at(mCK, i, yhigh) ? missing : b_grid[mCKb_pol[i, yhigh][1].idx]  
        for i in eachindex(b_grid)
    ]
    out_mid_CK = [ 
        all_default_at(mCK, i, yi) ? missing : b_grid[mCKb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]

    plot(b_grid, out_mid, lw = 2)
    plot!(b_grid, out_low, lw = 2)
    plot!(b_grid, out_high, lw = 2)

    plot!([0, 1.5], [0, 1.5], lw = 1, ls = :dot)
    plot!(title = "Policy functions: EG (solid) CK (dashed)")

    plot!(b_grid, out_mid_CK, lw = 2, ls = :dash)
    plot!(b_grid, out_low_CK, lw = 2, ls = :dash)
    plot!(b_grid, out_high_CK, lw = 2, ls = :dash)

end

# -

let mEG = mEGLB, mCK = mCKLB
    b_grid = get_b_grid(mEG)
    qEG = get_q(mEG)
    qCK = get_q(mCK)
    
    plot(b_grid, qEG[:, 100], lw = 3) 
    plot!(b_grid, qEG[:, 1], lw = 3, title = "Prices: EG (solid) CK (dashed)", xlabel = "b") 
    plot!(b_grid, qEG[:, end], lw = 3) 
    
    plot!(b_grid, qCK[:, 100], lw = 2, ls = :dash) 
    plot!(b_grid, qCK[:, 1], lw =  2, ls = :dash) 
    plot!(b_grid, qCK[:, end], lw = 2, ls = :dash) 
end 


# ## No persistent income shocks

# ### One period

detEG, detCK = let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 1, ρ = 0.948503, std = 0.0, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.0001, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    
    bond = Bond(n = 350, min = 0.0, max = 3.0, κ = 0.03, λ = 1.0)

    eg = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R
    )

    ck = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R, 
        η = 0.05
    )

    generate_workspace(eg), generate_workspace(ck)
end;

@time solve!(detEG; max_iters = 2000, g = 1.0, err = 1e-7, print_every = 50)
@time solve!(detCK; max_iters = 2000, g = 1.0, err = 1e-7, print_every = 50)


# +
let yi = 1, ylow = 1, mCK = detCK, mEG = detEG
    b_grid = get_b_grid(mCK) 
    mCKb_pol = get_b_pol(mCK)
    mEGb_pol = get_b_pol(mEG)

    outEG_high = [ 
        all_default_at(mEG, i, yi) ? missing : b_grid[mEGb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]
    outEG_low = [ 
        all_default_at(mEG, i, yi) ? missing : b_grid[mEGb_pol[i, yi][1].idx]  
        for i in eachindex(b_grid)
    ]


    outCK_high = [ 
        all_default_at(mCK, i, yi) ? missing : b_grid[mCKb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]

    outCK_low = [ 
        all_default_at(mCK, i, yi) ? missing : b_grid[mCKb_pol[i, yi][1].idx]  
        for i in eachindex(b_grid)
    ]

    plot(b_grid, outEG_high, lw = 2)
    plot!(b_grid, outEG_low, lw = 2)

    plot!(b_grid, outCK_high, lw = 2, ls = :dash)
    plot!(b_grid, outCK_low, lw = 2, ls = :dash)

    plot!([0, last(b_grid)], [0, last(b_grid)], lw = 1, ls = :dot)
    plot!(title = "Policy function: EG vs CK (only idd shocks)\nOne period bond")

end

# -

let  mCK = detCK, mEG = detEG
    b_grid = get_b_grid(mEG)
    qEG = get_q(mEG)
    qCK = get_q(mCK)
    
    plot(b_grid, qEG[:, 1], lw = 3, title = "Prices: EG (solid) CK (dashed)", xlabel = "b") 
    plot!(b_grid, qEG[:, end], lw = 3) 
    
    plot!(b_grid, qCK[:, 1], lw =  2, ls = :dash) 
    plot!(b_grid, qCK[:, end], lw = 2, ls = :dash) 
end 

# ### Long Bond

detEGLB, detCKLB = let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 1, ρ = 0.948503, std = 0.0, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 10.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    
    bond = Bond(n = 350, min = 0.0, max = 2.0, κ = 0.03, λ = 0.05)

    eg = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R
    )

    ck = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bond, 
        def_costs = penalty, 
        R = R, 
        η = 0.05
    )

    generate_workspace(eg), generate_workspace(ck)
end;

@time solve!(detEGLB; max_iters = 2000, g = 1.0, err = 1e-7, print_every = 50)
@time solve!(detCKLB; max_iters = 10_000, g = 0.1, err = 1e-7, print_every = 50)


# +
let yi = 1, ylow = 1, mCK = detCKLB, mEG = detEGLB
    b_grid = get_b_grid(mCK) 
    mCKb_pol = get_b_pol(mCK)
    mEGb_pol = get_b_pol(mEG)
    
    outEG_high = [ 
        all_default_at(mEG, i, yi) ? missing : b_grid[mEGb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]
    outEG_low = [ 
        all_default_at(mEG, i, yi) ? missing : b_grid[mEGb_pol[i, yi][1].idx]  
        for i in eachindex(b_grid)
    ]


    outCK_high = [ 
        all_default_at(mCK, i, yi) ? missing : b_grid[mCKb_pol[i, yi][end].idx]  
        for i in eachindex(b_grid)
    ]

    outCK_low = [ 
        all_default_at(mCK, i, yi) ? missing : b_grid[mCKb_pol[i, yi][1].idx]  
        for i in eachindex(b_grid)
    ]

    plot(b_grid, outEG_high, lw = 2)
    plot!(b_grid, outEG_low, lw = 2)

    plot!(b_grid, outCK_high, lw = 2, ls = :dash)
    plot!(b_grid, outCK_low, lw = 2, ls = :dash)

    plot!([0, last(b_grid)], [0, last(b_grid)], lw = 1, ls = :dot)
    plot!(title = "Policy function: EG vs CK (only iid shocks)\nLong Bonds")

end

# -

let  mCK = detCKLB, mEG = detEGLB
    b_grid = get_b_grid(mEG)
    qEG = get_q(mEG)
    qCK = get_q(mCK)
    
    plot(b_grid, qEG[:, 1], lw = 3, title = "Prices: EG (solid) CK (dashed)", xlabel = "b") 
    plot!(b_grid, qEG[:, end], lw = 3) 
    
    plot!(b_grid, qCK[:, 1], lw =  2, ls = :dash) 
    plot!(b_grid, qCK[:, end], lw = 2, ls = :dash) 
end 


