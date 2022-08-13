# -*- coding: utf-8 -*-
# # Simulations and moments
#
# This script computes the simulations and associated moments of several version of the model, with and without self-fulfilling shocks, with short and long debt, and with the floating rate bond.

# +
import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using LTBonds
using Random 
using PrettyTables
# -


SAVE_MOMENTS = false # set to true to save the moments to file. 

benchmark_parameters =  let
    R = 1.01
    β = 0.9540232420
    pref = Preferences(β = β, u = make_CRRA(ra = 2))
    y = discretize(YProcess(n = 50, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
    m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
    penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
    η=0.1

    (R=R, pref=pref,  y=y, m=m, penalty=penalty, η=η )
end

mFR, mFRlowκ, mEGLT, mCKLT, mEGST, mCKST = let
    R, pref, y, m, penalty, η = benchmark_parameters

    bondFR = FloatingRateBond(;n = 350, min = 0.0, max = 1.5, λ = 0.05, κbar = 1.0) 
    bondFRlowκ = FloatingRateBond(;n = 350, min = 0.0, max = 1.5, λ = 0.05, κbar = 0.015) 
    bondLT = Bond(n = 350, min = 0.0, max = 1.5, κ = R - 1, λ = 0.05)  
    bondST = Bond(n = 350, min = 0.0, max = 1.5, κ = R - 1, λ = 1.0)  

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

    egLT = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bondLT, 
        def_costs = penalty, 
        R = R,
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

    egST = LTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bondST, 
        def_costs = penalty, 
        R = R,
    )

    ckST = CKLTBondModel(
        y = y,
        m = m, 
        preferences = pref, 
        bond = bondST, 
        def_costs = penalty, 
        R = R,
        η = η
    )

    (generate_workspace(fr), generate_workspace(frlowκ), generate_workspace(egLT), generate_workspace(ckLT), generate_workspace(egST), generate_workspace(ckST) )
end;

for m ∈ (mFR, mFRlowκ, mEGST, mCKST, mEGLT, mCKLT)
    @time solve!(m; print_every = 100, max_iters = 5000)
end

big_T = 20_000 
big_N = 1_000
rng = Random.seed!(1234)

@time shocks, paths = create_shocks_paths(mCKST, big_T, big_N; rng) 

##EGST
@time simulation!(paths, shocks, mEGST; n = big_T, trim = 1000, trim_def = 20)
@time moments_EGST = moments(paths, mEGST)

##CKST
@time simulation!(paths, shocks, mCKST; n = big_T, trim = 1000, trim_def = 20)
@time moments_CKST = moments(paths, mCKST)

##EGLT
@time simulation!(paths, shocks, mEGLT; n = big_T, trim = 1000, trim_def = 20)
@time moments_EGLT = moments(paths, mEGLT)

##CKLT
@time simulation!(paths, shocks, mCKLT; n = big_T, trim = 1000, trim_def = 20)
@time moments_CKLT = moments(paths, mCKLT)

##FR
@time simulation!(paths, shocks, mFR; n = big_T, trim = 1000, trim_def = 20)
@time moments_FR = moments(paths, mFR)

##FRlowκ
@time simulation!(paths, shocks, mFRlowκ; n = big_T, trim = 1000, trim_def = 20)
@time moments_FRlowκ = moments(paths, mFRlowκ)

pretty_table(
    [
        pairs(moments_EGST),
        pairs(moments_CKST),
        pairs(moments_EGLT),
        pairs(moments_CKLT),
        pairs(moments_FR),
        pairs(moments_FRlowκ)
    ],
    row_names = ["EGST", "CKST", "EGLT", "CKLT", "FR", "FRlowκ"]
)

SAVE_MOMENTS && open(joinpath(@__DIR__,"..","output","moments.txt"), "w") do f
    pretty_table(f,
    [
        pairs(moments_EGST),
        pairs(moments_CKST),
        pairs(moments_EGLT),
        pairs(moments_CKLT),
        pairs(moments_FR),
        pairs(moments_FRlowκ)
    ],
    row_names = ["EGST", "CKST", "EGLT", "CKLT", "FR", "FRlowκ"]
)
end
