

abstract type AbstractModel{B} end


struct LTBondModel{B<:AbstractBond, T0<:Preferences, T1<:YDiscretized, T2<:MShock, T3<:DefCosts, T4<:Real} <: AbstractModel{B}
    preferences::T0
    y::T1
    m::T2   
    bond::B
    def_costs::T3
    R::T4
end

LTBondModel(;y, m, preferences, bond, def_costs, R) = LTBondModel(preferences, y, m, bond, def_costs, R)
get_base_pars(model::LTBondModel) = model 


struct CKLTBondModel{B<:AbstractBond, T0<:Preferences, T1<:YDiscretized, T2<:MShock, T3<:DefCosts, T4<:Real, T5<:Real} <: AbstractModel{B}
    preferences::T0
    y::T1
    m::T2   
    bond::B
    def_costs::T3
    R::T4
    η::T5
end

CKLTBondModel(;y, m, preferences, bond, def_costs, R, η) = CKLTBondModel(preferences, y, m, bond, def_costs, R, η)
get_base_pars(model::CKLTBondModel) = model


get_bond(model) = get_base_pars(model).bond
get_m(model) = get_base_pars(model).m
get_preferences(model) = get_base_pars(model).preferences
get_R(model) = get_base_pars(model).R
get_y(model) = get_base_pars(model).y
get_η(model) = get_base_pars(model).η
get_b_grid(model) = get_bond(model).grid
get_y_grid(model) = get_y(model).grid
get_u(model) = get_preferences(model).u


mutable struct WorkSpace{M0, M1, M2, M3, M4}
    model::M0
    new::M1
    current::M2
    policies::M3
    cache::M4
end 

get_preferences(ws::WorkSpace) = get_preferences(ws.model)
get_u(ws::WorkSpace) = get_u(ws.model)
get_b_grid(ws::WorkSpace) = get_b_grid(ws.model)
get_y_grid(ws::WorkSpace) = get_y_grid(ws.model)
get_bond(ws::WorkSpace) = get_bond(ws.model)

get_new(ws) = ws.new
get_current(ws) = ws.current
get_policies(ws) = ws.policies
get_cache(ws) = ws.cache
get_base_pars(ws) = get_base_pars(ws.model)

get_v(ws) = ws.current.v
get_q(ws) = ws.current.q
get_κ(ws) = ws.current.κ

get_vD(ws) = ws.cache.vD
get_Ev(ws) = ws.cache.Ev

get_d_pol(ws) = ws.policies.d
get_b_pol(ws) = ws.policies.b
get_ck_d_pol(ws) = ws.policies.ck_d
get_bond_return(ws) = ws.cache.bond_return 
get_decay(ws) = ws.cache.decay_b
get_repay(ws) = ws.cache.repay
get_v_new(ws) = ws.new.v

struct Iterables{T0, T1, T2, T3}
    v::T0
    vD1::T1
    q::T2
    κ::T3
end 

Base.copy(i::Iterables) = Iterables(copy(i.v), copy(i.vD1), copy(i.q), copy(i.κ))

struct Policies{T0, T1, T2}
    b::T0
    d::T1
    ck_d::T2
end 

Base.copy(p::Policies) = Policies(copy(p.b), copy(p.d), copy(p.ck_d))


struct Cache{T1, T2, T3, T4, T5, T6, T7, T8}
    vD::T1
    Ev::T2
    cdef::T3
    ucdef::T4
    tmp_vD::T5
    bond_return::T6
    decay_b::T7
    repay::T8
end 

Base.copy(c::Cache) = Cache(copy(c.vD), copy(c.Ev), copy(c.cdef), copy(c.ucdef), copy(c.tmp_vD), copy(c.bond_return), copy(c.decay_b), copy(c.repay)) 

abstract type AbstractRuns end 

struct WithRuns <: AbstractRuns end 
struct NoRuns <: AbstractRuns end 

has_runs(::LTBondModel) = NoRuns()
has_runs(::CKLTBondModel) = WithRuns()
