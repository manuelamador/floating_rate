
abstract type AbstractBond end 
abstract type AbstractFixedRateBond <: AbstractBond end
abstract type AbstractFloatingRateBond <: AbstractBond end


struct Bond{F<:Real, N<:Integer, T1, T2} <: AbstractFixedRateBond
    min::F
    max::F
    n::N
    κ::F  # coupon rate
    λ::F  # inverse maturity
    grid::T1 
    zero::T2 #location of zero debt
end

struct BondCE2012{F<:Real, N<:Integer, T1, T2} <: AbstractFixedRateBond
    min::F
    max::F
    n::N
    κ::F  # coupon rate
    λ::F  # inverse maturity
    grid::T1 
    zero::T2 #location of zero debt
end

struct FloatingRateBond{F<:Real, N<:Integer, T1, T2} <: AbstractFloatingRateBond
    min::F
    max::F
    n::N
    κbar::F  # maximum coupon rate
    λ::F  # inverse maturity
    grid::T1 
    zero::T2 #location of zero debt
end


get_λ(b::AbstractBond) = b.λ
get_grid(b::AbstractBond) = b.grid
get_zero_index(b::AbstractBond) = b.zero

get_κ(b::AbstractFixedRateBond) = b.κ
get_κbar(b::AbstractFloatingRateBond) = b.κbar


function _bond_constructor(;min, max, n, κ, λ)
    grid = collect(LinRange(min, max, n))
    # making sure that 0 is on the grid
    zeroidx = findfirst(x -> x >= 0, grid)
    grid[zeroidx] = zero(eltype(grid))
    return (min, max, n, κ, λ, grid, zeroidx)
end 

Bond(;min, max, n, κ, λ) = Bond(_bond_constructor(;min, max, n, κ, λ)...)
BondCE2012(;min, max, n, κ, λ) = BondCE2012(_bond_constructor(;min, max, n, κ, λ)...)
FloatingRateBond(;min, max, n, κbar, λ) = FloatingRateBond(_bond_constructor(;min, max, n, κ = κbar, λ)...)


# Maturing coupon traits 

abstract type FinalCoupon end 

struct MaturingCoupon <: FinalCoupon end 
struct NoMaturingCoupon <: FinalCoupon end 

has_maturing_coupon(::Bond) = MaturingCoupon()
has_maturing_coupon(::BondCE2012) = NoMaturingCoupon()
has_maturing_coupon(::FloatingRateBond) = MaturingCoupon()