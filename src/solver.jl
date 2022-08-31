
# Main solver -- iterates on prices and values until convergence or max_iters is reached. 
function solve!(m, a;
    method = DivideAndConquer(), max_iters = 10_000, g = 1.0, err = 1e-10, print_every = 50
)
    i = 0
    converged = false
    while true
        update_Ev!(m, a)
        update_all_vD!(m, a)
        update_v_and_policies!(m, a; method)
        update_q!(m, a)
        swap!(m, a, g)

        dis = distances(m, a)
        (max(dis...) < err) && (converged = true)

        if print_every != 0 && (mod(i, print_every) == 0 || (i + 1 >= max_iters)) || converged
            println("$(i+1): $dis")
        end
        converged && (println("Converged."); break)
    
        i += 1
        (i >= max_iters) && (@error("Did not converged"); break)
    end 
end 

solve!(a; kwargs...) = solve!(a.model, a; kwargs...)


# Updates the cotinuation under default 
function update_Ev!(_, a)  
    Ev, v = get_Ev(a), get_v(a)
    model = get_base_pars(a)
    mul!(Ev, v, get_y(model).Π)
end 

# Updates the continuation value under default for the floating rate bond model. 
# The difference here is that the floating rate bond carries two exogenous staets (yi, yi0)
function update_Ev!(::AbstractModel{B}, a) where B<:AbstractFloatingRateBond
    Ev, v = get_Ev(a), get_v(a)
    model = get_base_pars(a)
    Π = get_y(model).Π
    @batch for idx in CartesianIndices((axes(v,1), axes(v,2)))
        bi, yi = Tuple(idx)
        runningsum = zero(eltype(Ev))
        for yiprime in axes(v, 2)
            @inbounds runningsum += v[bi, yiprime, yi] * Π[yiprime, yi]
        end
        @inbounds Ev[bi,yi] = runningsum
    end
end


# Updates the default values  
function update_all_vD!(model, a)
    (; vD, Ev, ucdef, cdef) = get_cache(a)
    vD1_new = get_new(a).vD1
    vD1 = get_current(a).vD1

    (; β, u) = get_preferences(model)
    izero = get_zero_index(get_bond(model))
    θ = get_base_pars(model).def_costs.reentry
    Π = get_y(model).Π
    m_max = get_m(model).m_max

    # updates the next period default value 
    @batch for yi in eachindex(get_y_grid(model))
        @inbounds vD1_new[yi] = ucdef[yi] + β * θ * Ev[izero, yi] + β * (1 - θ) * dot(view(Π, :, yi), vD1)
    end

    # updates the current value -- note the role of the m-shock here
    @batch for yi in eachindex(get_y_grid(model))
        @inbounds vD[yi] = u(cdef[yi] - m_max) + β * θ * Ev[izero, yi] + β * (1 - θ) * dot(view(Π, :, yi), vD1_new)
    end
end 


# Main function that optimizes the bellman equation updating the values and policies 
function update_v_and_policies!(m, a; method = DivideAndConquer()) 
    @batch for yi in axes(get_v(a), 2)
        bellman_given_yi!(m, a, (yi,),  method)
    end
end 


# Specializes to the Floating Rate Bond model. The key differences is that the state 
# now contains not just current yi but also previous period yi
function update_v_and_policies!(m::AbstractModel{B}, a; method = DivideAndConquer()) where B<:AbstractFloatingRateBond
    v = get_v(a)
    Threads.@threads for x in CartesianIndices((axes(v, 2), axes(v,3)))
        bellman_given_yi!(m, a, Tuple(x),  method)
    end
end 


# Updates the prices
function update_q!(model, a)
    q_new = get_new(a).q
    d_pol = get_d_pol(a)
    bond_return = get_bond_return(a)

    Π = get_y(model).Π
    m_min = get_m(model).m_min
    
    @batch for idx in CartesianIndices((axes(q_new, 1), axes(q_new, 2)))
        bi, yi = Tuple(idx)
        @inbounds begin 
            tmp = zero(eltype(q_new))
            for yiprime in axes(q_new, 2)
                mD = d_pol[bi, yiprime]
                if mD > m_min 
                    tmp += bond_return[bi, yiprime] * Π[yiprime, yi]
                end
            end
            q_new[bi, yi] = tmp / get_R(model)
        end 
    end
end


# Updates the prices for the Floating rate bond model. In this case, the 
# coupon rate is also updated. 
function update_q!(model::AbstractModel{B}, a) where B<:AbstractFloatingRateBond
    q_new = get_new(a).q
    κ_new = get_new(a).κ
    d_pol = get_d_pol(a)
    bond_return = get_bond_return(a)
    repay = get_cache(a).repay

    κbar = get_κbar(get_bond(model))
    m_min = get_m(model).m_min
    R = get_R(model)
    Π = get_y(model).Π
 
    @batch for idx in CartesianIndices((axes(κ_new, 1), axes(κ_new, 2)))
        bi, yi = Tuple(idx)
        q_S = zero(eltype(q_new))
        tmp = zero(eltype(q_new))
        @inbounds for yiprime in axes(κ_new, 2)
            mD = d_pol[bi, yiprime, yi]
            if mD > m_min 
                q_S += repay[bi, yiprime, yi] * Π[yiprime, yi]
                tmp += bond_return[bi, yiprime, yi] * Π[yiprime, yi]
            end
        end
        q_S = q_S / R
        @inbounds κ_new[bi, yi] = min((1/q_S - 1), κbar)  # coupon on maturing bond as well
        @inbounds q_new[bi, yi] = tmp / R
    end
end


# Swap the old values for the new ones, and smooth with factor g
function swap!(m, a, g)
    cur = get_current(a)
    new = get_new(a)

    if g == 1 
        a.current, a.new = a.new, a.current
    else
        cur.q .= (1 - g) .* cur.q .+ g .* new.q
        cur.v .= (1 - g) .* cur.v .+ g .* new.v
        cur.vD1 .= (1 - g) .* cur.vD1 .+ g .* new.vD1
    end 
end 


# Specializes to the Floating Rate Bond Model 
function swap!(::AbstractModel{B}, a, g) where B<:AbstractFloatingRateBond
    cur = get_current(a)
    new = get_new(a)

    if g == 1 
        a.current, a.new = a.new, a.current
    else
        cur.q .= (1 - g) .* cur.q .+ g .* new.q
        cur.v .= (1 - g) .* cur.v .+ g .* new.v
        cur.vD1 .= (1 - g) .* cur.vD1 .+ g .* new.vD1
        cur.κ .= (1 - g) .* cur.κ .+ g .* new.κ
    end 
end 


function distances(m, a)
    cur, new = get_current(a), get_new(a)
    s1, s2, s3 = zero(eltype(cur.v)), zero(eltype(cur.q)), zero(eltype(cur.vD1))

    @tturbo for j in eachindex(cur.v, cur.q)
        s1 = max(s1, abs(cur.v[j] - new.v[j]))
        s2 = max(s2, abs(cur.q[j] - new.q[j]))
    end

    @tturbo for j in eachindex(cur.vD1, new.vD1)
        s3 = max(s3, abs(cur.vD1[j] - new.vD1[j]))
    end

    return (v = s1, q = s2, vD = s3)
end 


function distances(::AbstractModel{B}, a) where B<:AbstractFloatingRateBond
    cur, new = get_current(a), get_new(a)
    s1, s2, s3 = zero(eltype(cur.v)), zero(eltype(cur.q)), zero(eltype(cur.κ))
    s4 = zero(eltype(cur.vD1))


    @tturbo for j in eachindex(cur.v, new.v)
        s1 = max(s1, abs(cur.v[j] - new.v[j]))
    end

    @tturbo for j in eachindex(cur.κ, new.κ)
        s2 = max(s2, abs(cur.κ[j] - new.κ[j]))
        s3 = max(s3,abs(cur.q[j] - new.q[j]))
    end

    @tturbo for j in eachindex(cur.vD1, new.vD1)
        s4 = max(s4, abs(cur.vD1[j] - new.vD1[j]))
    end
    return (v = s1, κ = s2, q = s3, vD = s4)
end 

