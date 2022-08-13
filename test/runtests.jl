using Test
using LTBonds


@testset verbose=true "preferences" begin
    u = CRRA2() 
    @test u(2.0) == - 1 / 2.0    
    @test inv_u(u, 2.0) == - 1 / 2.0    

    # generic version of CRRA
    u_float = CRRA(2.0)  
    @test u(2.0) == - 1 / 2.0    
    @test inv_u(u, 2.0) == - 1 / 2.0    

    p = Preferences(β = 0.9, u = u)
    @test p.β == 0.9
    @test u(2.0) == - 1 / 2.0

    @test LTBonds.find_m_root(u, 1.0, 2.0, 1.0) ≈ 0.3819660112501051

    # testing with the generic version 
    @test LTBonds.find_m_root(u_float, 1.0, 2.0, 1.0) ≈ LTBonds.find_m_root(u, 1.0, 2.0, 1.0) 

    # checking that log solution matches approximately the CRRA one
    @test LTBonds.find_m_root(LTBonds.Log(), 1.0, 2.0, 1.0) ≈ 0.41802329313067355
end


@testset verbose=true "tauchen" begin 

    out = LTBonds.tauchen((std = 2.0, ρ = 0.8, μ = 0.0, span = 2.0, n = 3, tails = false)) 
    @test out[2][2] ≈ 0.1602205276394577
    @test out[2][end] ≈ 0.8397720561441128

end


@testset verbose=true "long bonds" begin

    model = let
        R = 1.01
        β = 0.9540232420
        pref = Preferences(β = β, u = make_CRRA(ra = 2))
        y = discretize(YProcess(n = 20, ρ = 0.948503, std = 0.027092, μ = 0.0, span = 3.0, tails = false))
        m = MTruncatedNormal(; std = 0.01, span = 2.0, quadN = 100)
        penalty = DefCosts(pen1 = -0.1881927550, pen2 = 0.2455843389, quadratic = true, reentry = 0.0385)
        bond = Bond(n = 350, min = 0.0, max = 1.5, κ = R - 1, λ = 0.05)  
    
        LTBondModel(
            y = y,
            m = m, 
            preferences = pref, 
            bond = bond, 
            def_costs = penalty, 
            R = R,
        )
    end

    ws = generate_workspace(model);
    solve!(ws; err = 1e-10, print_every = 0)

    @test get_q(ws)[100, 10] ≈ 0.8684817791047236
    @test get_v(ws)[100, 10] ≈ -22.008065357421078
    @test get_vD(ws)[10] ≈ -22.56604433683637
end