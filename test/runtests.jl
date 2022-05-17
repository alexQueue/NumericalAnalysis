using Test
using LinearAlgebra
include("../Beam1D.jl")


@testset "Parameter validations" begin
  @test_throws AssertionError Beam1D.make_BC_from_dict(Dict("x_0"=>1, "xprime_0" => 1, "Q_0" => 0, "Q_L" => 0, "M_0" => 0))
  @test_throws AssertionError Beam1D.make_BC_from_dict(Dict("x_0"=>1))
  # @test Main.Beam1D.BoundaryConditions(1.0, 1.0, nothing, 0.0, nothing, nothing, nothing, 0.0) == Main.Beam1D.BoundaryConditions(1.0, 1.0, nothing, 0.0, nothing, nothing, nothing, 0.0)
  # @test Beam1D.make_BC_from_dict(Dict("x_0"=>1, "xprime_0" => 1, "Q_0" => 0, "Q_L" => 0)) == Beam1D.BoundaryConditions(1, 1, nothing, 0, nothing, nothing, nothing, 0)
end

@testset "Static build integration tests" begin
  L_0 = 0.0
  L = 1.0
  h = 0.01

  # BC using dictionary instead
  BC = Dict("x_0"=>1,"xprime_0"=>2,"M_0"=>3,"Q_0"=>4)
  BoundaryConditions = Beam1D.make_BC_from_dict(BC)

  x_grid = collect(L_0:h:L)
  q(x) = 1
  EI(x) = 1
  mu(x) = 1

  par = Beam1D.Parameters(mu,EI,q,BoundaryConditions)
  sys = Beam1D.build(x_grid,par)
  sol, u = Beam1D.solve_st(sys)


  # x_0 is correct for BCs
  @test u[1] ≈ 1
  # x' is correct for BCs
  @test u[2] ≈ 2
  
  # These two are not working correctly.
  # x'''
  @test dot(u[1:4], -[6/h^2,   4/h,     -6/h^2,  2/h]) ≈ 3
  # # x''''
  @test dot(u[1:4], -[12/h^3,   6/h^2,   -12/h^3, 6/h^2]) ≈ 4
end

touch(joinpath(ENV["HOME"], "julia-runtest"))