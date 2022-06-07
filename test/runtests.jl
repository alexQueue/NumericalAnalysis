using Test
using LinearAlgebra
include("../src/Beam1D.jl")


@testset "Parameter validations" begin
  pars = (mu=x->1, EI=x->1, q=x->-(x<0.5))

  BCs  = Dict((0,'H')=>1,
              (0,'G')=>2,
              (1,'M')=>3,
              (1,'Q')=>4)
  n = 10
  grid = collect(LinRange(0,1,n))
  h = grid[2] - grid[1]

  @test_throws AssertionError Beam1D.Problem(
    pars,
    Dict((1,'M')=>3,
         (1,'Q')=>4
        ),
    grid
  )
  @test_throws AssertionError Beam1D.Problem(
    pars,
    Dict((0,'H')=>1,
        (0,'G')=>2,
        (1,'H')=>2,
        (1,'G')=>2,
        (1,'M')=>3,
        (1,'Q')=>4
    ),
    grid
  )

  # Successful path
  prob = Beam1D.Problem(
    pars,
    Dict((0,'H')=>1,
        (0,'G')=>2,
        (1,'M')=>3,
        (1,'Q')=>4
    ),
    grid
  )

  @test prob.BCs[0, 'H'] == 1
end

@testset "Static build integration tests" begin
  pars = (mu=x->1, EI=x->1, q=x->-(x<0.5))

  BCs  = Dict((0,'H')=>1,
              (0,'G')=>2,
              (1,'M')=>3,
              (1,'Q')=>4)
  n = 10
  grid = collect(LinRange(0,1,n))
  h = grid[2] - grid[1]

  sys   = Beam1D.build(Beam1D.Problem(pars,BCs,grid))
  static_sol, u = Beam1D.solve_st(sys)

  # x_0 is correct for BCs
  @test u[1] ≈ 1
  # x' is correct for BCs
  @test u[2] ≈ 2
  
  # These two are not working correctly.
  
  # x'''
  @test dot(u[1:4], -[6/h^2,   4/h,     -6/h^2,  2/h]) ≈ 3 broken=true
  # x''''
  @test dot(u[1:4], -[12/h^3,   6/h^2,   -12/h^3, 6/h^2]) ≈ 4 broken=true
end