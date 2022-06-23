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

  sys   = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
  static_sol = Beam1D.solve_st(sys)
  @test !isnothing(static_sol)
end