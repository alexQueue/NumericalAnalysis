module Diagnostics
  import LinearAlgebra

  function print_diagnostics(u::Vector{Float64}, h::Float64)
    println("At position 0:")
    println(u[1:4])
    println(string("X is ", u[1]))
    println(string("X' is ", u[2]))
    println(string("X'' is ", LinearAlgebra.dot(u[1:4], -[6/h^2,   4/h,     -6/h^2,  2/h])))
    println(string("X''' is ", LinearAlgebra.dot(u[1:4], -[12/h^3,   6/h^2,   -12/h^3, 6/h^2])))

    println("-------------")
    println("At position L:")
    println(string("X is ", u[end-1]))
    println(string("X' is ", u[end]))
    println(string("X'' is ", LinearAlgebra.dot(u[end-3:end], [6/h^2,   2/h,  -6/h^2,     4/h])))
    println(string("X''' is ", LinearAlgebra.dot(u[end-3:end], [12/h^3,     6/h^2, -12/h^3,     6/h^2])))
  end
end