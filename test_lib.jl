using BenchmarkTools

include("iterative_lib.jl")



# initialize the matrix
A =  [10. -1. 2. 0.;-1. 11. -1. 3.;2. -1. 10. -1.;0.0 3. -1. 8.]
# initialize the RHS vector
b = [6.; 25.; -11.; 15.]
# initialise the first guess
x = [0.; 0.; 0.; 0.]

w = 1

# x1 = @btime jacobi($A, $b, $x)
# x2 = @btime gauss_seidal($A, $b, $x)
# x3 = @btime weighted_jacobi($A, $b, $x, $w)
# x4 = @btime SOR($A, $b, $x, $w)
# println(x1)
# println(x2)
# println(x3)
# println(x4)
# println(x5)
# println("Solution:")
# println(x)
# error = A\b - x
# println("Error:")
# println(error)
# Profile.print(maxdepth=11)

x = @btime A\b
println(x)

# println(x)
