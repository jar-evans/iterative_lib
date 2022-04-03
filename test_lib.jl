using BenchmarkTools

include("iterative_lib.jl")



# initialize the matrix
A =  [10. -1. 2. 0.;-1. 11. -1. 3.;2. -1. 10. -1.;0.0 3. -1. 8.]
# initialize the RHS vector
b = [6.; 25.; -11.; 15.]
# initialise the first guess
x = [0.; 0.; 0.; 0.]

w = 1

x = @btime gauss_seidal(x, A, b)

println(x)