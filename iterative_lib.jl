using LinearAlgebra
using StaticArrays
using BenchmarkTools

const MAX_ITER = 100
const EPS = 1e-10


function jacobi(A, b, x)

    N = length(b)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)

        @inbounds @views for j = 1:N
            if j == 1
                L = 0.0
                U = A[j, j+1:N]' * x[j+1:N]
            elseif j == N
                L = A[j, 1:j-1]' * x[1:j-1]
                U = 0.0
            else
                L = A[j, 1:j-1]' * x[1:j-1]
                U = A[j, j+1:N]' * x[j+1:N]
            end

            x_new[j] = (b[j] .- U .- L) / A[j, j]

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        e = x - x_new
        (e' * e)^0.5 < EPS && break

        x .= x_new
        i += 1

    end
    x

end

function gauss_seidal(A, b, x)

    N = length(b)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)

        @inbounds @views for j = 1:N
            if j == 1
                L = 0
                U = A[j, j+1:N]' * x[j+1:N]
            elseif j == N
                L = A[j, 1:j-1]' * x_new[1:j-1]
                U = 0
            else
                L = A[j, 1:j-1]' * x_new[1:j-1]
                U = A[j, j+1:N]' * x[j+1:N]
            end

            x_new[j] = (b[j] - U - L) / A[j, j]

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        norm(x - x_new) < EPS && break

        x .= x_new
        i += 1

    end

    x

end

function weighted_jacobi(A, b, x, w)

    N = length(b)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)

        @inbounds @views for j = 1:N
            if j == 1
                L = 0
                U = A[j, j+1:N]' * x[j+1:N]
            elseif j == N
                L = A[j, 1:j-1]' * x[1:j-1]
                U = 0
            else
                L = A[j, 1:j-1]' * x[1:j-1]
                U = A[j, j+1:N]' * x[j+1:N]
            end

            W = (1 - w) * x[j]

            x_new[j] = (w * (b[j] - L - U) / A[j, j]) + W

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        norm(x - x_new) < EPS && break

        x .= x_new
        i += 1

    end

    x

end

function SOR(A, b, x, w)

    N = length(b)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)

        @inbounds @views for j = 1:N
            if j == 1
                L = 0
                U = A[j, j+1:N]' * x[j+1:N]
            elseif j == N
                L = A[j, 1:j-1]' * x_new[1:j-1]
                U = 0
            else
                L = A[j, 1:j-1]' * x_new[1:j-1]
                U = A[j, j+1:N]' * x[j+1:N]
            end

            W = (1 - w) * x[j]

            x_new[j] = (w * (b[j] - L - U) / A[j, j]) + W

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        norm(x - x_new) < EPS && break

        x .= x_new
        i += 1

    end

    x

end

function CG(A, b, x)

    N = length(b)

    r = b - A * x
    p = r
    r_old = r' * r

    i = 0

    while i < MAX_ITER

        Ap = A * p
        alpha = r_old / (p' * Ap)

        x = x + alpha * p
        r = r - alpha * Ap

        r_new = r' * r

        r_new^0.5 < EPS && break

        p = r + (r_new / r_old) * p
        r_old = r_new

        i += 1

    end

    x

end
