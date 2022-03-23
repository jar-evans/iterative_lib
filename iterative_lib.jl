using LinearAlgebra
using BenchmarkTools

const MAX_ITER = 100

function norm(x)
    norm = 0
    for i in eachindex(x)
        @inbounds norm += x[i]^2
    end
    norm ^0.5
end

function jacobi(A, b, x)
    
    N = size(A, 1)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)        
    
        for j in 1:N

            L = A[j, 1:j-1]'*x[1:j-1]
            U = A[j, j+1:N]'*x[j+1:N]

            x_new[j] = (b[j] - L - U) / A[j, j]

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        norm(x - x_new) < 1e-10 && break
    
        x = x_new
        i += 1

    end

    x

end

function gauss_seidal(A, b, x)

    N = size(A, 1)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)        
    
        for j in 1:N

            L = A[j, 1:j-1]'*x_new[1:j-1]
            U = A[j, j+1:N]'*x[j+1:N]

            x_new[j] = (b[j] - L - U) / A[j, j]

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        norm(x - x_new) < 1e-10 && break
    
        x = x_new
        i += 1

    end

    x

end

function weighted_jacobi(A, b, x, w)

    N = size(A, 1)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)        
    
        for j in 1:N

            L = A[j, 1:j-1]'*x[1:j-1]
            U = A[j, j+1:N]'*x[j+1:N]
            W = (1 - w) * x[j]

            x_new[j] = (w * (b[j] - L - U) / A[j, j]) + W

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        norm(x - x_new) < 1e-10 && break
    
        x = x_new
        i += 1

    end

    x

end

function SOR(A, b, x, w)

    N = size(A, 1)
    i = 0

    while i < MAX_ITER

        x_new = zeros(N)

        for j in 1:N

            L = A[j, 1:j-1]'*x_new[1:j-1]
            U = A[j, j+1:N]'*x[j+1:N]
            W = (1 - w) * x[j]

            x_new[j] = (w * (b[j] - L - U) / A[j, j]) + W

            (j > 1) && x_new[j] == x_new[j-1] && break

        end

        norm(x - x_new) < 1e-10 && break
    
        x = x_new
        i += 1

    end

    x

end

function CG(A, b, x)

    N = size(A, 1)

    r = b-A*x;
    p = r;
    r_old = r'*r;

    i = 0

    while i < MAX_ITER

        Ap = A*p
        alpha = r_old/(p'*Ap)

        x = x+alpha*p
        r = r-alpha*Ap

        r_new = r'*r

        r_new^0.5 < 1e-10 && break

        p = r+(r_new/r_old)*p
        r_old = r_new

        i += 1

    end

    x

end
