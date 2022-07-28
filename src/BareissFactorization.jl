module BareissFactorization

using LinearAlgebra, SparseArrays

export left_looking_bareiss

function left_looking_bareiss(A)
    Base.require_one_based_indexing(A)
    m, n = size(A)
    P = Matrix(I(m))
    T = eltype(A)
    L = zeros(T, m, m)
    U = zeros(T, m, n)
    ds = zeros(Rational{T}, m)
    dp = 1
    for k in 1:n
        Ll = Rational{T}[L[:, 1:k-1] * Diagonal(@view ds[1:k-1]) [zeros(k-1, m-k+1); I(m-k+1)]]
        for j in k:m
            Ll[j, j] = 1//dp
        end
        x = T.(Ll \ (P * A[:, k]))
        U[1:k-1, k] = @view x[1:k-1]
        i = findfirst(!iszero, view(x, k:m))
        if i === nothing
            dc = dp
            if k <= m
                L[k, k] = dp
            end
        else
            i += k - 1
            L[i, :], L[k, :] = L[k, :], L[i, :]
            P[i, :], P[k, :] = P[k, :], P[i, :]
            x[i], x[k] = x[k], x[i]
            dc = L[k, k] = U[k, k] = x[k]
            L[k+1:m, k] = x[k+1:m]
        end
        if k <= m
            ds[k] = 1//(dp * dc)
        end
        dp = dc
    end
    P, L, Diagonal(ds), U
end

end
