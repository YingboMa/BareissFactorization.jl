module BareissFactorization

using LinearAlgebra, SparseArrays

export left_looking_bareiss

function interchange_rows!(A, i, i′)
    i == i′ && return nothing
    for j in axes(A, 2)
        A[i′, j], A[i, j] = A[i, j], A[i′, j]
    end
    return nothing
end

function find_pivot(x, k)
    idx = findfirst(!iszero, view(x, k:length(x)))
    idx === nothing && return nothing
    idx + k - 1
end

function left_looking_bareiss(A; find_pivot = find_pivot)
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
        i = find_pivot(x, k)
        if i === nothing
            dc = dp
            if k <= m
                L[k, k] = dp
            end
        else
            interchange_rows!(L, i, k)
            interchange_rows!(P, i, k)
            interchange_rows!(x, i, k)
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
