module BareissFactorization

using Base: require_one_based_indexing

using LinearAlgebra, SparseArrays

export left_looking_bareiss

function scaled_ldiv!(A::LowerTriangular, ds, b)
    require_one_based_indexing(A, b)
    isempty(b) && return b
    L = A.data
    n = size(L, 2) # LowerTriangular guarantees squareness.
    if !(n == length(b)) || n > length(ds)
        throw(DimensionMismatch("size(L) = $(size(L)), size(b) = $(size(b)), length(ds) = $(length(ds))"))
    end
    p = one(eltype(ds))
    @inbounds for i in 2:n
        d = ds[i-1]
        p *= d
        b[i] *= p
    end
    p = one(eltype(ds))
    @inbounds for j in 1:n
        if j > 2
            p *= ds[j - 2]
        end
        #iszero(L[j, j]) && throw(SingularException(j))
        xj, r = divrem(b[j], p)
        @assert iszero(r)
        b[j] = xj
        xj *= p
        for i in j+1:n
            b[i] -= L[i, j] * xj
            if i != n
                xj *= ds[i]
            end
        end
    end
    return b
end

function interchange_rows!(A, i, i′)
    i == i′ && return nothing
    for j in axes(A, 2)
        A[i′, j], A[i, j] = A[i, j], A[i′, j]
    end
    return nothing
end

find_pivot(x) = findfirst(!iszero, x)

function left_looking_bareiss(A; find_pivot = find_pivot)
    Base.require_one_based_indexing(A)
    m, n = size(A)
    P = Matrix(I(m))
    T = eltype(A)
    L = zeros(T, m, m)
    U = zeros(T, m, n)
    ds = zeros(Rational{T}, m)
    ps = zeros(T, m)
    dp = 1
    for k in 1:n
        #Ll = Rational{T}[L[:, 1:k-1] * Diagonal(@view ds[1:k-1]) [zeros(k-1, m-k+1); I(m-k+1)]]
        #for j in k:m
        #    Ll[j, j] = 1//dp
        #end
        Pb = P * A[:, k]
        #x = T.(Ll \ (Pb))
        x1 = scaled_ldiv!(LowerTriangular(L[1:k-1, 1:k-1]), ps, Pb[1:k-1])
        x2 = Pb[k:m] * dp - Int.(L[k:m, 1:k-1] * Diagonal(@view ds[1:k-1]) * x1 * dp)

        U[1:k-1, k] = x1
        i = find_pivot(x2)
        if i === nothing
            dc = dp
            if k <= m
                L[k, k] = dp
            end
        else
            i′ = i + k - 1
            interchange_rows!(L, i′, k)
            interchange_rows!(P, i′, k)
            interchange_rows!(x2, i, 1)
            dc = L[k, k] = U[k, k] = x2[1]
            L[k+1:m, k] = x2[2:m-k+1]
        end
        if k <= m
            ds[k] = 1//(dp * dc)
            ps[k] = dc
        end
        dp = dc
    end
    P, L, Diagonal(ds), U
end

end
