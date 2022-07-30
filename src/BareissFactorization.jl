module BareissFactorization

using Base: require_one_based_indexing

using LinearAlgebra, SparseArrays

export left_looking_bareiss

# Compute A D \ b
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

# Compute y = A D x
function gaxpy!(A, ds, x::AbstractVector, y::AbstractVector, tmp)
    require_one_based_indexing(A, x, y, tmp)
    m, n = size(A)
    n == 0 && return y
    if (m != length(y) || n != length(x) || n != length(ds) || length(ds) > length(tmp))
        throw(DimensionMismatch(""))
    end
    p = first(ds)
    for i in 1:min(2, n)
        tmp[i] = one(p)
    end
    for i in 3:n
        tmp[i] = p
        p *= ds[i - 1]
    end
    p_nolast = n == 1 ? one(p) : p
    p = n == 1 ? one(p) : ds[n-1]
    for i in n-2:-1:1
        tmp[i] *= p
        if i != 1
            p *= ds[i]
        end
    end
    #=
    ll = similar(ds)
    for i in 1:n
        p = one(p)
        for j in 1:n-1
            (j == i || j == i - 1) && continue
            p *= ds[j]
        end
        ll[i] = p
    end
    @assert ll == @view tmp[1:n]
    tmp = ll
    =#
    for j in 1:n-1
        xj = -tmp[j] * x[j]
        for i in axes(A, 1)
            y[i] += A[i, j] * xj
        end
    end
    pl = ds[end]
    xn = tmp[n] * x[n]
    for i in eachindex(y)
        # This overflows quite often
        og = y[i] * pl - A[i, n] * xn
        yi, r = divrem(og, p_nolast)
        @assert iszero(r)
        y[i] = yi
    end
    y
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
    tmp = zeros(T, max(0, m - 1))
    ps = zeros(T, m)
    dp = 1
    for k in 1:n
        Pb = P * A[:, k]
        #Ll = Rational{T}[L[:, 1:k-1] * Diagonal(@view ds[1:k-1]) [zeros(k-1, m-k+1); I(m-k+1)]]
        #for j in k:m
        #    Ll[j, j] = 1//dp
        #end
        #x = T.(Ll \ (Pb))
        x1 = scaled_ldiv!(LowerTriangular(L[1:k-1, 1:k-1]), ps, Pb[1:k-1])
        #x2′ = Pb[k:m] * dp - Int.(L[k:m, 1:k-1] * Diagonal(@view ds[1:k-1]) * x1 * dp)
        L21, p_ps, y, tmp = (@view L[k:m, 1:k-1]), (@view ps[1:k-1]), zeros(T, m - k + 1), tmp
        x2 = dp * Pb[k:m] + gaxpy!(L21, p_ps, x1, y, tmp)

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
