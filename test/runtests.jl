using BareissFactorization
using SparseArrays, LinearAlgebra
using Test

@testset "Left looking Bareiss" begin
    for i in 1:1000
        for m in 1:9, n in 1:9
            A = rand([-2, -1, 0, 0, 0, 1, 2, 3], 5, 4)
            P, L, D, U = left_looking_bareiss(A)
            isvalid = P * A == Int.(L * D * U)
            if !isvalid
                @info "A = " A
            end
            @test isvalid
        end
    end
end
