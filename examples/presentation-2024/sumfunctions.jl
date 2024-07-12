function fast_sum(x::Matrix)
    s = zero(eltype(x))
    m, n = size(x)
    @inbounds @fastmath @simd for j in 1:n
        for i in 1:m
            s += x[i,j]
        end
    end
    return s
end



function slow_sum(x::Matrix)
    s = zero(eltype(x))
    m, n = size(x)
    @inbounds @fastmath @simd for i in 1:m
        for j in 1:n
            s += x[i,j]
        end
    end
    return s
end



precompile(fast_sum, (Matrix{Float64},));

precompile(slow_sum, (Matrix{Float64},));