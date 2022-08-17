function push_range!(collection::Vector{Float64}, vec::Vector{Float64})
    [push!(collection, i) for i in vec]
end

function pushfirst_range!(collection::Vector{Float64}, vec::Vector{Float64})
    [pushfirst!(collection, i) for i in reverse(vec)]
end