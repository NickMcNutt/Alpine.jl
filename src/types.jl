typealias AtomID Int
typealias Group Set{Int}
typealias Indices Vector{Int}
typealias Coords Matrix{Float64}
typealias AbstractIndices AbstractVector{Int}
typealias AbstractCoords{T} AbstractMatrix{T}

type Cluster
    num::Int
    max_num::Int
    indices::Indices
    r::Vector{Float64}
    coords::Coords


    Cluster(n::Int) = new(zero(Int), n, Indices(n), Vector{Float64}(n), Coords(3, n))
end

typealias Clusters Dict{Int, Cluster}
