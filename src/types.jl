import Base: length, size, sizeof, getindex

#Indices = Vector{Int}
AtomID = Int

type Container{T}
    object::T
    items::Int
    size::Int
end

#==Container{T<:Real}(::Type{Matrix{T}}, n::Int) = Container(Matrix{T}(3, n), 0, n)
getindex{T<:Real}(c::Container{Matrix{T}}, i::Int) = c.object[:, i]

Container(Cluster, num_li * num_frames * num_neighbor_types)==#

immutable Clusters
    num::Int
    indices_central::Indices

    num_neighbors::Int
    indices_neighbors::Matrix{Int}

    dist::Matrix{Float64}
    coords::Array{Float64, 3}

    function Clusters(num_neighbors::Int, num::Int)
        nn, nc = num_neighbors, num
        new(nc, Vector{Int}(nc),
            nn, Matrix{Int}(nn, nc),
            Matrix{Float64}(nn, nc),
            Array{Float64}(3, nn, nc))
    end

    function Clusters(num_neighbors::Int, indices_central::Indices)
        nn, nc = num_neighbors, length(indices_central)
        new(nc, copy(indices_central),
            nn, Matrix{Int}(nn, nc),
            Matrix{Float64}(nn, nc),
            Array{Float64}(3, nn, nc))
    end
end

length(clusters::Clusters) = clusters.num_neighbors * clusters.num
size(clusters::Clusters) = clusters.num_neighbors, clusters.num
sizeof(clusters::Clusters) = sum(x -> sizeof(clusters.(x)), fieldnames(Clusters))
