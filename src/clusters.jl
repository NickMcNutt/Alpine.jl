#allocate_clusters(indices::Indices, atom_type_counts::Dict{Symbol, Int}) = Dict(Pair{Symbol, Clusters}[t => Clusters(n, indices) for (t, n) in atom_type_counts])

type Cluster{T<:Real}
    id::AtomID
    ids::Vector{AtomID}
    coords::Matrix{T}
end


