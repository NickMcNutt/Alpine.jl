using VPTrees

function ndf{T}(metric::Function, point_indices::Vector{Int}, radius::T, Δr::T)
    num_points::Int = length(point_indices)
    num_bins::Int = ceil(Int, radius / Δr)
    bins = zeros(UInt64, num_bins)

    tree = VPTree(metric, point_indices)

    neighbors = NeighborList{T}(num_points)
    @inbounds for i in 1:num_points
        rangesearch!(neighbors, tree, i, radius)

        k = length(neighbors)
        for j in 1:k
            i == neighbors.indices[j] && continue
            r::T = neighbors.distances[j]
            b::Int = floor(Int, r / Δr) + 1
            bins[b] += 1
        end
    end

    return bins
end

function ndf{T}(pair_count::Function, metric::Function, point_indices::Vector{Int}, radius::T, Δr::T)
    num_points::Int = length(point_indices)
    tree = VPTree(metric, point_indices)
    neighbors = NeighborList{T}(num_points)

    @inbounds for i in 1:num_points
        rangesearch!(neighbors, tree, i, radius)

        k = length(neighbors)
        for j in 1:k
            ni = neighbors.indices[j]
            i == ni && continue
            r::T = neighbors.distances[j]
            b::Int = floor(Int, r / Δr) + 1
            pair_count(b, i, ni)
        end
    end
end

function distance{T}(coords::Matrix{T}, xbw::T, ybw::T, zbw::T)
    xhbw = xbw / 2
    yhbw = ybw / 2
    zhbw = zbw / 2

    function(i::Int, j::Int)
        @inbounds x = coords[1, i] - coords[1, j]
        @inbounds y = coords[2, i] - coords[2, j]
        @inbounds z = coords[3, i] - coords[3, j]

        if x < -xhbw x += xbw elseif x > xhbw x -= xbw end
        if y < -yhbw y += ybw elseif y > yhbw y -= ybw end
        if z < -zhbw z += zbw elseif z > zhbw z -= zbw end

        return sqrt(x^2 + y^2 + z^2)
    end
end

function rdf{T}(atoms::Atoms, xbw::T, ybw::T, zbw::T, radius::T, Δr::T)
    num_bins::Int = ceil(Int, radius / Δr)
    bins = Matrix{T}(num_bins, 2)

    metric = distance(atoms[:coords], xbw, ybw, zbw)
    num_atoms = length(atoms)
    point_indices = collect(1:num_atoms)
    int_bins = ndf(metric, point_indices, radius, Δr)

    for i in 1:num_bins
        r_mid = i*Δr - Δr/2
        divisor = (num_atoms * 4π * r_mid * r_mid * Δr) * (num_atoms / (xbw * ybw * zbw))
        bins[i, 1] = r_mid
        bins[i, 2] = int_bins[i] / divisor
    end

    return bins
end

function rdf{T}(f::Function, nf::Int, atoms::Atoms, xbw::T, ybw::T, zbw::T, radius::T, Δr::T)
    num_bins::Int = ceil(Int, radius / Δr)
    bins = Matrix{T}(nf, num_bins)
    int_bins = zeros(UInt64, nf, num_bins)

    metric = distance(atoms[:coords], xbw, ybw, zbw)
    num_atoms = length(atoms)
    point_indices = collect(1:num_atoms)
    ndf((b, i, j) -> f(atoms, int_bins, b, i, j), metric, point_indices, radius, Δr)

    for i in 1:num_bins
        r_mid = i*Δr - Δr/2
        divisor = (num_atoms * 4π * r_mid^2 * Δr) * (num_atoms / (xbw * ybw * zbw))
        for g in 1:nf
            bins[g, i] = int_bins[g, i] / divisor
        end
    end

    return bins
end

@inline function rdf{T}(pair_count::Function, nf::Int, frame::Frame, radius::T, Δr::T)
    rdf(pair_count, nf, frame[:atoms], box_dims(frame)..., radius, Δr)
end

function rdf{T}(pair_count::Function, nf::Int, frames::Vector{Frame}, radius::T, Δr::T)
    mean(pmap(frames) do frame
        rdf(pair_count, nf, frame[:atoms], box_dims(frame)..., radius, Δr)
    end)
end

function rdf_axis{T}(radius::T, Δr::T)
    num_bins = ceil(Int, radius / Δr)::Int
    T[i*Δr - Δr/2 for i in 1:num_bins]
end

function rdf{N, T}(frames::Vector{Frame}, element_types::NTuple{N, Symbol}, max_radius::T, Δr::T)
    list = map(frames) do frame
        indices = findin(frame[:types], element_types)
        c = frame[indices][:coords]

        (box_dims(frame)..., c, max_radius, Δr)
    end

    mean(pmap(x -> rdf(x...), list))
end

function rdf{N, T}(frames::AbstractArray{Int}, types::NTuple{N, Symbol}, bw::Vector{Tuple{T, T, T}}, max_radius::T, Δr::T, num_atoms::Int, atom_coords::AbstractMatrix{T}, atom_types::Vector{Symbol})
    list = map(frames) do f
        indices = get_atoms(f, types, atom_types, num_atoms)
        coords = atom_coords[:, indices]
        
        (bw[f][1], bw[f][2], bw[f][3], coords, max_radius, Δr)
    end
    
    mean(pmap(x -> rdf(x...), list))
end

function rdf{N, T}(frames::AbstractArray{Int}, types::NTuple{N, Symbol}, box_width::T, max_radius::T, Δr::T, num_atoms::Int, atom_coords::AbstractMatrix{T}, atom_types::Vector{Symbol})
    list = map(frames) do f
        indices = get_atoms(f, types, atom_types, num_atoms)
        coords = atom_coords[:, indices]
        
        (box_width, coords, max_radius, Δr)
    end
    
    mean(pmap(x -> rdf(x...), list))
end

function rdf_slow{T}(xbw::T, ybw::T, zbw::T, C::Matrix{T}, max_radius::T, Δr::T)
    num_atoms::Int = size(C, 2)
    num_bins::Int = floor(Int, max_radius / Δr)
    int_bins = zeros(UInt64, num_bins)
    bins = Matrix{Float64}(num_bins, 2)

    @inbounds for i in 2:num_atoms
        for j in 1:i-1
            r::T = sqrt(distancesq(xbw, ybw, zbw, C, i, j))
            if r < max_radius
                p::Int = floor(Int, r / Δr) + 1
                int_bins[p] += 1
            end
        end
    end

    @inbounds for i in 1:num_bins
        r_mid = i*Δr - Δr / 2
        divisor = (num_atoms * 4π * r_mid * r_mid * Δr) * (num_atoms / (xbw * ybw * zbw))
        bins[i, 1] = r_mid
        bins[i, 2] = 2 * int_bins[i] / divisor
    end
    
    return bins
end
