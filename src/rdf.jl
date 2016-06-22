using VPTrees

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

function rdf{T}(xbw::T, ybw::T, zbw::T, coords::Matrix{T}, radius::T, Δr::T)
    xhbw = xbw / 2
    yhbw = ybw / 2
    zhbw = zbw / 2

    function metric(i::Int, j::Int)
        @inbounds x = coords[1, i] - coords[1, j]
        @inbounds y = coords[2, i] - coords[2, j]
        @inbounds z = coords[3, i] - coords[3, j]

        if x < -xhbw x += xbw elseif x > xhbw x -= xbw end
        if y < -yhbw y += ybw elseif y > yhbw y -= ybw end
        if z < -zhbw z += zbw elseif z > zhbw z -= zbw end

        return sqrt(x^2 + y^2 + z^2)
    end

    num_bins::Int = ceil(Int, radius / Δr)
    bins = Matrix{T}(num_bins, 2)

    num_atoms = size(coords, 2)
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

function rdf{N, T}(frames::AbstractArray{Int}, types::NTuple{N, Symbol}, bw::Vector{Tuple{T, T, T}}, max_radius::T, Δr::T, num_atoms::Int, atom_coords::AbstractMatrix{T}, atom_types::Vector{Symbol})
    list = map(frames) do f
        indices = get_atoms(f, types, atom_types, num_atoms)
        coords = atom_coords[:, indices]
        
        (bw[f][1], bw[f][2], bw[f][3], coords, max_radius, Δr)
    end
    
    mean(pmap(x -> rdf(x...), list))
end
