function rdf{T}(bw::T, C::Matrix{T}, max_radius::T, dr::T)
    num_atoms::Int = size(C, 2)
    num_bins::Int = floor(Int, max_radius / dr)
    int_bins = zeros(UInt64, num_bins)
    bins = Matrix{Float64}(num_bins, 2)

    @inbounds for i in 2:num_atoms
        for j in 1:i-1
            r::T = sqrt(distancesq(bw, C, i, j))
            if r < max_radius
                p::Int = ceil(Int, r / dr) - 1
                int_bins[p] += 1
            end
        end
    end

    @inbounds for i in 1:num_bins
        r_mid = i*dr - dr / 2
        divisor = (num_atoms * 4π * r_mid * r_mid * dr) * (num_atoms / bw ^ 3)
        bins[i, 1] = r_mid
        bins[i, 2] = 2 * int_bins[i] / divisor
    end
    
    return bins
end

function rdf{N, T}(frames::AbstractArray{Int}, types::NTuple{N, Symbol}, box_width::T, max_radius::T, dr::T, num_atoms::Int, atom_coords::AbstractMatrix{T}, atom_types::Vector{Symbol})
    list = map(frames) do f
        indices = get_atoms(f, types, atom_types, num_atoms)
        coords = atom_coords[:, indices]
        
        (box_width, coords, max_radius, dr)
    end
    
    mean(pmap(x -> rdf(x...), list))
end

function rdf{T}(xbw::T, ybw::T, zbw::T, C::Matrix{T}, max_radius::T, dr::T)
    num_atoms::Int = size(C, 2)
    num_bins::Int = floor(Int, max_radius / dr)
    int_bins = zeros(UInt64, num_bins)
    bins = Matrix{Float64}(num_bins, 2)

    @inbounds for i in 2:num_atoms
        for j in 1:i-1
            r::T = sqrt(distancesq(xbw, ybw, zbw, C, i, j))
            if r < max_radius
                p::Int = ceil(Int, r / dr) - 1
                int_bins[p] += 1
            end
        end
    end

    @inbounds for i in 1:num_bins
        r_mid = i*dr - dr / 2
        divisor = (num_atoms * 4π * r_mid * r_mid * dr) * (num_atoms / (xbw * ybw * zbw))
        bins[i, 1] = r_mid
        bins[i, 2] = 2 * int_bins[i] / divisor
    end
    
    return bins
end


function rdf{N, T}(frames::AbstractArray{Int}, types::NTuple{N, Symbol}, bw::Vector{Tuple{T, T, T}}, max_radius::T, dr::T, num_atoms::Int, atom_coords::AbstractMatrix{T}, atom_types::Vector{Symbol})
    list = map(frames) do f
        indices = get_atoms(f, types, atom_types, num_atoms)
        coords = atom_coords[:, indices]
        
        (bw[f][1], bw[f][2], bw[f][3], coords, max_radius, dr)
    end
    
    mean(pmap(x -> rdf(x...), list))
end
