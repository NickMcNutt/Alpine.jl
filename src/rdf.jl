using VPTrees

# The distance metric for points in a periodic box
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

function distance{T}(coords::Matrix{T}, coords2::Matrix{T}, xbw::T, ybw::T, zbw::T)
    xhbw = xbw / 2
    yhbw = ybw / 2
    zhbw = zbw / 2

    function(i::Int, j::Int)
        @inbounds x = coords[1, i] - coords2[1, j]
        @inbounds y = coords[2, i] - coords2[2, j]
        @inbounds z = coords[3, i] - coords2[3, j]

        if x < -xhbw x += xbw elseif x > xhbw x -= xbw end
        if y < -yhbw y += ybw elseif y > yhbw y -= ybw end
        if z < -zhbw z += zbw elseif z > zhbw z -= zbw end

        return sqrt(x^2 + y^2 + z^2)
    end
end

# NDF for a single set of points
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

# NDF for component functions of a set of points
function ndf{T}(f::Function, metric::Function, point_indices::Vector{Int}, radius::T, Δr::T)
    num_points::Int = length(point_indices)
    tree = VPTree(metric, point_indices)
    neighbors = NeighborList{T}(num_points)

    @inbounds for i in 1:num_points
        rangesearch!(neighbors, tree, i, radius)

        k = length(neighbors)::Int
        for j in 1:k
            ni = neighbors.indices[j]::Int
            i == ni && continue
            r::T = neighbors.distances[j]::T
            b::Int = floor(Int, r / Δr)::Int + 1
            f(b, i, ni)
        end
    end
end

# NDF for component functions of one set of points against another set
function ndf{T}(f::Function, metric::Function, metric2::Function, point_indices::Vector{Int}, point_indices_2::Vector{Int}, radius::T, Δr::T)
    num_points::Int = length(point_indices)
    tree = VPTree(metric, point_indices)
    neighbors = NeighborList{T}(num_points)

    num_points_2 = length(point_indices_2)
    @inbounds for i in 1:num_points_2
        rangesearch!(metric2, neighbors, tree, i, radius)

        k = length(neighbors)::Int
        for j in 1:k
            ni = neighbors.indices[j]::Int
            r::T = neighbors.distances[j]::T
            b::Int = floor(Int, r / Δr)::Int + 1
            f(b, i, ni)
        end
    end
end

# PDF for a collection of atoms in a periodic box of size (xbw, ybw, zbw)
function rdf{T}(atoms::Atoms, xbw::T, ybw::T, zbw::T, radius::T, Δr::T)
    num_bins::Int = ceil(Int, radius / Δr)
    bins = Matrix{T}(num_bins, 2)

    metric = distance(collect(atoms[:coords]), xbw, ybw, zbw)
    num_atoms = length(atoms)::Int
    point_indices = collect(1:num_atoms)::Vector{Int}
    int_bins = ndf(metric, point_indices, radius, Δr)

    for i in 1:num_bins
        r_mid = i*Δr - Δr/2
        divisor = (num_atoms * 4π * r_mid * r_mid * Δr) * (num_atoms / (xbw * ybw * zbw))
        bins[i, 1] = r_mid
        bins[i, 2] = int_bins[i] / divisor
    end

    return bins
end

# PDF for component functions for a collection of atoms in a periodic box of size (xbw, ybw, zbw)
function rdf{T}(f::Function, nf::Int, atoms::Atoms, xbw::T, ybw::T, zbw::T, radius::T, Δr::T)
    num_bins::Int = ceil(Int, radius / Δr)
    bins = Matrix{T}(nf, num_bins)
    int_bins = zeros(UInt64, nf, num_bins)

    metric = distance(collect(atoms[:coords]), xbw, ybw, zbw)
    num_atoms = length(atoms)::Int
    point_indices = collect(1:num_atoms)::Vector{Int}
    ndf((b::Int, i::Int, j::Int) -> f(atoms, int_bins, b, i, j), metric, point_indices, radius, Δr)

    for i in 1:num_bins
        r_mid = i*Δr - Δr/2
        divisor = (num_atoms * 4π * r_mid^2 * Δr) * (num_atoms / (xbw * ybw * zbw))
        for g in 1:nf
            bins[g, i] = int_bins[g, i] / divisor
        end
    end

    return bins
end

# PDF for component functions of one set of points against another set
function rdf{T}(f::Function, nf::Int, atoms::Atoms, atoms2::Atoms, xbw::T, ybw::T, zbw::T, radius::T, Δr::T)
    num_bins::Int = ceil(Int, radius / Δr)
    bins = Matrix{T}(nf, num_bins)
    int_bins = zeros(UInt64, nf, num_bins)

    metric = distance(collect(atoms[:coords]), xbw, ybw, zbw)
    metric2 = distance(collect(atoms[:coords]), collect(atoms2[:coords]), xbw, ybw, zbw)

    num_atoms = length(atoms)::Int
    num_atoms2 = length(atoms2)::Int

    point_indices = collect(1:num_atoms)::Vector{Int}
    point_indices2 = collect(1:num_atoms2)::Vector{Int}

    ndf((b::Int, i::Int, j::Int) -> f(atoms, atoms2, int_bins, b, i, j), metric, metric2, point_indices, point_indices2, radius, Δr)

    for i in 1:num_bins
        r_mid = i*Δr - Δr/2
        divisor = (4π * r_mid^2 * Δr) * (num_atoms * num_atoms2) / (xbw * ybw * zbw)
        for g in 1:nf
            bins[g, i] = int_bins[g, i] / divisor
        end
    end

    return bins
end

# PDF for component functions for a collection of atoms in a frame
@inline function rdf{T}(f::Function, nf::Int, frame::Frame, radius::T, Δr::T)
    rdf(f, nf, frame[:atoms], box_dims(frame)..., radius, Δr)
end

# PDF for component functions of one set of points against another set.
# The frame supplies the box dimensions
@inline function rdf{T}(f::Function, nf::Int, frame::Frame, atoms1::Atoms, atoms2::Atoms, radius::T, Δr::T)
    rdf(f, nf, atoms1, atoms2, box_dims(frame)..., radius, Δr)
end

# Frame-averaged PDF for component functions for a collection of frames
function rdf{T}(f::Function, nf::Int, frames::Vector{Frame}, radius::T, Δr::T)
    bins = pmap(frames) do frame
        rdf(f, nf, frame, radius, Δr)::Matrix{T}
    end

    mean(bins)
end

# Frame-averaged PDF for component functions of one set of points against another set for a collection of frames
function rdf{T}(f::Function, nf::Int, frames::Vector{Frame}, atoms1::Vector{Atoms}, atoms2::Vector{Atoms}, radius::T, Δr::T)
    num_frames = length(frames)
    list = [(frames[i], atoms1[i], atoms2[i]) for i in 1:num_frames]

    bins = pmap(list) do item
        rdf(f, nf, item[1], item[2], item[3], radius, Δr)::Matrix{T}
    end

    mean(bins)
end

# Returns the x-axis for a PDF
function rdf_axis{T}(radius::T, Δr::T)
    num_bins = ceil(Int, radius / Δr)::Int
    T[i*Δr - Δr/2 for i in 1:num_bins]
end

# A slow version of the PDF for a collection of coords in a periodic box of size (xbw, ybw, zbw)
function rdf_slow{T}(xbw::T, ybw::T, zbw::T, coords::Matrix{T}, max_radius::T, Δr::T)
    num_atoms::Int = size(coords, 2)
    num_bins::Int = floor(Int, max_radius / Δr)
    int_bins = zeros(UInt64, num_bins)
    bins = Matrix{Float64}(num_bins, 2)

    @inbounds for i in 2:num_atoms
        for j in 1:i-1
            r::T = sqrt(distancesq(xbw, ybw, zbw, coords, i, j))
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
