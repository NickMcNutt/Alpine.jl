function push!(func::Function, neighbor_atoms::Cluster, n::Int)
    num = neighbor_atoms.num
    new_num = num + n

    if new_num > neighbor_atoms.max_num
        error("Too many atoms in neighborhood! Expand yo buffer fool.")
    end

    neighbor_atoms.num = new_num
    r = (num + 1):new_num
    func(sub(neighbor_atoms.indices, r),
         sub(neighbor_atoms.coords, :, r),
         sub(neighbor_atoms.r, r))
end

function push!(neighbor_atoms::Cluster, n::Int)
    num = neighbor_atoms.num
    new_num = num + n

    if new_num > neighbor_atoms.max_num
        error("Too many atoms in neighborhood! Expand yo buffer fool.")
    end

    neighbor_atoms.num = new_num
    r = (num + 1):new_num
    return sub(neighbor_atoms.indices, r), sub(neighbor_atoms.r, r), sub(neighbor_atoms.coords, :, r)
end

function push!(neighbor_atoms::Cluster, atom_index::Int, coords::Vector{Float64})
    push!(neighbor_atoms, 1) do indices, n_coords, r
        indices[1] = atom_index
        copy!(n_coords, coords)
        r[1] = norm3(coords)
    end
end

function push!(neighbor_atoms::Cluster, atom_indices::AbstractVector{Int}, coords::AbstractMatrix{Float64})
    n = length(atom_indices)

    push!(neighbor_atoms, n) do indices, n_coords, r
        copy!(indices, atom_indices)
        copy!(n_coords, coords)
        map!(i -> norm3(sub(coords, :, i)), r, 1:n)
    end
end

#nonzeroNeighbors(num_all_neighbors::Dict{Symbol, Int}) = filter((k, v) -> v > 0, num_all_neighbors)

#==function allocate_clusters(group::Group, neighbor_types::Dict{Symbol, Int})
    Dict(Pair{Symbol, Dict{AtomID, Cluster}}[t => Dict(Pair{AtomID, Cluster}[id => Cluster(k) for id in group]) for (t, k) in neighbor_types])
end==#

immutable Octree
    box_width::Float64
    cube_width::Float64
    cps::Int64
    NUM_MAX::Int64
    num::Array{Int64, 3}
    indices::Array{Int64, 4}
    coords::Array{Float64, 5}

    function Octree(box_width::Float64, cube_width::Float64, NUM_MAX::Int64 = 8)
        cps = ceil(Int64, box_width / cube_width)
        new(box_width,
            cube_width,
            cps,
            NUM_MAX,
            zeros(Int64, (cps, cps, cps)),
            zeros(Int64, (NUM_MAX, cps, cps, cps)),
            zeros(Float64, (3, NUM_MAX, cps, cps, cps))
        )
    end
end

immutable Parallelepiped{T<:AbstractFloat}
    periodic::NTuple{3, Bool}
    a::NTuple{3, T}
    b::NTuple{3, T}
    c::NTuple{3, T}

    Parallelepiped{T<:AbstractFloat}(periodic::NTuple{3, Bool}, v₁::AbstractVector{T}, v₂::AbstractVector{T}, v₃::AbstractVector{T}) = new(periodic, v₁, v₂, v₃)
end

Parallelepiped{T<:AbstractFloat}(basis::AbstractMatrix{T}) = Parallelepiped(basis[:, 1], basis[:, 2], basis[:, 3])
Parallelepiped{T<:AbstractFloat}(a::T, b::T, c::T, α::T, β::T, γ::T)

immutable Cluster
    num::Int
    max_num::Int
    indices::Indices
    coords::Coords
    r::Vector{Float64}

    Cluster(n::Int) = new(n, n, Indices(n), Vector{Float64}(n), Coords(3, n))
end

function nearestneighbors!{T<:AbstractFloat}(neighbors_i::AbstractVector{Int}, neighbors_r::AbstractVector{T}, origin::AbstractVector{T}, o::Octree, k::Int)
    nx, ny, nz = cube(o, origin)

    r_max = Inf
    fill!(neighbors_r, Inf)
    fill!(neighbors_i, 765)

    num_points = 0
    complete = false

    for s in 0:o.cps
        p = -s:max(s, 1):s
        for dx in p, dy in p, dz in p
            bx, by, bz = cube(o, nx + dx, ny + dy, nz + dz)

            n = o.num[bx, by, bz]
            println("dx: $dx\tdy: $dy\tdz: $dz")
            num_points += n

            for i in 1:n
                println("NP: ", num_points)
                coord = sub(o.coords, :, i, bx, by, bz)
                r = distanceSqPeriodic(origin, coord, o.box_width)

                if 0.000001 < r < r_max
                    j = 1
                    while r > neighbors_r[j]
                        j += 1
                    end

                    for m in k:-1:(j+1)
                        neighbors_r[m] = neighbors_r[m - 1]
                        neighbors_i[m] = neighbors_i[m - 1]
                    end

                    neighbors_r[j] = r
                    neighbors_i[j] = o.indices[i, bx, by, bz]

                    r_max = neighbors_r[k]
                end
            end
        end

        complete && break
        (num_points > k) && (complete = true)
    end

    return neighbors_i, neighbors_r
end

nearestneighbors{T<:AbstractFloat}(origin::AbstractVector{T}, octree::Octree, k::Int) = nearestneighbors!(Vector{Int}(k), Vector{Float64}(k), origin, octree, k)

function kNN!{T<:Float64}(neighbors_i::AbstractVector{Int}, neighbors_r::AbstractVector{T}, neighbors_d::AbstractCoords{T}, origin::AbstractVector{T}, o::Octree, coords::AbstractCoords{T}, k::Int, box_width::Float64)
    nx, ny, nz = cube(o, origin)

    r_max = Inf
    fill!(neighbors_r, Inf)

    num_points = 0
    complete = false

    for s in 0:o.cps
        t = Set()
        for a in (-1, 1), b in (-1, 1), c in (-1, 1), i in 0:s, j in 0:i
            union!(t, Set(collect(permutations([a*s, b*i, c*j]))))
        end

        for u in t
            bx, by, bz = cube(o, nx + u[1], ny + u[2], nz + u[3])

            n = o.num[bx, by, bz]
            num_points += n

            for i in 1:n
                coord = sub(o.coords, :, i, bx, by, bz)
                r = distanceSqPeriodic(origin, coord, box_width)

                if 0.0001 < r < r_max
                    j = 1
                    while r > neighbors_r[j]
                        j += 1
                    end

                    for m in k:-1:(j+1)
                        neighbors_r[m] = neighbors_r[m - 1]
                        neighbors_i[m] = neighbors_i[m - 1]
                    end

                    neighbors_r[j] = r
                    neighbors_i[j] = o.indices[i, bx, by, bz]

                    r_max = neighbors_r[k]
                end
            end
        end

        if complete
            break
        end

        if num_points > k
            complete = true
        end
    end

    distancesPeriodic!(neighbors_d, origin, sub(coords, :, neighbors_i[:]), box_width)

    return
end

@generated function cubefaceindices{N}(::Type{Val{N}})
    indices = Set{Vector{Int}}()
    @inbounds for px in (-1, 1), py in (-1, 1), pz in (-1, 1), d1 in 0:N, d2 in 0:d1
        union!(indices, permutations([px*N, py*d1, pz*d2]))
    end
    return [NTuple{3, Int}((i...)) for i in indices]
end

function kNN!(clusters::Clusters, atom_coords::Matrix{Float64}, o::Octree)
    bw = o.box_width
    num_neighbors = clusters.num_neighbors
    num_clusters = clusters.num_clusters

    ids = clusters.atom_ids_neighbors
    radii = clusters.r²

    @simd for i in eachindex(radii)
        @inbounds radii[i] = Float64(Inf)
    end

    @inbounds for c in 1:num_clusters
        id = clusters.atom_ids_central[c]

        xc = atom_coords[1, id]
        yc = atom_coords[2, id]
        zc = atom_coords[3, id]
        nx, ny, nz = cube(o, xc, yc, zc)

        r_max = Float64(Inf)
        num_points = zero(Int)
        complete = false

        @inbounds for d in 0:o.cps
            @inbounds for u in cubefaceindices(Val{d})
                bx, by, bz = cube(o, nx + u[1], ny + u[2], nz + u[3])

                n = o.num[bx, by, bz]
                num_points += n

                @inbounds for i in 1:n
                    r = distanceSqPeriodic(
                            bw,
                            xc, yc, zc,
                            o.coords[1, i, bx, by, bz], o.coords[2, i, bx, by, bz], o.coords[3, i, bx, by, bz]
                    )

                    if r < r_max && r > 0.0001
                        j = 1
                        @inbounds while r > radii[j, c] j += 1 end

                        @inbounds for m in num_neighbors:-1:(j+1)
                            radii[m, c] = radii[m - 1, c]
                            ids[m, c] = ids[m - 1, c]
                        end

                        radii[j, c] = r
                        ids[j, c] = o.indices[i, bx, by, bz]

                        r_max = radii[num_neighbors, c]
                    end
                end
            end

            complete && break
            num_points > num_neighbors && (complete = true)
        end

        #distancesPeriodic!(sub(clusters.coords, :, :, c), coord_central, sub(atom_coords, :, ids[:, c]), box_width)
    end
end

function kNN!(clusters::Clusters, atom_coords::Matrix{Float64}, box_width::Float64)
    num_atoms::Int = size(atom_coords, 2)
    num_neighbors = clusters.num_neighbors
    num_clusters = clusters.num_clusters
    ac = atom_coords
    ids = clusters.atom_ids_neighbors
    radii = clusters.r²

    @simd for i in eachindex(radii)
        @inbounds radii[i] = Float64(Inf)
    end

     @inbounds for c in 1:num_clusters
        id = clusters.atom_ids_central[c]

        xc = ac[1, id]
        yc = ac[2, id]
        zc = ac[3, id]

        r_max = Float64(Inf)

         @inbounds for i in 1:num_atoms
             r = distanceSqPeriodic(box_width, xc, yc, zc, ac[1, i], ac[2, i], ac[3, i])

             if r < r_max && r > 0.001
                 j = 1
                 while r > radii[j, c] j += 1 end

                 for m in num_neighbors:-1:(j+1)
                     radii[m, c] = radii[m - 1, c]
                     ids[m, c] = ids[m - 1, c]
                 end

                 radii[j, c] = r
                 ids[j, c] = i

                 r_max = radii[num_neighbors, c]
             end
        end
        
        #distancesPeriodic!(sub(clusters.coords, :, :, c), coord_central, sub(ac, :, ids[:, c]), box_width)
    end
end
    
function kNN!{T<:Float64}(neighbors_i::AbstractVector{Int}, neighbors_r::AbstractVector{T}, neighbors_d::AbstractCoords{T}, center::AbstractVector{T}, coords::AbstractCoords{T}, k::Int, box_width::Float64)
    r_max = Inf
    fill!(neighbors_r, Inf)

    @inbounds for i in 1:size(coords, 2)
        coord = sub(coords, :, i)
        r = distanceSqPeriodic(center, coord, box_width)

        if 0.0001 < r < r_max
            j = 1
            while r > neighbors_r[j]
                j += 1
            end

            for m in k:-1:(j+1)
                neighbors_r[m] = neighbors_r[m - 1]
                neighbors_i[m] = neighbors_i[m - 1]
            end

            neighbors_r[j] = r
            neighbors_i[j] = i

            r_max = neighbors_r[k]
        end
    end

    distancesPeriodic!(neighbors_d, center, sub(coords, :, neighbors_i[:]), box_width)

    return nothing
end

function kNN!{T<:AbstractFloat}(cluster::Cluster, center::AbstractVector{T}, coords::AbstractCoords{T}, k::Int, box_width::T)
    r_max = Inf
    fill!(cluster.r, Inf)

    @inbounds for i in 1:size(coords, 2)
        coord = sub(coords, :, i)
        r = distanceSqPeriodic(center, coord, box_width)

        if 0.0001 < r < r_max
            j = 1
            while r > cluster.r[j]
                j += 1
            end

            @inbounds for m in k:-1:(j+1)
                cluster.r[m] = cluster.r[m - 1]
                cluster.indices[m] = cluster.indices[m - 1]
            end

            cluster.r[j] = r
            cluster.indices[j] = i

            r_max = cluster.r[k]
        end
    end

    #distancesPeriodic!(cluster.coords, center, coords[:, cluster.indices], box_width)
    #distancesPeriodic!(cluster.coords, center, sub(coords, :, cluster.indices), box_width)

    return
end

function kNN!{T<:Float64}(distanceFunc, neighbors_i::AbstractVector{Int}, neighbors_r::AbstractVector{T}, center::AbstractVector{T}, coords::AbstractCoords{T}, k::Int)
    r_max = Inf
    fill!(neighbors_r, Inf)

    for i in 1:size(coords, 2)
        coord = sub(coords, :, i)
        #coord ≡ center && continue
        r = distanceFunc(center, coord)

        if r < r_max
            j = 1
            while r > neighbors_r[j]
                j += 1
            end

            for m in k:-1:(j+1)
                neighbors_r[m] = neighbors_r[m - 1]
                neighbors_i[m] = neighbors_i[m - 1]
            end

            neighbors_r[j] = r
            neighbors_i[j] = i

            r_max = neighbors_r[k]
        end
    end
end

function centroid{T<:AbstractFloat}(C::AbstractMatrix{T})
    n = size(C, 2)
    x, y, z = zero(T), zero(T), zero(T)

    @inbounds @simd for i in 1:n
        x += coords[1, i]
        y += coords[2, i]
        z += coords[3, i]
    end

    x /= n
    y /= n
    z /= n

    return x, y, z
end

function centroid{T<:AbstractFloat}(cs::PointSet{T}, bw::T, x::T, y::T, z::T)
    hbw = bw / 2
    hbwx = hbw + x
    nhbwx = -hbw + x
    hbwy = hbw + y
    nhbwy = -hbw + y
    hbwz = hbw + z
    nhbwz = -hbw + z

    n = length(cs)
    C = cs.coords
    sx, sy, sz = zero(T), zero(T), zero(T)
    cx, cy, cz = zero(Int), zero(Int), zero(Int)

    for i in 1:n
        sx += C[1, i]
        sy += C[2, i]
        sz += C[3, i]

        if C[1, i] > hbwx cx -= 1 elseif C[1, i] < nhbwx cx += 1 end
        if C[2, i] > hbwy cy -= 1 elseif C[2, i] < nhbwy cy += 1 end
        if C[3, i] > hbwz cz -= 1 elseif C[3, i] < nhbwz cz += 1 end
    end

    sx = (sx + cx*bw)/n
    sy = (sy + cy*bw)/n
    sz = (sz + cz*bw)/n

    return sx, sy, sz
end
