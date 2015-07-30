include("types.jl")

import Base.LinAlg.BLAS: gemm!

@generated function normN{N, T<:Real}(::Type{Val{N}}, v::AbstractVector{T})
    x = :(0)
    for i in 1:N
        x = :(v[$i] * v[$i] + $x)
    end
    x = :(sqrt($x))
    return x
end

norm3{T<:Real}(v::AbstractVector{T}) = normN(Val{3}, v)

normFrobenius{T<:FloatingPoint}(M::AbstractCoords{T}) = vecnorm(M)

function rmsd{T<:FloatingPoint}(M::AbstractCoords{T})
    total = zero(eltype(M))
    @inbounds @simd for i = 1:length(M)
        total += M[i] ^ 2
    end

    return sqrt(total)
end
    

function rmsd{T<:FloatingPoint}(M1::AbstractCoords{T}, M2::AbstractCoords{T})
    total = zero(eltype(M1))
    @inbounds @simd for i = 1:length(M1)
        total += (M1[i] - M2[i])^2
    end

    return sqrt(total)
end

function kNN{T<:Float64}(distanceFunc, center::AbstractVector{T}, coords::AbstractCoords{T}, k::Int)
    neighbors_i = zeros(Int, k)
    neighbors_r = zeros(Float64, k)
    kNN!(distanceFunc, neighbors_i, neighbors_r, center, coords, k)

    return neighbors_i, neighbors_r
end
    
function kNN!{T<:Float64}(neighbors_i::AbstractVector{Int}, neighbors_r::AbstractVector{T}, neighbors_d::AbstractCoords{T}, center::AbstractVector{T}, coords::AbstractCoords{T}, k::Int, box_width::Float64)
    r_max = Inf
    fill!(neighbors_r, Inf)

    for i in 1:size(coords, 2)
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

function centerOfMass{T<:FloatingPoint}(com::AbstractVector{T}, coords::AbstractCoords{T})
    sum = zeros(Float64, 3)

    @inbounds @simd for i in size(com, 2)
        x1 += coords[1, i]
        x2 += coords[2, i]
        x3 += coords[3, i]
    end

    x1 /= size(com, 2)
    x2 /= size(com, 2)
    x3 /= size(com, 2)
end

function rodriguesRotationMatrix{T<:Real}(axis::AbstractVector{T}, θ::T)
    v = axis / norm(axis)
    W = Float64[ 0 -v[3] v[2] ; v[3] 0 -v[1] ; -v[2] v[1] 0 ]
    I = eye(3)
    I + sind(θ)*W + (2*sind(θ/2)*sind(θ/2))*W*W
end

function rotateCoords!{T<:Real}(coords::AbstractCoords{T}, axis::AbstractVector{T}, θ::T)
    M = rodriguesRotationMatrix(axis, θ)
    gemm!('N', 'N', one(Float64), M, coords, zero(Float64), coords)
end

function rotateCoordsVector!(coords, vector1, vector2)
    axis = cross(vector1, vector2)
    θ = -acosd(dot(vector1, vector2) / (norm(vector1) * norm(vector2)))
    rotateCoords!(coords, axis, θ)
end

function unwrapCoords!(box_width::Float64, initial_coords, coords)
    for i = 1:size(coords, 2)
        for j = 1:3
            if coords[j, i] - initial_coords[j, i] > box_width/2.0
                coords[j, i] -= box_width
            elseif coords[j, i] - initial_coords[j, i] < -box_width/2.0
                coords[j, i] += box_width
            end
        end
    end
end

unwrapCoords!(i::Int, args...) = unwrapCoords!(Float64(i), args...)

function wrapCoord!(atom, index, num_frame, box_width)
    for d = 1:3
        x = atom[d, index, num_frame]
        if x > box_width/2
            x -= box_width
        elseif x < -box_width/2
            x += box_width
        end
        atom[d, index, num_frame] = x
    end
end

function wrapCoord!(atom, index, box_width)
    wrapCoord!(atom, index, 1, box_width)
end

function wrapCoords!(atoms, num_frame, box_width)
    for i = 1:size(atoms, 2)
        wrapCoord!(atoms, i, 1, box_width)
    end
end

function wrapCoords!(atoms, box_width)
    wrapCoords!(atoms, 1, box_width)
end

function alignWithAxes!(coords, center, box_width)
    # Removes three degrees of freedom from a set of coordinates.
    # Centers the set at (0, 0, 0), and rotates the set
    # so that the vector to coord1 is orthogonal to the y and z axes
    # and the vector to coord2 is orthogonal to the z axis.

    const xaxis = [1.0, 0.0, 0.0]

    transposeCoords!(coords, -center)
    #should be moved upa  level
    wrapCoords!(coords, box_width)

    coord1 = coords[:, 1]
    coord2 = coords[:, 2]

    rotateCoordsVector!(coord2, xaxis, coord1)
    θ = acosd(coord2[2] ./ sqrt(coord2[2]^2 + coord2[3]^2))
    rotateCoords!(coord2, xaxis, θ)
    if abs(coord2[3]) > 0.00001
        θ = -θ
    end

    rotateCoordsVector!(coords, xaxis, coord1)
    rotateCoords!(coords, xaxis, θ)
end

function randomRotate!(coords, center, box_width)

    transposeCoords!(coords, -center)
    wrapCoords!(coords, box_width)


    rand_vec = randn(3)
    rand_vec /= norm(rand_vec)

    rotateCoords!(coords, rand_vec, 2π*rand())
end

function minMSD(array_coords)
    
end

function rdf(atom_coords::Array{Float64, 3}, box_width::Float64, max_radius::Float64, dr::Float64)
    num_frames::Int = size(atom_coords, 3)
    num_atoms::Int = size(atom_coords, 2)
    num_bins::Int = floor(Int, max_radius / dr)
    int_bins = zeros(Uint64, num_bins)
    bins = Array{Float64}(num_bins, 2)

    for f = 1:num_frames
        for i = 2:num_atoms
            for j = 1:(i-1)
                r = sqrt(distanceSqPeriodic(atom_coords, i, atom_coords, j, f, box_width))
                if r < max_radius
                    int_bins[floor(Int, r / dr) + 1] += 1
                end
            end
        end
    end

    for i = 1:num_bins
        r_mid = i*dr - dr/2.0
        divisor = (num_atoms * 4.0π * r_mid * r_mid * dr) * (num_atoms / (box_width*box_width*box_width))
        bins[i, 1] = r_mid
        bins[i, 2] = 2.0 * int_bins[i] / divisor / num_frames
    end
    
    return bins
end

macro distanceSqPeriodic(box_width::Real)
    N = convert(Float64, box_width)
    M = N / 2.0
    ex = :(s = $(zero(Float64)))
    for i = 1:3
        x = symbol(:x, i)
        ex = :(begin
            $ex
            $x = c1[$i] - c2[$i]
            if $x > $M
                $x -= $N
            elseif $x < -$M
                $x += $N
            end
            s += $x * $x
        end)
    end

    ex = :(begin
        function distanceSqPeriodic{T<:FloatingPoint}(c1::AbstractVector{T}, c2::AbstractVector{T})
            @inbounds begin
                $ex
                return s
            end
        end
    end)

    return ex
end

distanceSqPeriodic(box_width::Real) = eval(:(@distanceSqPeriodic $box_width))

function distanceSqPeriodic{T<:Real}(c1::AbstractVector{T}, c2::AbstractVector{T}, box_width::Float64)
    @inbounds begin
        x = c1[1] - c2[1]
        if x > box_width/2
            x -= box_width
        elseif x < -box_width/2
            x += box_width
        end

        y = c1[2] - c2[2]
        if y > box_width/2
            y -= box_width
        elseif y < -box_width/2
            y += box_width
        end

        z = c1[3] - c2[3]
        if z > box_width/2
            z -= box_width
        elseif z < -box_width/2
            z += box_width
        end

        return x*x + y*y + z*z
    end
end


distancePeriodic{T<:Real}(c1::AbstractVector{T}, c2::AbstractVector{T}, box_width::Float64) = sqrt(distanceSqPeriodic(c1, c2, box_width))

function distancesSqPeriodic!{T<:Real}(d::AbstractVector{T}, c1::AbstractVector{T}, c2::AbstractVector{T}, box_width::Float64)
    @inbounds begin
        x = c1[1] - c2[1]
        if x > box_width/2
            x -= box_width
        elseif x < -box_width/2
            x += box_width
        end
        d[1] = x

        y = c1[2] - c2[2]
        if y > box_width/2
            y -= box_width
        elseif y < -box_width/2
            y += box_width
        end
        d[2] = y

        z = c1[3] - c2[3]
        if z > box_width/2
            z -= box_width
        elseif z < -box_width/2
            z += box_width
        end
        d[3] = z

        return x*x + y*y + z*z
    end
end

distancesPeriodic!{T<:Real}(d::AbstractVector{T}, c1::AbstractVector{T}, c2::AbstractVector{T}, box_width::Float64) = sqrt(distancesSqPeriodic!(d, c1, c2, box_width))

function distancesPeriodic!{T<:Real}(r::AbstractVector{T}, d::AbstractMatrix{T}, c1::AbstractVector{T}, c2::AbstractMatrix{T}, box_width::Float64)
    for i in 1:length(r)
        r[i] = distancesPeriodic!(sub(d, :, i), c1, sub(c2, :, i), box_width)
    end
end

function distancesPeriodic!{T<:Real}(d::AbstractMatrix{T}, c1::AbstractVector{T}, c2::AbstractMatrix{T}, box_width::Float64)
    const bw = box_width
    const hbw = box_width/2
    @inbounds @simd for i in 1:size(d, 2)
        x = c1[1] - c2[1, i]
        if x > hbw
            x -= bw
        elseif x < -hbw
            x += bw
        end
        d[1, i] = x

        y = c1[2] - c2[2, i]
        if y > hbw
            y -= bw
        elseif y < -hbw
            y += bw
        end
        d[2, i] = y

        z = c1[3] - c2[3, i]
        if z > hbw
            z -= bw
        elseif z < -hbw
            z += bw
        end
        d[3, i] = z
    end
end
