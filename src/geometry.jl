function centerOfMass(atoms, atom_indices)
    sum = zeros(Float64, 3)
    for i in atom_indices
        sum += atoms[:, i]
    end
    sum /= length(atom_indices)
end

function rodriguesRotationMatrix(axis, θ)
    v = axis / norm(axis)
    W = [ 0 -v[3] v[2] ; v[3] 0 -v[1] ; -v[2] v[1] 0 ]
    I = eye(3)
    I + sind(θ)*W + (2*sind(θ/2)*sind(θ/2))*W*W
end

function rotateCoords!(coords, axis, θ)
    M = rodriguesRotationMatrix(axis, θ)
    for i = 1:size(coords, 2)
        coords[:, i] = M*coords[:, i]
    end
end

function rotateCoordsVector!(coords, vector1, vector2)
    axis = cross(vector1, vector2)
    θ = -acosd(dot(vector1, vector2) / (norm(vector1) * norm(vector2)))
    rotateCoords!(coords, axis, θ)
end

function transposeCoords!(coords, vector)
    for i = 1:size(coords, 2)
        @devec coords[:, i] += vector[:]
    end
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

unwrapAtoms!(i::Int, args...) = unwrapAtoms!(Float64(i), args...)

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

function distanceSqPeriodic(atom1, i1, atom2, i2, num_frame, box_width)
    @inbounds x = atom2[1, i2, num_frame] - atom1[1, i1, num_frame]
    if x > box_width/2
        x -= box_width
    elseif x < -box_width/2
        x += box_width
    end

    @inbounds y = atom2[2, i2, num_frame] - atom1[2, i1, num_frame]
    if y > box_width/2
        y -= box_width
    elseif y < -box_width/2
        y += box_width
    end

    @inbounds z = atom2[3, i2, num_frame] - atom1[3, i1, num_frame]
    if z > box_width/2
        z -= box_width
    elseif z < -box_width/2
        z += box_width
    end

    return x*x + y*y + z*z
end

function weylMatrix!(M, coords)
    for i = 1:size(coords, 2)
        for j = 1:size(coords, 2)
            M[i, j] = dot(coords[:, i], coords[:, j])
        end
    end
end


function distanceSqPeriodic(atom1, i1, atom2, i2, box_width)
    distanceSqPeriodic(atom1, i1, atom2, i2, 1, box_width)
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
