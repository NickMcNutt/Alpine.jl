# Constructs a Gram matrix and stores it in M
function gramMatrix!(M, coords)
    for i = 1:size(coords, 2)
        for j = 1:size(coords, 2)
            M[i, j] = dot(coords[:, i], coords[:, j])
        end
    end
end

# Computes the d_ref similarity measure between two Gram matrices
function d_ref(coords1, coords2)
    num_coords = size(coords1, 2)

    gram1 = zeros(Float64, (num_coords, num_coords))
    gram2 = zeros(Float64, (num_coords, num_coords))

    gramMatrix!(gram1, coords1)
    indices_order = [1:num_coords]
    for i = 1:factorial(num_coords)
        gramMatrix!(gram2, coords2[:, nthperm(indices_order, i)])
        chol 
        println(norm(gram1 - gram2))
    end
end

# Rotates an array of coordinates to align to the basis axes
# (used to reverse Gram matrix back into coordinates)
function alignWithAxes!(coords, center, box_width)
    # Removes three degrees of freedom from a set of coordinates.
    # Centers the set at (0, 0, 0), and rotates the set
    # so that the vector to coord1 is orthogonal to the y and z axes
    # and the vector to coord2 is orthogonal to the z axis.

    const xaxis = [1.0, 0.0, 0.0]

    transposeCoords!(coords, -center)
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

# Returns a rotation matrix that corresponds to a "θ" degree rotation around vector "axis"
function rodriguesRotationMatrix(axis, θ)
    v = axis / norm(axis)
    W = [ 0 -v[3] v[2] ; v[3] 0 -v[1] ; -v[2] v[1] 0 ]
    I = eye(3)
    I + sind(θ)*W + (2*sind(θ/2)*sind(θ/2))*W*W
end

# Rotates an array of coordinates using a Rodrigues Rotation Matrix.  Modifies "coords".
function rotateCoords!(coords, axis, θ)
    M = rodriguesRotationMatrix(axis, θ)
    for i = 1:size(coords, 2)
        coords[:, i] = M*coords[:, i]
    end
end

# Rotates an array of coordinates in the same way that "vector1" rotates to "vector2". Modifies "coords".
function rotateCoordsVector!(coords, vector1, vector2)
    axis = cross(vector1, vector2)
    θ = -acosd(dot(vector1, vector2) / (norm(vector1) * norm(vector2)))
    rotateCoords!(coords, axis, θ)
end

# Shifts coordinates by a vector. Modifies "coords".
function transposeCoords!(coords, vector)
    for i = 1:size(coords, 2)
        @devec coords[:, i] += vector[:]
    end
end

# Wrap an atom coordinate that spans across a periodic boundary condition
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

# Wraps a single atom coordinate (only one frame)
function wrapCoord!(atom, index, box_width)
    wrapCoord!(atom, index, 1, box_width)
end

# Wraps many atom coordinates
function wrapCoords!(atoms, num_frame, box_width)
    for i = 1:size(atoms, 2)
        wrapCoord!(atoms, i, 1, box_width)
    end
end

# Wraps many atom coordinates (only one frame)
function wrapCoords!(atoms, box_width)
    wrapCoords!(atoms, 1, box_width)
end

# Computes distance squared between two atoms for a periodic box width
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
