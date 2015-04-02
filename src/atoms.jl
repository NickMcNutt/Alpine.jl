function whichEnergyGroup(energy_groups, atoms, index, num_frame)
    energy_group = 0
    
    for i = 1:length(energy_groups)
        @inbounds if energy_groups[i][1] <= atoms[4, index, num_frame] < energy_groups[i][2]
            energy_group = i
            break
        end
    end

    return energy_group
end

function atomIndicesToCoords(atoms, indices::Array{Int, 1}, num_frame)
    atoms[1:3, indices, num_frame]
end

function atomIndexToCoords(atoms, index::Int, num_frame)
    atoms[1:3, index, num_frame]
end

function atomIndicesToCoords!(new_atoms, atoms, indices::Array{Int, 1}, num_frame)
    for i = 1:length(indices)
        @simd for d = 1:3
            new_atoms[d, i] = atoms[d, indices[i], num_frame]
        end
    end
end

function atomIndexToCoords!(new_atom, atoms, index::Int, num_frame)
    @simd for d = 1:3
        new_atom[d] = atoms[d, index, num_frame]
    end
end
