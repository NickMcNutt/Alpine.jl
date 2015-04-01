function readXYZFile(filename)
    file = open(filename, "r")

    num_atoms = parse(Int, readline(file))
    readline(file)
    lines = readlines(file)

    parts = split(lines[1], r" +|\t")
    num_parts = length(parts)
    num_misc_fields = num_parts - 4
    comment_present = false

    if (num_parts > 4) && (parts[5] == "#")
        num_misc_fields -= 1
        comment_present = true
    end

    if num_misc_fields < 0
        num_misc_fields = 0
    end

    atom_types = Array(String, num_atoms)
    atom_coords = Array(Float64, (3, num_atoms))
    atom_misc = Array(Any, (num_misc_fields, num_atoms))

    for i = 1:length(lines)
        line = chomp(lines[i])
        parts = split(line, r" +|\t")

        atom_types[i] = string(parts[1])
        atom_coords[1:3, i] = map(x -> parse(Float64, x), parts[2:4])
        if num_misc_fields > 0
            if comment_present
                atom_misc[:, i] = parts[6:(6 + num_misc_fields - 1)]
            else
                atom_misc[:, i] = parts[5:(5 + num_misc_fields - 1)]
            end
        end
    end

    return (size(atom_coords, 2), atom_types, atom_coords, atom_misc)
end

function countFrames!(file_xyz)
    seekstart(file_xyz)
    num_atoms_total = int(readline(file_xyz))
    seekstart(file_xyz)
    num_lines = countlines(file_xyz)

    return int(num_lines / (num_atoms_total + 2))
end

function countFramesTRJ!(file_trj)
    seekstart(file_trj)
    for i = 1:3
        readline(file_trj)
    end

    num_atoms_total = int(readline(file_trj))
    seekstart(file_trj)
    num_lines = countlines(file_trj)

    return int(num_lines / (num_atoms_total + 9))
end

function countAtoms!(file, atom_type)
    seekstart(file)
    num_atoms_total = int(readline(file))
    readline(file)

    num_atoms_type = 0
    for i = 1:num_atoms_total
        line = readline(file)
        if lowercase(split(chomp(line), ' ')[1]) == atom_type
            num_atoms_type += 1
        end
    end

    return num_atoms_type
end

function countAtomsTRJ!(file)
    seekstart(file)
    for i = 1:3
        readline(file)
    end

    num_atoms_total = int(readline(file))

    return num_atoms_total
end

function countAtomsTRJ!(file, atom_type)
    seekstart(file)
    for i = 1:3
        readline(file)
    end

    num_atoms_total = int(readline(file))
    
    for i = 1:5
        readline(file)
    end

    num_atoms_type = 0
    for i = 1:num_atoms_total
        line = readline(file)
        if lowercase(split(chomp(line), ' ')[1]) == atom_type
            num_atoms_type += 1
        end
    end

    return num_atoms_type
end

function createReadersTRJ(file, atom_type)
    local atoms::Array{Float64, 3}
    #local atom_ids_to_delete::Array{Bool, 1}

    function createData!(file, num_frames)
        num_atoms_type = countAtomsTRJ!(file, atom_type)
        atoms = Array(Float64, (4, num_atoms_type, num_frames))
        #num_atoms_total = countAtomsTRJ!(file)
        #atom_ids_to_delete = Array(Bool, num_atoms_total)
        #for i = 1:num_atoms_total
        #    atom_ids_to_delete[i] = false
        #end

        return atoms
    end

    function skipFrame!(file)
        for i = 1:3
            readline(file)
        end

        num_atoms_total = int(readline(file))

        for i = 1:(num_atoms_total + 5)
            readline(file)
        end
    end

    function readFrame!(file, current_frame)
        for i = 1:3
            readline(file)
        end

        num_atoms_total = int(readline(file))

        for i = 1:5
            readline(file)
        end
        
        current_atom = 1
        for i = 1:num_atoms_total
            line = readline(file)
            atom = split(chomp(line), ' ')
            if lowercase(atom[1]) == atom_type
                atoms[1:3, current_atom, current_frame] = float(atom[2:4])
                atoms[4, current_atom, current_frame] = float(atom[6])
                current_atom += 1
            end
        end
    end

    return (createData!, skipFrame!, readFrame!)
end

function createReadersXYZ(file, atom_type)
    local atoms::Array{Float64, 3}

    function createData!(file, num_frames)
        num_atoms_type = countAtoms!(file, atom_type)
        atoms = Array(Float64, (3, num_atoms_type, num_frames))
    end

    function skipFrame!(file)
        num_atoms_total = int(readline(file))
        for i = 1:(num_atoms_total + 1)
            readline(file)
        end
    end

    function readFrame!(file, current_frame)
        num_atoms_total = int(readline(file))
        readline(file)
        
        current_atom = 1
        for i = 1:num_atoms_total
            line = readline(file)
            atom = split(chomp(line), ' ')
            if lowercase(atom[1]) == atom_type
                atoms[1:3, current_atom, current_frame] = float(atom[2:4])
                current_atom += 1
            end
        end
    end

    return (createData!, skipFrame!, readFrame!)
end

function createReadersEnergies(file, num_atoms_type)
    local energies::Array{Float64, 2}

    function createData!(file, num_frames)
        energies = Array(Float64, (num_atoms_type, num_frames))
    end

    function skipFrame!(file)
        for i = 1:3
            readline(file)
        end

        num_atoms = int(readline(file))

        for i = 1:5+num_atoms
            readline(file)
        end
    end

    function readFrame!(file, current_frame)
        for i = 1:3
            readline(file)
        end

        num_atoms = int(readline(file))

        for i = 1:5
            readline(file)
        end

        for i = 1:num_atoms
            energies[i, current_frame] = float(chomp(readline(file)))
        end
    end

    return (createData!, skipFrame!, readFrame!)
end

function readFramesFile!(file, start_frame_num, num_frames, createData!, skipFrame!, readFrame!)
    data = createData!(file, num_frames)

    seekstart(file)
    for current_frame = 1:(start_frame_num - 1)
        skipFrame!(file)
        print("Skipping frame: $(current_frame)\r")
    end

    for current_frame = 1:num_frames
        readFrame!(file, current_frame)
        print("Loading frame: $(start_frame_num + current_frame - 1) \r")
    end
    println()

    return data
end

function writeTRJFile(filename, frames::Array{Int, 1}, box_width, atom_names::Array{String, 1}, atoms)
    file = open(filename, w)
        for f in frames
            println(file, "ITEM: TIMESTEP")
            println(file, f)
            println(file, "ITEM: NUMBER OF ATOMS")
            println(file, size(atoms, 2))
            println(file, "ITEM: BOX BOUNDS pp pp pp")
            println(file, "-$(box_width/2) $(box_width/2)")
            println(file, "-$(box_width/2) $(box_width/2)")
            println(file, "-$(box_width/2) $(box_width/2)")
            println(file, "ITEM: ATOMS element xu yu zu")
            for i = 1:size(atoms, 2)
                @printf(file, "C %.4f %.4f %.4f\n", atoms[1, i, f], atoms[2, i, f], atoms[3, i, f])
            end
        end
    close(file)
end

function printTRJ(file, f, box_width, atoms)
    println(file, "ITEM: TIMESTEP")
    println(file, f)
    println(file, "ITEM: NUMBER OF ATOMS")
    println(file, size(atoms, 2) + 1)
    println(file, "ITEM: BOX BOUNDS pp pp pp")
    println(file, "-$(box_width/2) $(box_width/2)")
    println(file, "-$(box_width/2) $(box_width/2)")
    println(file, "-$(box_width/2) $(box_width/2)")
    println(file, "ITEM: ATOMS element xu yu zu")
    println(file, "Li 0.0 0.0 0.0")
    for i = 1:size(atoms, 2)
        @printf(file, "C %.4f %.4f %.4f\n", atoms[1, i, 1], atoms[2, i, 1], atoms[3, i, 1])
    end
end
