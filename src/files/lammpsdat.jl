type LAMMPSDAT <: AtomFileType end

typealias DAT LAMMPSDAT

file_exts(::Type{LAMMPSDAT}) = (".dat",)

function Atoms(::Type{LAMMPSDAT}, num_atoms::Int)
    Atoms(num_atoms,
        :molecule => Vector{Int}(num_atoms),
        :type => Vector{Symbol}(num_atoms),
        :charge => Vector{Float64}(num_atoms),
        :coords => Matrix{Float64}(3, num_atoms)
    )
end

function read_frame_header{T}(::Type{T}, ::Type{LAMMPSDAT}, chunk::Chunk)
    num_atoms_regex = r"(\d+)\s*atoms"
    num_bonds_regex = r"(\d+)\s*bonds"
    x_dims_regex = r"(\S+)\s+(\S+)\s+xlo\s+xhi"
    y_dims_regex = r"(\S+)\s+(\S+)\s+ylo\s+yhi"
    z_dims_regex = r"(\S+)\s+(\S+)\s+zlo\s+zhi"
    masses_regex = r"Masses"
    atoms_regex = r"Atoms"
    bonds_regex = r"Bonds"

    frame_start = chunk.i
    
    num_atoms = parse(Int, read_through_regex(chunk, num_atoms_regex)[1])
    num_bonds = parse(Int, read_through_regex(chunk, num_bonds_regex)[1])

    x_dims = read_through_regex(chunk, x_dims_regex)
    x_lo = parse(Float64, x_dims[1])
    x_hi = parse(Float64, x_dims[2])

    y_dims = read_through_regex(chunk, y_dims_regex)
    y_lo = parse(Float64, y_dims[1])
    y_hi = parse(Float64, y_dims[2])

    z_dims = read_through_regex(chunk, z_dims_regex)
    z_lo = parse(Float64, z_dims[1])
    z_hi = parse(Float64, z_dims[2])

    read_through_regex(chunk, atoms_regex)
    read_newline(chunk)

    atom_data_start = chunk.i
    read_newlines(chunk, num_atoms)
    atom_data_size = chunk.i - atom_data_start

    read_through_regex(chunk, bonds_regex)
    read_newline(chunk)

    bond_data_start = chunk.i
    read_newlines(chunk, num_bonds)
    bond_data_size = chunk.i - bond_data_start

    frame_end = endof(chunk.data)
    frame_size = frame_end - frame_start + 1
    chunk.i = frame_end
    
    frame_header = FrameHeader(
        frame_start,
        frame_size,
        num_atoms,
        atom_data_start,
        atom_data_size,
        num_bonds,
        bond_data_start,
        bond_data_size
    )

    frame = Frame(
        :box_min => (x_lo, y_lo, z_lo),
        :box_max => (x_hi, y_hi, z_hi)
    )

    return frame_header, frame
end

function read_atoms{T}(::Type{T}, ::Type{LAMMPSDAT}, chunk::Chunk, atoms::Atoms, indices::UnitRange{Int})
    molecules = atoms.props[:molecule]::Vector{Int}
    types = atoms.props[:type]::Vector{Symbol}
    charges = atoms.props[:charge]::Vector{T}
    coords = atoms.props[:coords]::Matrix{T}

    for i in indices
        read_whitespace(chunk)
        read_nonwhitespace(chunk)
        read_whitespace(chunk)

        molecule = read_int(chunk)
        read_whitespace(chunk)

        atom_type = read_symbol(chunk)
        read_whitespace(chunk)

        charge = read_float(T, chunk)
        read_whitespace(chunk)
        
        x = read_float(T, chunk)
        read_whitespace(chunk)
        
        y = read_float(T, chunk)
        read_whitespace(chunk)
        
        z = read_float(T, chunk)
        read_newline(chunk)

        molecules[i] = molecule
        types[i] = atom_type
        charges[i] = charge
        coords[1, i] = x
        coords[2, i] = y
        coords[3, i] = z
    end
end

function read_bonds(::Type{LAMMPSDAT}, chunk::Chunk, bonds::Matrix{Int}, indices::UnitRange{Int})
    num_bonds = size(bonds, 2)

    for i in indices
        read_whitespace(chunk)
        read_nonwhitespace(chunk)
        read_whitespace(chunk)

        read_nonwhitespace(chunk)
        read_whitespace(chunk)

        atom1 = read_int(chunk)
        read_whitespace(chunk)

        atom2 = read_int(chunk)
        read_newline(chunk)

        bonds[1, i] = atom1
        bonds[2, i] = atom2
    end
end

function write_atoms(::Type{LAMMPSDAT}, io::IO, frame::Frame)
    haskey(frame, :atoms) || throw(ArgumentError("Need frame atoms"))
    haskey(frame, :box_min) || throw(ArgumentError("Need frame box min"))
    haskey(frame, :box_max) || throw(ArgumentError("Need frame box max"))
    atoms = frame[:atoms]

    haskey(atoms, :coords) || throw(ArgumentError("Need atom coords"))
    haskey(atoms, :molecule) || throw(ArgumentError("Need atom molecule IDs)"))
    haskey(atoms, :type) || throw(ArgumentError("Need atom types"))
    haskey(atoms, :charge) || throw(ArgumentError("Need atom charges)"))
    haskey(atoms, :coords) || throw(ArgumentError("Need atom coords"))

    num_atoms = length(atoms)
    unique_elements = sort!(unique(atoms[:type]))
    element_ids = Dict{Symbol, Int}(zip(unique_elements, eachindex(unique_elements)))
    element_masses = Dict(:C => 12.0107, :N => 12.0107, :Li => 6.941, :H => 1.00794)
    xlo, ylo, zlo = frame[:box_min]
    xhi, yhi, zhi = frame[:box_max]

    println(io,
        """
        LAMMPS data file
        $num_atoms atoms
        $(length(unique_elements)) atom types
        $xlo $xhi xlo xhi
        $ylo $yhi ylo yhi
        $zlo $zhi zlo zhi

        Masses
        """
    )

    for (id, element) in enumerate(unique_elements)
        println(io, id, ' ', element_masses[element], " # ", element)
    end
    println(io)

    println(io, "Atoms")
    println(io)

    for (i, atom) in enumerate(atoms)
        print(io,
            i, ' ',
            atom[:molecule], ' ',
            element_ids[atom[:type]], ' ',
            atom[:charge], ' ',
            atom[:coords][1], ' ',
            atom[:coords][2], ' ',
            atom[:coords][3], '\n'
        )
    end
end
