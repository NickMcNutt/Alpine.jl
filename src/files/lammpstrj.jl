type LAMMPSTRJ <: AtomFileType end

typealias TRJ LAMMPSTRJ

file_exts(::Type{TRJ}) = (".lammpstrj", ".trj")

function Atoms(::Type{TRJ}, num_atoms::Int)
    Atoms(num_atoms,
        :type => Vector{Symbol}(num_atoms),
        :molecule => Vector{Int}(num_atoms),
        :coords => Matrix{Float64}(3, num_atoms),
        :charge => Vector{Float64}(num_atoms),
        :energy => Vector{Float64}(num_atoms)
    )
end

function read_frame_header{T}(::Type{T}, ::Type{TRJ}, chunk::Chunk)
    frame_start = chunk.i
    
    read_newline(chunk)
    
    timestep = read_int(chunk)
    read_newlines(chunk, 2)
    
    num_atoms = read_int(chunk)
    read_newlines(chunk, 2)
    
    xl = read_float(T, chunk)
    read_whitespace(chunk)
    xh = read_float(T, chunk)
    read_newline(chunk)
    
    yl = read_float(T, chunk)
    read_whitespace(chunk)
    yh = read_float(T, chunk)
    read_newline(chunk)
    
    zl = read_float(T, chunk)
    read_whitespace(chunk)
    zh = read_float(T, chunk)
    read_newlines(chunk, 2)

    atom_data_start = chunk.i
    read_newlines(chunk, num_atoms)
    atom_data_size = chunk.i - atom_data_start

    frame_size = chunk.i - frame_start
    
    frame_header = FrameHeader(
        frame_start,
        frame_size,
        num_atoms,
        atom_data_start,
        atom_data_size
    )

    frame = Frame(
        :timestep => timestep,
        :box_min => (xl, yl, zl),
        :box_max => (xh, yh, zh),
    )

    return frame_header, frame
end

function read_atoms{T}(::Type{T}, ::Type{TRJ}, chunk::Chunk, atoms::Atoms, indices::UnitRange{Int})
    types = atoms.props[:type]::Vector{Symbol}
    molecules = atoms.props[:molecule]::Vector{Int}
    coords = atoms.props[:coords]::Matrix{T}
    charges = atoms.props[:charge]::Vector{T}
    energies = atoms.props[:energy]::Vector{T}

    for i in indices
        atom_type = read_symbol(chunk)
        read_whitespace(chunk)

        molecule = read_int(chunk)
        read_whitespace(chunk)
        
        x = read_float(T, chunk)
        read_whitespace(chunk)
        
        y = read_float(T, chunk)
        read_whitespace(chunk)
        
        z = read_float(T, chunk)
        read_whitespace(chunk)
        
        charge = read_float(T, chunk)
        read_whitespace(chunk)
        
        energy = read_float(T, chunk)
        read_newline(chunk)

        types[i] = atom_type
        molecules[i] = molecule
        coords[1, i] = x
        coords[2, i] = y
        coords[3, i] = z
        charges[i] = charge
        energies[i] = energy
    end
end
