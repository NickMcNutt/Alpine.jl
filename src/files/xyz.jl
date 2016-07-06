type XYZ <: AtomFileType end

file_exts(::Type{XYZ}) = (".xyz",)

function Atoms(::Type{XYZ}, num_atoms::Int)
    Atoms(num_atoms,
        :type => Vector{Symbol}(num_atoms),
        :coords => Matrix{Float64}(3, num_atoms),
        :molecule => Vector{Int}(num_atoms)
    )
end

function read_frame_header{T}(::Type{T}, ::Type{XYZ}, chunk::Chunk)
    frame_start = chunk.i
    
    num_atoms = read_int(chunk)
    read_newlines(chunk, 2)
    
    atom_data_start = chunk.i
    read_newlines(chunk, num_atoms)
    atom_data_size = chunk.i - atom_data_start

    frame_size = chunk.i - frame_start

    frame_header = FrameHeader(frame_start, frame_size, num_atoms, atom_data_start, atom_data_size)

    frame = Frame()

    return frame_header, frame
end

function read_atoms{T}(::Type{T}, ::Type{XYZ}, chunk::Chunk, atom::Atoms, indices::UnitRange{Int})
    types = atom.props[:type]::Vector{Symbol}
    coords = atom.props[:coords]::Matrix{T}
    molecules = atom.props[:molecule]::Vector{Int}

    for i in indices
        atom_type = read_symbol(chunk)
        read_whitespace(chunk)
        
        x = read_float(T, chunk)
        read_whitespace(chunk)
        
        y = read_float(T, chunk)
        read_whitespace(chunk)
        
        z = read_float(T, chunk)
        read_whitespace(chunk)

        molecule = read_int(chunk)
        read_newline(chunk)
        
        types[i] = atom_type
        coords[1, i] = x
        coords[2, i] = y
        coords[3, i] = z
        molecules[i] = molecule
    end
end
