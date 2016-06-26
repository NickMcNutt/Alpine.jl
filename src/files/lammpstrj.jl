function read_frame_header{T}(chunk::Chunk, ::Type{T})
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
    frame_end = chunk.i - 1
    
    Frame(num_atoms, timestep, (xh - xl, yh - yl, zh - zl), frame_start, frame_end, atom_data_start)
end

function readlineatoms{T}(chunk::Chunk, atom_data::AtomData{T}, i::Int)
    atom_type = read_symbol(chunk)
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
    
    atom_data.types[i] = atom_type
    atom_data.coords[1, i] = x
    atom_data.coords[2, i] = y
    atom_data.coords[3, i] = z
    atom_data.charges[i] = charge
    atom_data.energies[i] = energy
end

function read_frame{T}(chunk::Chunk, atom_data::AtomData{T}, frame_index::Int)
    start_index = atom_data.indices[frame_index]
    end_index = start_index + atom_data.frames[frame_index].num_atoms - 1
    for i in start_index:end_index
        readlineatoms(chunk, atom_data, i)
    end
end

