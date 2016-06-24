import Base: show, display, getindex

###########################
#    Types
###########################

type Chunk
    data::Vector{UInt8}
    i::Int64
    
    Chunk(num_bytes::Int64) = new(Vector{UInt8}(num_bytes), 1)
    Chunk(data::Vector{UInt8}) = new(data, 1)
end

immutable Frame{T}
    num_atoms::Int
    timestep::Int
    box_dims::Tuple{T, T, T}
    frame_start::Int64
    frame_end::Int64
    atom_data_start::Int64
end

immutable AtomData{T}
    frames::Vector{Frame{T}}
    indices::Vector{Int64}
    types::Vector{Symbol}
    coords::Matrix{T}
    charges::Vector{T}
    energies::Vector{T}

    function AtomData(frames::Vector{Frame{T}}, indices::Vector{Int64}, total_atoms::Int64)
        new(frames, indices, Vector{Symbol}(total_atoms), Matrix{T}(3, total_atoms), Vector{T}(total_atoms), Vector{T}(total_atoms))
    end
end

function AtomData{T}(frames::Vector{Frame{T}}, indices::Vector{Int64}, total_atoms::Int64)
    AtomData{T}(frames, indices, total_atoms)
end

###########################
#    Display functions
###########################

show(io::IO, frame::Frame) = @printf(io, "t = %d  %d atoms  box = %f %f %f", frame.timestep, frame.num_atoms, frame.box_dims[1], frame.box_dims[2], frame.box_dims[3])
display(frame::Frame) = show(frame)

function show(io::IO, ad::AtomData)
    println(io, "$(length(ad.frames)) frames with atom types $(tuple(unique(ad.types)...))\n")
    for frame in ad.frames
        show(io, frame)
        println(io)
    end
end

display(ad::AtomData) = show(ad)

###########################
#    Atom data functions
###########################

function get_atom_indices{T}(ad::AtomData{T}, frame_indices::AbstractVector{Int})
    num_frames = length(frame_indices)
    atom_indices = Vector{UnitRange{Int}}(num_frames)
    for (i, f) in enumerate(frame_indices)
        start_index = ad.indices[f]
        end_index = start_index + ad.frames[f].num_atoms - 1
        atom_indices[i] = start_index:end_index
    end

    return atom_indices
end

function get_atom_indices{T}(ad::AtomData{T}, frame_index::Int)
    start_index = ad.indices[frame_index]
    end_index = start_index + ad.frames[frame_index].num_atoms - 1
    start_index:end_index
end

function coords(ad::AtomData, frame_indices::AbstractVector{Int})
    atom_indices = get_atom_indices(ad, frame_indices)
    ad.coords[:, vcat(atom_indices...)]
end

coords(ad::AtomData, frame_index::Int) = ad.coords[:, get_atom_indices(ad, frame_index)]
coords(ad::AtomData) = length(ad.frames) == 1 ? coords(ad, 1) : coords(ad, 1:length(ad.frames))

function types(ad::AtomData, frame_indices::AbstractVector{Int})
    atom_indices = get_atom_indices(ad, frame_indices)
    ad.types[vcat(atom_indices...)]
end

types(ad::AtomData, frame_index::Int) = ad.types[get_atom_indices(ad, frame_index)]
types(ad::AtomData) = length(ad.frames) == 1 ? types(ad, 1) : types(ad, 1:length(ad.frames))

function charges(ad::AtomData, frame_indices::AbstractVector{Int})
    atom_indices = get_atom_indices(ad, frame_indices)
    ad.charges[vcat(atom_indices...)]
end

charges(ad::AtomData, frame_index::Int) = ad.charges[get_atom_indices(ad, frame_index)]
charges(ad::AtomData) = length(ad.frames) == 1 ? charges(ad, 1) : charges(ad, 1:length(ad.frames))

function energies(ad::AtomData, frame_indices::AbstractVector{Int})
    atom_indices = get_atom_indices(ad, frame_indices)
    ad.energies[vcat(atom_indices...)]
end

energies(ad::AtomData, frame_index::Int) = ad.energies[get_atom_indices(ad, frame_index)]
energies(ad::AtomData) = length(ad.frames) == 1 ? energies(ad, 1) : energies(ad, 1:length(ad.frames))

num_atoms(ad::AtomData, frame_indices::AbstractVector{Int}) = [frame.num_atoms for frame in ad.frames[frame_indices]]
num_atoms(ad::AtomData, frame_index::Int) = ad.frames[frame_index].num_atoms
num_atoms(ad::AtomData) = length(ad.frames) == 1 ? num_atoms(ad, 1) : num_atoms(ad, 1:length(ad.frames))

timestep(ad::AtomData, frame_indices::AbstractVector{Int}) = [frame.timestep for frame in ad.frames[frame_indices]]
timestep(ad::AtomData, frame_index::Int) = ad.frames[frame_index].timestep
timestep(ad::AtomData) = length(ad.frames) == 1 ? timestep(ad, 1) : timestep(ad, 1:length(ad.frames))

box_dims(ad::AtomData, frame_indices::AbstractVector{Int}) = [frame.box_dims for frame in ad.frames[frame_indices]]
box_dims(ad::AtomData, frame_index::Int) = ad.frames[frame_index].box_dims
box_dims(ad::AtomData) = length(ad.frames) == 1 ? box_dims(ad, 1) : box_dims(ad, 1:length(ad.frames))

################################
#    File reading functions
################################

function readbyte(chunk::Chunk)
    b = chunk.data[chunk.i]
    chunk.i += 1
    return b
end

previewbyte(chunk::Chunk) = chunk.data[chunk.i]

nextbyte(chunk::Chunk) = chunk.i += 1

readnewline(chunk::Chunk) = while readbyte(chunk) != 0x0a end

function readnewlines(chunk::Chunk, num_newlines::Int)
    s = 0
    while s < num_newlines
        if readbyte(chunk) == 0x0a
            s += 1
        end
    end
end

function readfloat{T}(::Type{T}, chunk::Chunk)
    a = T(0)
    b = T(0)
    p = 1
    
    c = previewbyte(chunk)
    if c == 0x2d
        p = -1
        nextbyte(chunk)
    end
    
    while true
        c = previewbyte(chunk)
        (c < 0x30 || c > 0x39) && break
        nextbyte(chunk)
        a = 10a + (c - 0x30)
    end
    
    c != 0x2e && return p*a
    nextbyte(chunk)
    
    i = 1
    while true
        c = previewbyte(chunk)
        (c < 0x30 || c > 0x39) && break
        nextbyte(chunk)
        a = 10a + (c - 0x30)
        i *= 10
    end
    
    c != 0x45 && c != 0x65 && return (p*a/i)
    nextbyte(chunk)
    
    b = readint(chunk)
    return (p*a/i) * T(10)^b
end

function readint(chunk::Chunk)
    n = 0
    p = 1
    
    c = previewbyte(chunk)
    if c == 0x2d
        p = -1
        nextbyte(chunk)
    end
    
    while true
        c = previewbyte(chunk)
        (c < 0x30 || c > 0x39) && break
        nextbyte(chunk)
        n = 10n + (c - 0x30)
    end
    
    return p*n
end

function readsymbol(chunk::Chunk, atom_types::Matrix{UInt8}, i::Int)
    n = size(atom_types, 1)

    j = 1
    while j <= n
        c = previewbyte(chunk)
        (c == 0x20 || c == 0x09) && break
        atom_types[j, i] = c
        nextbyte(chunk)
        j += 1
    end

    while j <= n
        atom_types[j, i] = 0x20
        j += 1
    end
end

function readsymbol(chunk::Chunk)
    i1 = chunk.i
    while true
        c = previewbyte(chunk)
        if (c == 0x20 || c == 0x09)
            return ccall(:jl_symbol_n, Ref{Symbol}, (Ptr{UInt8}, Int32), pointer(chunk.data, i1), Int32(chunk.i - i1))
        end
        nextbyte(chunk)
    end
end

function readnonblank(chunk::Chunk)
    while true
        c = previewbyte(chunk)
        (c == 0x20 || c == 0x09 || c == 0x0a) && return
        nextbyte(chunk)
    end
end

function readblank(chunk::Chunk)
    while true
        c = previewbyte(chunk)
        c != 0x20 && c != 0x09 && return
        nextbyte(chunk)
    end
end

function readframeheader{T}(chunk::Chunk, ::Type{T})
    frame_start = chunk.i
    
    readnewline(chunk)
    
    timestep = readint(chunk)
    readnewlines(chunk, 2)
    
    num_atoms = readint(chunk)
    readnewlines(chunk, 2)
    
    xl = readfloat(T, chunk)
    readblank(chunk)
    xh = readfloat(T, chunk)
    readnewline(chunk)
    
    yl = readfloat(T, chunk)
    readblank(chunk)
    yh = readfloat(T, chunk)
    readnewline(chunk)
    
    zl = readfloat(T, chunk)
    readblank(chunk)
    zh = readfloat(T, chunk)
    readnewlines(chunk, 2)

    atom_data_start = chunk.i
    readnewlines(chunk, num_atoms)
    frame_end = chunk.i - 1
    
    Frame(num_atoms, timestep, (xh - xl, yh - yl, zh - zl), frame_start, frame_end, atom_data_start)
end

function readlineatoms{T}(chunk::Chunk, atom_data::AtomData{T}, i::Int)
    atom_type = readsymbol(chunk)
    readnonblank(chunk)
    readblank(chunk)
    
    x = readfloat(T, chunk)
    readblank(chunk)
    
    y = readfloat(T, chunk)
    readblank(chunk)
    
    z = readfloat(T, chunk)
    readblank(chunk)
    
    charge = readfloat(T, chunk)
    readblank(chunk)
    
    energy = readfloat(T, chunk)
    readnewline(chunk)
    
    atom_data.types[i] = atom_type
    atom_data.coords[1, i] = x
    atom_data.coords[2, i] = y
    atom_data.coords[3, i] = z
    atom_data.charges[i] = charge
    atom_data.energies[i] = energy
end

function readframeatoms{T}(chunk::Chunk, atom_data::AtomData{T}, frame_index::Int)
    start_index = atom_data.indices[frame_index]
    end_index = start_index + atom_data.frames[frame_index].num_atoms - 1
    for i in start_index:end_index
        readlineatoms(chunk, atom_data, i)
    end
end

function read_atomfile{T}(filename::AbstractString, frames::Vector{Frame{T}})
    max_buffer_size = 1*1024*1024*1024

    num_frames = length(frames)
    frame_sizes = [frame.frame_end - frame.atom_data_start + 1 for frame in frames]
    max_frame_size = maximum(frame_sizes)
    
    max_buffer_size < max_frame_size && error("Enlarge buffer size.")

    num_chunks = floor(Int, max_buffer_size / max_frame_size)
    num_chunks = min(num_chunks, num_frames)

    chunk_size = max_frame_size
    chunks = [Chunk(chunk_size) for _ in 1:num_chunks]

    indices = cumsum(map(f -> Int64(f.num_atoms), frames))
    total_atoms = indices[end]
    indices .-= indices[1] - 1

    atom_data = AtomData(frames, indices, total_atoms)

    ios = open(filename, "r")

    for (i, frame) in enumerate(frames)
        j = mod1(i, num_chunks)
        seek(ios, frame.atom_data_start - 1)
        chunks[j].i = 1
        readbytes!(ios, chunks[j].data, frame_sizes[i], all = false)
        readframeatoms(chunks[j], atom_data, i)
    end

    close(ios)

    return atom_data
end

function read_atomfile_structure(filename::AbstractString)
    file_size = stat(filename).size
    chunk = Chunk(file_size)

    frames = Vector{Frame{Float64}}()

    ios = open(filename, "r")
    readbytes!(ios, chunk.data, file_size, all = false)
    close(ios)
    
    while chunk.i < file_size
        frame = readframeheader(chunk, Float64)
        push!(frames, frame)
    end

    return frames
end
