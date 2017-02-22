import Base: show, display, getindex, read

###########################
#    Types
###########################

abstract type AtomFileType end

immutable FrameHeader
    frame_start::Int64
    frame_size::Int64

    num_atoms::Int64
    atom_data_start::Int64
    atom_data_size::Int64

    num_bonds::Int64
    bond_data_start::Int64
    bond_data_size::Int64
end

FrameHeader(frame_start::Int64, frame_size::Int64, num_atoms::Int64, atom_data_start::Int64, atom_data_size::Int64) = FrameHeader(frame_start, frame_size, num_atoms, atom_data_start, atom_data_size, 0, 0, 0)

immutable AtomFile{F <: AtomFileType}
    filename::AbstractString
    frame_headers::Vector{FrameHeader}
    frames::Vector{Frame}

    AtomFile(filename::AbstractString) = new(filename, Vector{FrameHeader}(), Vector{Frame}())
end

function AtomFile(filename::AbstractString)
    base, ext = splitext(filename)
    for file_type in subtypes(AtomFileType)
        if ext ∈ file_exts(file_type)
            return AtomFile{file_type}(filename)
        end
    end

    throw(ArgumentError("Unknown file type $ext"))
end

AtomFile{F <: AtomFileType}(::Type{F}, filename::AbstractString) = AtomFile{F}(filename)

type Chunk
    data::Vector{UInt8}
    i::Int64
    
    Chunk(num_bytes::Int64) = new(Vector{UInt8}(num_bytes), 1)
    Chunk(data::Vector{UInt8}) = new(data, 1)
end

################################
#    File reading functions
################################

function read_byte(chunk::Chunk)
    b = chunk.data[chunk.i]
    chunk.i += 1
    return b
end

preview_byte(chunk::Chunk) = chunk.data[chunk.i]

next_byte(chunk::Chunk) = chunk.i += 1

read_newline(chunk::Chunk) = while read_byte(chunk) != 0x0a end

function read_newlines(chunk::Chunk, num_newlines::Int)
    s = 0
    while s < num_newlines
        if read_byte(chunk) == 0x0a
            s += 1
        end
    end
end

function read_float{T}(::Type{T}, chunk::Chunk)
    a = T(0)
    b = T(0)
    p = 1
    
    c = preview_byte(chunk)
    if c == 0x2d
        p = -1
        next_byte(chunk)
    end

    c == 0x2b && next_byte(chunk)
    
    while true
        c = preview_byte(chunk)
        (c < 0x30 || c > 0x39) && break
        next_byte(chunk)
        a = 10a + (c - 0x30)
    end
    
    c != 0x2e && return p*a
    next_byte(chunk)
    
    i = 1
    while true
        c = preview_byte(chunk)
        (c < 0x30 || c > 0x39) && break
        next_byte(chunk)
        a = 10a + (c - 0x30)
        i *= 10
    end
    
    c != 0x45 && c != 0x65 && return (p*a/i)
    next_byte(chunk)
    
    b = read_int(chunk)
    return (p*a/i) * T(10)^b
end

function read_int(chunk::Chunk)
    n = 0
    p = 1
    
    c = preview_byte(chunk)
    if c == 0x2d
        p = -1
        next_byte(chunk)
    end

    c == 0x2b && next_byte(chunk)
    
    while true
        c = preview_byte(chunk)
        (c < 0x30 || c > 0x39) && break
        next_byte(chunk)
        n = 10n + (c - 0x30)
    end
    
    return p*n
end

function read_symbol(chunk::Chunk, atom_types::Matrix{UInt8}, i::Int)
    n = size(atom_types, 1)

    j = 1
    while j <= n
        c = preview_byte(chunk)
        (c == 0x20 || c == 0x09) && break
        atom_types[j, i] = c
        next_byte(chunk)
        j += 1
    end

    while j <= n
        atom_types[j, i] = 0x20
        j += 1
    end
end

function read_symbol(chunk::Chunk)
    i1 = chunk.i
    while true
        c = preview_byte(chunk)
        if (c == 0x20 || c == 0x09)
            return ccall(:jl_symbol_n, Ref{Symbol}, (Ptr{UInt8}, Int32), pointer(chunk.data, i1), Int32(chunk.i - i1))
        end
        next_byte(chunk)
    end
end

function read_nonwhitespace(chunk::Chunk)
    while true
        c = preview_byte(chunk)
        (c == 0x20 || c == 0x09 || c == 0x23 || c == 0x0a) && return
        next_byte(chunk)
    end
end

function read_whitespace(chunk::Chunk)
    while true
        c = preview_byte(chunk)
        c != 0x20 && c != 0x09 && c != 0x23 && return
        next_byte(chunk)
    end
end

function read_line(chunk::Chunk)
    i = chunk.i
    read_newline(chunk)
    #ascii(pointer(chunk.data, i), chunk.i - i - 1)
    String(chunk.data[i:chunk.i - 1])
end

function read_through_regex(chunk::Chunk, r::Regex)
    while true
        s = read_line(chunk)
        ismatch(r, s) && return match(r, s)
    end
end
    
function write_atoms{F <: AtomFileType}(file::AtomFile{F}, frames::Vector{Frame})
    filename = file.filename
    io = open(filename, "w")

    for frame in frames
        write_atoms(F, io, frame)
    end

    close(io)
end

write_atoms{F <: AtomFileType}(file::AtomFile{F}, frame::Frame) = write_atoms(file, [frame])
write_atoms(filename::AbstractString, frames::Vector{Frame}) = write_atoms(AtomFile(filename), frames)
write_atoms(filename::AbstractString, frame::Frame) = write_atoms(AtomFile(filename), frame)

function read_atoms(filename::AbstractString, end_frames::Int = 0)
    file = AtomFile(filename)
    num_frames = length(read_frame_headers(file))
    end_frames == 0 && (end_frames = num_frames)
    frames = read_atoms(file, frames = num_frames-end_frames+1:num_frames)
    end_frames == 1 ? frames[1] : frames
end

function read_atoms{F <: AtomFileType}(file::AtomFile{F}; frames::AbstractVector{Int} = eachindex(file.frame_headers))
    max_buffer_size = 1*1024*1024*1024

    length(file.frame_headers) == 0 && read_frame_headers(file)

    if length(frames) == 0
        frames = eachindex(file.frame_headers)
    end

    filename = file.filename
    frame_headers = file.frame_headers[frames]
    frame_data = file.frames[frames]

    num_frames = length(frame_headers)
    max_frame_size = maximum(f -> f.frame_size, frame_headers)
    
    max_buffer_size >= max_frame_size || error("Enlarge buffer size.")

    num_chunks = floor(Int, max_buffer_size / max_frame_size)
    num_chunks = min(num_chunks, num_frames)

    chunks = [Chunk(max_frame_size) for _ in 1:num_chunks]

    atom_indices = cumsum([f.num_atoms for f in frame_headers])
    total_atoms = atom_indices[end]
    atom_indices .-= atom_indices[1] - 1

    atoms = Atoms(F, total_atoms)

    ios = open(filename, "r")

    for (i, f) in enumerate(frame_headers)
        atom_start = atom_indices[i]
        atom_end = atom_start + f.num_atoms - 1
        atom_range = atom_start:atom_end

        frame_data[i][:atoms] = atoms[atom_range]
        j = mod1(i, num_chunks)
        # Seek indexing is 0-based
        seek(ios, f.atom_data_start - 1)
        chunks[j].i = 1
        readbytes!(ios, chunks[j].data, f.atom_data_size, all = false)
        read_atoms(Float64, F, chunks[j], atoms, atom_range)

        if frame_headers[i].num_bonds > 0
            bond_range = 1:f.num_bonds
            bonds = Matrix{Int}(2, f.num_bonds)
            frame_data[i][:bonds] = bonds

            seek(ios, f.bond_data_start - 1)
            chunks[j].i = 1
            readbytes!(ios, chunks[j].data, f.bond_data_size, all = false)
            read_bonds(F, chunks[j], bonds, bond_range)
        end
    end

    close(ios)

    return frame_data
end

function read_frame_headers{F <: AtomFileType}(file::AtomFile{F})
    file_size = filesize(file.filename)
    chunk = Chunk(file_size)

    ios = open(file.filename, "r")
    readbytes!(ios, chunk.data, file_size, all = true)
    close(ios)
    
    empty!(file.frame_headers)
    empty!(file.frames)
    while chunk.i < file_size
        frame_header, frame = read_frame_header(Float64, F, chunk)
        push!(file.frame_headers, frame_header)
        push!(file.frames, frame)
    end

    return file.frames
end
