import Base: open, close, length, endof, first, last, start, next, done, push!

type AtomFile{T}
    name::String
end

type AtomStream{T}
    ios::IOStream
    num_frames::Int64
    frame_positions::Array{Tuple{Int64, Int64}, 1}
    num_atoms::Array{Int64, 1}

    AtomStream(ios::IOStream) = new(ios, zero(Int64), Int64[], Int64[])
end

#convert{T}(::Type{T:<AtomFile}, x::String) = T(x)
start(s::AtomStream) = 1
next(s::AtomStream, state::Int) = (getindex(s, state), state + 1)
done(s::AtomStream, state::Int) = (state > s.num_frames)
length(s::AtomStream) = s.num_frames
endof(s::AtomStream) = length(s)
first(s::AtomStream) = getindex(s, 1)
last(s::AtomStream) = getindex(s, endof(s))
getindex(s::AtomStream, i::Int) = (properties...) -> read(s, i, properties)

function open(f::Function, file::AtomFile{:xyz})
    ios = open(file.name)
    stream = AtomStream{:xyz}(ios)
    while !eof(ios)
        num_entries::Int64 = parse(Int64, readline(ios))
        push!(stream.num_atoms, num_entries)
        readline(ios)
        start_position::Int64 = position(ios)
        for i = 1:num_entries
            readline(ios)
        end
        push!(stream.frame_positions, (start_position, position(ios) - start_position - 1))
    end

    stream.num_frames = length(stream.frame_positions)

    result = f(stream)
    close(ios)
    
    return result
end

close(stream::AtomStream) = close(stream.ios)

function open(f::Function, file::AtomFile{:trj})
    stream = open(file)
    result = f(stream)
    close(stream)

    return result
end

function open(file::AtomFile{:trj})
    ios = open(file.name)
    stream = AtomStream{:trj}(ios)

    regex_num_atoms = r"item: *number of atoms"i
    regex_atoms = r"item: *atoms"i

    while !eof(ios)
        while !eof(ios)
            if ismatch(regex_num_atoms, readline(ios))
                break
            end
        end

        num_atoms::Int64 = parse(Int64, readline(ios))
        push!(stream.num_atoms, num_atoms)

        while !eof(ios)
            if ismatch(regex_atoms, readline(ios))
                break
            end
        end

        start_position::Int64 = position(ios)
        for i = 1:num_atoms
            readline(ios)
        end
        push!(stream.frame_positions, (start_position, position(ios) - start_position - 1))
    end

    stream.num_frames = length(stream.frame_positions)

    return stream
end
        
function read(stream::AtomStream, frame::Int, properties)
    num_atoms = stream.num_atoms[frame]
    seek(stream.ios, stream.frame_positions[frame][1])
    num_bytes = stream.frame_positions[frame][2]
    bytes = Vector{Uint8}(num_bytes)
    bytes_read = readbytes!(stream.ios, bytes, num_bytes)
    if num_bytes != bytes_read
        error("num_bytes != bytes_read")
    end
    lines = split(bytestring(bytes), "\n")

    parts = split(strip(lines[1]), r" +|\t")
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

    props_allocate = Dict{Symbol, Function}(
        :types => () -> Array{Symbol}(num_atoms),
        :coords => () -> Array{Float64}(3, num_atoms),
        :misc => () -> Array{Any}(num_misc_fields, num_atoms)
    )

    props_mem = Dict{Symbol, Any}([p => f() for (p, f) in props_allocate])

    for i = 1:length(lines)
        line = strip(chomp(lines[i]))
        parts = split(line, r" +|\t")

        if :num ∈ properties
            props_mem[:num] = num_atoms
        end

        if :types ∈ properties
            props_mem[:types][i] = symbol(string(parts[1]))
        end

        if :coords ∈ properties
            props_mem[:coords][1:3, i] = map(x -> parse(Float64, x), parts[2:4])
        end

        if :misc ∈ properties
            if num_misc_fields > 0
                if comment_present
                    props_mem[:misc][:, i] = parts[6:(6 + num_misc_fields - 1)]
                else
                    props_mem[:misc][:, i] = parts[5:(5 + num_misc_fields - 1)]
                end
            end
        end
    end

    return map(k -> props_mem[k], properties)
end

function readPOSFile(filename)
    # Obtain number of data entries.  4 bytes per Float32 and 4 numerical values per entry
    num_entries::Int = filesize(filename_input) / 4 / 4

    # Read file into a (num_entries x 4) matrix
    file_input = open(filename_input, "r")
    data = read(file_input, Float32, (4, num_entries))'
    close(file_input)

    # POS data is in big-endian format.  Who's idea was this?  Convert to little-endian.
    map!((x) -> ntoh(x), data)

    atom_coords = data[:, 1:3]
    atom_misc = data[:, 4]

    return (num_entries, atom_coords, atom_misc)
end

function readPOSFile(filename, atom_types)
    readPOSFile(filename)
    # finish this
end


newFilename(filename, new_extension) = match(r"(.*[.])", filename).captures[1] * new_extension

