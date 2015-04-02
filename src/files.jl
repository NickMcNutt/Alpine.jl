import Base: open, length, start, next, done

type XYZFile
    name::String
end

type XYZStream
    ios::IOStream
    num_frames::Int64
    frame_positions::Array{(Int64, Int64), 1}
    num_atoms::Array{Int64, 1}

    XYZStream(ios::IOStream) = new(ios, zero(Int64), Int64[], Int64[])
end

convert(::Type{XYZFile}, x::String) = XYZFile(x)
start(s::XYZStream) = 1
next(s::XYZStream, state::Int) = (() -> readFrame(s, state), state + 1)
done(s::XYZStream, state::Int) = (state > s.num_frames)
length(s::XYZStream) = s.num_frames
getindex(s::XYZStream, i::Int) = () -> readFrame(s, i)

function open(f::Function, xyz_file::XYZFile)
    ios = open(xyz_file.name)
    stream = XYZStream(ios)
    while !eof(ios)
        num_entries::Int64 = parse(Int, readline(ios))
        push!(stream.num_atoms, num_entries)
        readline(ios)
        start_position::Int64 = position(ios)
        for i = 1:num_entries
            readline(ios)
        end
        push!(stream.frame_positions, (start_position, position(ios) - start_position - 1))
    end

    stream.num_frames = length(stream.frame_positions)

    f(stream)
    close(ios)
end

function readFrame(stream::XYZStream, frame::Int)
    num_atoms = stream.num_atoms[frame]
    seek(stream.ios, stream.frame_positions[frame][1])
    num_bytes = stream.frame_positions[frame][2]
    bytes = Vector{Uint8}(num_bytes)
    bytes_read = readbytes!(stream.ios, bytes, num_bytes)
    if num_bytes != bytes_read
        error("num_bytes != bytes_read")
    end
    lines = split(bytestring(bytes), "\n")

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

