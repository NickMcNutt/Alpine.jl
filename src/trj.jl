# A quick hack to read in a LAMMPS trajectory file
function readtrj(filename::AbstractString)
    ios = open(filename, "r")

    num_lines = countlines(ios)
    seekstart(ios)

    regex_num_atoms = r"item: *number of atoms"i
    regex_atoms = r"item: *atoms"i

    while !eof(ios)
        ismatch(regex_num_atoms, readline(ios)) && break
    end

    num_atoms = parse(Int, readline(ios))
    num_frames = div(num_lines, (num_atoms + 9))
    total_atoms = num_atoms * num_frames
    #println("$num_frames frames")
    #println("$num_atoms atoms")

    seekstart(ios)

    box = Vector{Tuple{Float64, Float64, Float64}}(num_frames)
    atom_types = Vector{Symbol}(total_atoms)
    atom_coords = Matrix{Float64}(3, total_atoms)
    atom_charges = Vector{Float64}(total_atoms)
    atom_energies = Vector{Float64}(total_atoms)
    
    @inbounds for f in 1:num_frames
        for i in 1:5
            readline(ios)
        end

        xlo, xhi = map(x -> parse(Float64, x), split(readline(ios))[1:2])
        ylo, yhi = map(x -> parse(Float64, x), split(readline(ios))[1:2])
        zlo, zhi = map(x -> parse(Float64, x), split(readline(ios))[1:2])

        box[f] = (xhi - xlo, yhi - ylo, zhi - zlo)

        readline(ios)

        @inbounds for i in 1:num_atoms
            j = (f-1) * num_atoms + i
            s = split(readline(ios), ' ', limit = 12, keep = false)
            atom_types[j] = Symbol(s[1])
            atom_coords[1, j] = parse(Float64, s[2])
            atom_coords[2, j] = parse(Float64, s[3])
            atom_coords[3, j] = parse(Float64, s[4])
            atom_charges[j] = parse(Float64, s[5])
            atom_energies[j] = parse(Float64, s[6])
        end

        print('\r')
        print("$f/$num_frames")
        flush(STDOUT)
    end

    close(ios)

    return Dict(:num_frames => num_frames,
                :num_atoms => num_atoms,
                :box => box,
                :types => atom_types,
                :coords => atom_coords,
                :charges => atom_charges,
                :energies => atom_energies)
end
