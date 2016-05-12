function section{T<:AbstractFloat}(x::T, b::T, n::Int)
    fld(n * (x + b/2), b)
end

function density{T<:AbstractFloat}(file_out::IO, C::AbstractMatrix{T}, n::Int, bw::T = zero(Float64))
    num_atoms = size(C, 2)

    b::T = 2 * maxabs(C)
    println(b)

    if (bw != zero(T))
        bw < b && error("Your box specification is not big enough to contain all atoms.")

        b = bw
    else
        b *= (1 + 1/n)
    end

    println(file_out, "density")
    println(file_out, "$b   "^3, "90.00000   "^3)
    println(file_out, "$n   "^3)

    indices, counts = Int[], Int[]
    @inbounds for i in 1:num_atoms
        xs = section(C[1, i], b, n)
        ys = section(C[2, i], b, n)
        zs = section(C[3, i], b, n)
        index = zs + n*ys + n*n*xs + 1
        push!(indices, index)
        push!(counts, 1)
    end
    ρ = sparsevec(indices, counts)
    full_ρ = full(ρ)

    remaining_zeros = n^3 - length(full_ρ)

    print_joined(file_out, full_ρ, '\n')
    println(file_out, "\n0" ^ remaining_zeros)
end

function get_atoms{N}(frame::Int, types::NTuple{N, Symbol}, atom_types::Vector{Symbol}, num_atoms::Int)
    indices = 1 + (frame - 1) * num_atoms : frame * num_atoms
    findin(sub(atom_types, indices), types) .+ (first(indices) - 1)
end

function rdf{T}(bw::T, C::Matrix{T}, max_radius::T, dr::T)
    num_atoms::Int = size(C, 2)
    num_bins::Int = floor(Int, max_radius / dr)
    int_bins = zeros(UInt64, num_bins)
    bins = Matrix{Float64}(num_bins, 2)

    @inbounds for i in 2:num_atoms
        for j in 1:i-1
            r::T = sqrt(distancesq(bw, C, i, j))
            if r < max_radius
                p::Int = ceil(Int, r / dr) - 1
                int_bins[p] += 1
            end
        end
    end

    @inbounds for i in 1:num_bins
        r_mid = i*dr - dr / 2
        divisor = (num_atoms * 4π * r_mid * r_mid * dr) * (num_atoms / bw ^ 3)
        bins[i, 1] = r_mid
        bins[i, 2] = 2 * int_bins[i] / divisor
    end
    
    return bins
end

function rdf{N, T}(frames::AbstractArray{Int}, types::NTuple{N, Symbol}, box_width::T, max_radius::T, dr::T, num_atoms::Int, atom_coords::AbstractMatrix{T}, atom_types::Vector{Symbol})
    list = map(frames) do f
        indices = get_atoms(f, types, atom_types, num_atoms)
        coords = atom_coords[:, indices]
        
        (box_width, coords, max_radius, dr)
    end
    
    mean(pmap(x -> rdf(x...), list))
end

function density2{T<:AbstractFloat}(file_out::IO, C::AbstractMatrix{T}, n::Int, bw::T = zero(Float64))
    num_atoms = size(C, 2)

    b::T = 2 * maxabs(C)
    println(b)

    if (bw != zero(T))
        bw < b && error("Your box specification is not big enough to contain all atoms.")

        b = bw
    else
        b *= (1 + 1/n)
    end

    println(file_out, "density")
    println(file_out, "$b   "^3, "90.00000   "^3)
    println(file_out, "$n   "^3)

    indices, counts = Int[], Int[]
    @inbounds for i in 1:num_atoms
        xs = section(C[1, i], b, n)
        ys = section(C[2, i], b, n)
        zs = section(C[3, i], b, n)
        index = zs + n*ys + n*n*xs + 1
        push!(indices, index)
        push!(counts, 1)
    end
    ρ = sparsevec(indices, counts)
    full_ρ = full(ρ)

    remaining_zeros = n^3 - length(full_ρ)

    print_joined(file_out, full_ρ, '\n')
    println(file_out, "\n0" ^ remaining_zeros)
end

