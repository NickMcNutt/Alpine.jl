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

function rdf{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, max_radius::T, dr::T)
    num_atoms::Int = size(C, 2)
    num_bins::Int = floor(Int, max_radius / dr)
    int_bins = zeros(UInt64, num_bins)
    bins = Matrix{Float64}(num_bins, 2)

    @inbounds for i in 2:num_atoms, j in 1:i-1
        r = sqrt(distancesq(bw, C, i, j))
        r < max_radius && (int_bins[fld(r, dr) + 1] += 1)
    end

    @inbounds for i in 1:num_bins
        r_mid = i*dr - dr / 2
        divisor = (num_atoms * 4π * r_mid * r_mid * dr) * (num_atoms / bw ^ 3)
        bins[i, 1] = r_mid
        bins[i, 2] = 2 * int_bins[i] / divisor
    end
    
    return bins
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

function rdf{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, max_radius::T, dr::T)
    num_atoms::Int = size(C, 2)
    num_bins::Int = floor(Int, max_radius / dr)
    int_bins = zeros(UInt64, num_bins)
    bins = Matrix{Float64}(num_bins, 2)

    @inbounds for i in 2:num_atoms, j in 1:i-1
        r = sqrt(distancesq(bw, C, i, j))
        r < max_radius && (int_bins[fld(r, dr) + 1] += 1)
    end

    @inbounds for i in 1:num_bins
        r_mid = i*dr - dr / 2
        divisor = (num_atoms * 4π * r_mid * r_mid * dr) * (num_atoms / bw ^ 3)
        bins[i, 1] = r_mid
        bins[i, 2] = 2 * int_bins[i] / divisor
    end
    
    return bins
end
