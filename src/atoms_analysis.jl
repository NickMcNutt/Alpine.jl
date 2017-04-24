function group{T}(::Type{T}, atoms::Atoms, prop::Symbol)
    props = Dict{T, Atoms}()
    
    for i in atoms.indices
        k = (atoms.props[prop]::Array{T})[i]
        
        if !haskey(props, k)
            props[k] = Atoms(Int[], atoms.props)
        end
        
        push!(props[k].indices, i)
    end
    
    return props
end

function group(atoms::Atoms, prop::Symbol)
    T = eltype(atoms.props[prop])
    group(T, atoms, prop)
end

box_dims(f::Frame) = ntuple(i -> f[:box_max][i] - f[:box_min][i], 3)

@inline function unscale{T}(bl::T, bu::T, bw::T, coord::T)
    coord = bw * coord + bl
    
    if coord >= bu
        coord -= bw
    end

    if coord < bl
        coord = clamp(coord, bl, bu)
    end
    
    return coord
end

function unscale!{T, U <: AbstractMatrix{T}}(bl::Vector{T}, bu::Vector{T}, bw::Vector{T}, coords::U)
    ncols = size(coords, 2)

    for c in 1:ncols, r in 1:3
        @inbounds coords[r, c] = unscale(bl[r], bu[r], bw[r], coords[r, c])
    end

    return coords
end

unscale!(frame::Frame) = (unscale!(box_lower(frame), box_upper(frame), box_widths(frame), frame[:atoms][:coords]) ; frame)

function unwrap!{T}(atoms::Atoms, box_width::T, max_sep::T = 3.0)
    sep_sq = max_sep ^ 2
    bw = box_width
    hbw = bw / 2

    coords = atoms[:coords]

    function distance_sq(i::Int, j::Int)
        x = coords[1, i] - coords[1, j]
        y = coords[2, i] - coords[2, j]
        z = coords[3, i] - coords[3, j]

        x^2 + y^2 + z^2
    end

    num_atoms = length(atoms)

    b = fld(num_atoms, 2)
    path = Int[b]
    inpath = fill(false, num_atoms)
    inpath[b] = true

    prev_path_length = 1
    while true
        for i in 1:num_atoms
            inpath[i] && continue

            for j in path
                if distance_sq(i, j) < sep_sq
                    push!(path, i)
                    inpath[i] = true
                    break
                end
            end
        end

        length(path) == prev_path_length && break
        prev_path_length = length(path)
    end

    while true
        for i in 1:num_atoms
            inpath[i] && continue

            for j in path
                for k in 1:3
                    d = coords[k, i] - coords[k, j]
                    if d > hbw
                        coords[k, i] -= bw
                    elseif d < -hbw
                        coords[k, i] += bw
                    end
                end

                if distance_sq(i, j) < sep_sq
                    push!(path, i)
                    inpath[i] = true
                    break
                end
            end
        end

        length(path) == prev_path_length && break
        prev_path_length = length(path)
    end

    all(inpath) || error("Unable to correctly unwrap all atoms in fragment")
end
