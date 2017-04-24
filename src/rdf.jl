box_upper(frame::Alpine.Frame) = Float64[frame[:box_max]...]
box_lower(frame::Alpine.Frame) = Float64[frame[:box_min]...]
box_widths(frame::Alpine.Frame) = box_upper(frame) - box_lower(frame)

function distance_sq{T}(xbw::T, ybw::T, zbw::T, coords::Matrix{T}, i1::Int, i2::Int)
    xhbw = xbw / 2
    yhbw = ybw / 2
    zhbw = zbw / 2
    
    @inbounds xd = coords[1, i1] - coords[1, i2]
    @inbounds yd = coords[2, i1] - coords[2, i2]
    @inbounds zd = coords[3, i1] - coords[3, i2]
    
    if xd < -xhbw xd += xbw elseif xd > xhbw xd -= xbw end
    if yd < -yhbw yd += ybw elseif yd > yhbw yd -= ybw end
    if zd < -zhbw zd += zbw elseif zd > zhbw zd -= zbw end

    return xd*xd + yd*yd + zd*zd
end

function wrap(bw, z)
    hbw = bw / 2
    
    if z < -hbw
        z += bw
    elseif z > hbw
        z -= bw
    end
    
    return z
end

function distance_cells(box_width, cell1_center, cell1_width, cell2_center, cell2_width)
    bw = box_width
    hcw1, hcw2 = cell1_width / 2, cell2_width / 2
    cc = cell2_center - cell1_center
    
    d1 = abs(wrap(bw, cc + hcw2 - hcw1))
    d2 = abs(wrap(bw, cc + hcw2 + hcw1))
    d3 = abs(wrap(bw, cc - hcw2 - hcw1))
    d4 = abs(wrap(bw, cc - hcw2 + hcw1))
    
    return min(d1, d2, d3, d4), max(d1, d2, d3, d4)
end

function neighbor_cell_offsets{T}(frame::Alpine.Frame, r_cutoff::T, num_cells::Int)
    bw = box_widths(frame)
    cw = bw / num_cells
    
    xbw, ybw, zbw = bw[1], bw[2], bw[3]
    xcw, ycw, zcw = cw[1], cw[2], cw[3]
    
    neighbor_offsets = Vector{Int}[]
    
    for ix in 0:num_cells-1, iy in 0:num_cells-1, iz in 0:num_cells-1
        min_x, max_x = distance_cells(xbw, T(0), xcw, ix*xcw, xcw)
        min_y, max_y = distance_cells(ybw, T(0), ycw, iy*ycw, ycw)
        min_z, max_z = distance_cells(zbw, T(0), zcw, iz*zcw, zcw)

        if min_x^2 + min_y^2 + min_z^2 <= r_cutoff^2
            push!(neighbor_offsets, [ix, iy, iz])
        end
    end
    
    return neighbor_offsets
end

cell_index{T}(xlo::T, xcw::T, x::T) = floor(Int, (x - xlo) / xcw)::Int + 1

function sort_atoms_into_cells(frame::Alpine.Frame, num_cells::Int)
    num_atoms = length(frame[:atoms])::Int
    coords = collect(frame[:atoms][:coords])::Matrix{Float64}
    
    xlo::Float64, ylo::Float64, zlo::Float64 = box_lower(frame)
    xcw::Float64, ycw::Float64, zcw::Float64 = box_widths(frame) / num_cells
    
    cells = [Int[] for ix in 1:num_cells, iy in 1:num_cells, iz in 1:num_cells]

    for i in 1:num_atoms
        @inbounds ix = cell_index(xlo, xcw, coords[1, i])
        @inbounds iy = cell_index(ylo, ycw, coords[2, i])
        @inbounds iz = cell_index(zlo, zcw, coords[3, i])
        
        @inbounds push!(cells[ix, iy, iz], i)
    end
    
    return cells
end

function ndf!{T, U <: AbstractVector{Int}}(bins::Vector{UInt64}, xbw::T, ybw::T, zbw::T, coords::Matrix{T}, indices1::U, indices2::U, r_cutoff_sq::T, Δr::T)
    for i1 in indices1, i2 in indices2
        i1 == i2 && continue
        r_sq = distance_sq(xbw, ybw, zbw, coords, i1, i2)
        if r_sq < r_cutoff_sq
            r = sqrt(r_sq)
            b = floor(Int, r / Δr) + 1
            @inbounds bins[b] += 1
        end
    end
    
    return bins
end

function func_ndf{T}(frame::Alpine.Frame, r_cutoff::T, Δr::T, num_cells::Int)
    r_cutoff_sq = r_cutoff^2
    
    coords = collect(frame[:atoms][:coords])::Matrix{T}
    xbw::T, ybw::T, zbw::T = box_widths(frame)
    
    num_bins = ceil(Int, r_cutoff / Δr)
    dim_cells = (num_cells, num_cells, num_cells)
    neighbor_offsets = neighbor_cell_offsets(frame, r_cutoff, num_cells)
    
    return function (cells1::Array{Vector{Int}, 3}, cells2::Array{Vector{Int}, 3})
        bins = zeros(UInt64, num_bins)
        
        for i1x in 1:num_cells, i1y in 1:num_cells, i1z in 1:num_cells
            i1 = sub2ind(dim_cells, i1x, i1y, i1z)

            for offset in neighbor_offsets
                @inbounds i2x = mod1(i1x + offset[1], num_cells)
                @inbounds i2y = mod1(i1y + offset[2], num_cells)
                @inbounds i2z = mod1(i1z + offset[3], num_cells)

                i2 = sub2ind(dim_cells, i2x, i2y, i2z)

                if i1 >= i2
                    @inbounds indices1 = cells1[i1x, i1y, i1z]
                    @inbounds indices2 = cells2[i2x, i2y, i2z]

                    ndf!(bins, xbw, ybw, zbw, coords, indices1, indices2, r_cutoff_sq, Δr)
                end
            end
        end
        
        return bins
    end
end

function ndf_to_rdf{T}(bins_ints::Vector{UInt64}, ρ::T, Δr::T)
    num_bins = length(bins_ints)
    bins = zeros(T, num_bins)
    
    for i in 1:num_bins
        r_mid = i*Δr - Δr/2
        bins[i] = bins_ints[i] / (4π * r_mid^2 * Δr * ρ)
    end
    
    return bins
end

function split_by_types(frame::Alpine.Frame)
    types = group(frame[:atoms], :type)
    Dict(k => Alpine.Frame(
        :atoms => v,
        :timestep => frame[:timestep],
        :box_min => frame[:box_min],
        :box_max => frame[:box_max]
    ) for (k, v) in types)
end

function rdf_components{T}(frames::Vector{Alpine.Frame}, component_pairs::Vector{Vector{Symbol}}, r_cutoff::T, Δr::T, num_cells::Int)
    num_component_pairs = length(component_pairs)
    components = sort(unique(vcat(component_pairs...)))
    num_bins = ceil(Int, r_cutoff / Δr)
    
    bins_total = zeros(T, (num_bins, num_component_pairs))
    for (f, frame) in enumerate(frames)
        bins = zeros(T, (num_bins, num_component_pairs))
        
        xbw, ybw, zbw = box_widths(frame)
        volume = xbw * ybw * zbw
        
        frame_components = split_by_types(frame)
        num_atoms = Dict(component => length(frame_components[component][:atoms])::Int for component in components)
        cells = Dict(component => sort_atoms_into_cells(frame_components[component], num_cells) for component in components)
        
        ndf = func_ndf(frame, r_cutoff, Δr, num_cells)
        
        for (i, component_pair) in enumerate(component_pairs)
            c1, c2 = component_pair[1], component_pair[2]
            bins_ints = ndf(cells[c1], cells[c2])
            ρ = (num_atoms[c1] * num_atoms[c2]) / volume
            bins[:, i] .= ndf_to_rdf(bins_ints, ρ, Δr)
        end
        
        bins_total += bins

        println(f)
        flush(STDOUT)
    end
    
    rdfs = Dict(join(component_pair, '-') => bins_total[:, i] for (i, component_pair) in enumerate(component_pairs))
    
    return rdfs
end

function rdf_components_all{T}(frames::Vector{Alpine.Frame}, r_cutoff::T, Δr::T, num_cells::Int)
    components = sort(unique(frames[1][:atoms][:type]))::Vector{Symbol}
    component_pairs = unique([sort([c1, c2]) for c1 in components, c2 in components])
    
    return rdf_components(frames, component_pairs, r_cutoff, Δr, num_cells)
end
