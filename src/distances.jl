# Compute total squared displacement of a set of coords
function sd{T<:AbstractFloat}(C::AbstractArray{T})
    t = zero(T)
    @inbounds @simd for i in eachindex(C)
        t += C[i] ^ 2
    end
    t
end

# Compute total squared displacement between two sets of coords
function sd{T<:AbstractFloat}(C1::AbstractArray{T}, C2::AbstractArray{T})
    t = zero(T)
    @inbounds @simd for i in eachindex(C1)
        t += (C1[i] - C2[i])^2
    end
    t
end

# Compute total squared displacement between two sets of coords in a periodic box
function sd{T<:AbstractFloat}(bw::T, C1::AbstractArray{T}, C2::AbstractArray{T})
    t = zero(T)
    @inbounds @simd for i in size(C1, 2)
        t += distancesq(bw, C1[1, i], C1[2, i], C1[3, i], C2[1, i], C2[2, i], C2[3, i])
    end
    t
end

# Compute mean squared displacements
msd{T<:AbstractFloat}(C::AbstractArray{T}) = sd(C) / length(C)
msd{T<:AbstractFloat}(C1::AbstractArray{T}, C2::AbstractArray{T}) = sd(C1, C2) / length(C1)
msd{T<:AbstractFloat}(bw::T, C1::AbstractArray{T}, C2::AbstractArray{T}) = sd(bw, C1, C2) / length(C1)

# Compute root mean squared deviations
rmsd{T<:AbstractFloat}(C::AbstractArray{T}) = sqrt(msd(C))
rmsd{T<:AbstractFloat}(C1::AbstractArray{T}, C2::AbstractArray{T}) = sqrt(msd(C1, C2))
rmsd{T<:AbstractFloat}(bw::T, C1::AbstractArray{T}, C2::AbstractArray{T}) = sqrt(msd(bw, C1, C2))

# Compute the magnitude of a vector
mag{T<:AbstractFloat}(C::AbstractVector{T}) = sqrt(sd(C))

# Wrap a set of coords into a periodic box
function wrapcoords!{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T})
    hbw = bw / 2
    @inbounds @simd for i in 1:size(C, 2)
        if C[1] < -hbw C[1] += bw elseif C[1] > hbw C[1] -= bw end
        if C[2] < -hbw C[2] += bw elseif C[2] > hbw C[2] -= bw end
        if C[3] < -hbw C[3] += bw elseif C[3] > hbw C[3] -= bw end
    end
end

# Unwrap each coord in C into its image in IC
function unwrapcoords!{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, IC::AbstractMatrix{T})
    hbw = bw / 2
    @inbounds @simd for i in 1:size(C, 2)
        if (C[1, i] - IC[1, i]) < hbw C[1, i] += bw elseif C[1, i] - IC[1, i] > hbw C[1, i] -= bw end
        if (C[2, i] - IC[2, i]) < hbw C[2, i] += bw elseif C[2, i] - IC[2, i] > hbw C[2, i] -= bw end
        if (C[3, i] - IC[3, i]) < hbw C[3, i] += bw elseif C[3, i] - IC[3, i] > hbw C[3, i] -= bw end
    end
end

# Compute the squared distance between two coords
function distancesq{T<:AbstractFloat}(bw::T, x1::T, y1::T, z1::T, x2::T, y2::T, z2::T)
    hbw = bw / 2
    xd = x1 - x2
    yd = y1 - y2
    zd = z1 - z2
    if xd < -hbw xd += bw elseif xd > hbw xd -= bw end
    if yd < -hbw yd += bw elseif yd > hbw yd -= bw end
    if zd < -hbw zd += bw elseif zd > hbw zd -= bw end

    return xd*xd + yd*yd + zd*zd
end

distancesq{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, i::Int, j::Int) = distancesq(bw, C[1, i], C[2, i], C[3, i], C[1, j], C[2, j], C[3, j])
distancesq{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, i::Int, x::T, y::T, z::T) = distancesq(bw, C[1, i], C[2, i], C[3, i], x, y, z)
distancesq{T<:AbstractFloat}(bw::T, C1::AbstractMatrix{T}, C2::AbstractMatrix{T}, i1::Int, i2::Int) = distancesq(bw, C1[1, i1], C1[2, i1], C1[3, i1], C2[1, i2], C2[2, i2], C2[3, i2])

# Compute the distances between a set of coords and a particular coord
function distances!{T<:AbstractFloat}(bw::T, D::AbstractVector{T}, C::AbstractMatrix{T}, x::T, y::T, z::T)
    hbw = bw / 2
    @inbounds @simd for i in 1:size(D, 2)
        xd = C[1, i] - x
        yd = C[2, i] - y
        zd = C[3, i] - z
        if xd < -hbw xd += bw elseif xd > hbw xd -= bw end
        if yd < -hbw yd += bw elseif yd > hbw yd -= bw end
        if zd < -hbw zd += bw elseif zd > hbw zd -= bw end
        D[i] = sqrt(xd^2 + yd^2 + zd^2)
    end
end

function distances!{T<:AbstractFloat}(bw::T, D::AbstractArray{T, 3}, ci::Int, C::AbstractMatrix{T}, indices::AbstractMatrix{Int}, i::Int)
    hbw = bw / 2
    @inbounds begin
        x = C[1, i]
        y = C[2, i]
        z = C[3, i]
        @simd for j in 1:size(indices, 1)
            xd = C[1, indices[j, ci]] - x
            yd = C[2, indices[j, ci]] - y
            zd = C[3, indices[j, ci]] - z
            if xd < -hbw xd += bw elseif xd > hbw xd -= bw end
            if yd < -hbw yd += bw elseif yd > hbw yd -= bw end
            if zd < -hbw zd += bw elseif zd > hbw zd -= bw end
            D[1, j, ci] = xd
            D[2, j, ci] = yd
            D[3, j, ci] = zd
        end
    end
end
