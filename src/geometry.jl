# Generate a rotation matrix for a given axis of rotation v and angle θ
function rotation{T<:AbstractFloat}(v::AbstractVector{T}, θ::T)
    @inbounds begin
        K = Matrix{T}(3, 3)
        m = mag(v)
        v1 = v[1] / m ; v2 = v[2] / m ; v3 = v[3] / m

        K[1, 1] =  0  ; K[1, 2] = -v3 ; K[1, 3] =  v2
        K[2, 1] =  v3 ; K[2, 2] =  0  ; K[2, 3] = -v1
        K[3, 1] = -v2 ; K[3, 2] =  v1 ; K[3, 3] =  0

        R = I + sin(θ)*K + (1 - cos(θ))*K^2

        return R
    end
end

# Generate the Rodrigues rotation matrix R that rotates vector v1 to vector v2
function rotation{T<:AbstractFloat}(v1::AbstractVector{T}, v2::AbstractVector{T})
    axis = cross(v1, v2)
    θ = -acos(dot(v1, v2) / (mag(v1) * mag(v2)))
    rotation(axis, θ)
end

# Given two vectors v1 and v2, generate the rotation matrix that
# aligns v1 with the x-axis and makes v2 orthogonal to the z-axis
function rotation{T<:AbstractFloat}(v1::AbstractVector{T}, v2::AbstractVector{T})
    xaxis = T[1, 0, 0]
    v = Vector{T}(3)

    R1 = rotation(v1, xaxis)
    A_mul_B!(v, R1, v2)
    θ = acos(v[3] / sqrt(v[2]^2 + v[3]^2))
    R2 = rotation(xaxis, θ)
    A_mul_B!(R2, R2, R1)
end

# Generate a random rotation matrix
randomrotation() = rotation(randn(3), 2π * rand())
