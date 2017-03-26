function rename_param(param)
    mappings = Dict(
        "atom" => "name",
        "r1" => "r₁",
        "theta0" => "θ₀",
        "x1" => "x₁",
        "d1" => "d₁",
        "zeta" => "ζ",
        "z1" => "z₁",
        "vi" => "vᵢ",
        "uj" => "uⱼ",
        "xi" => "xᵢ"
    )

    get(mappings, param, param)
end

atom_type_number(atom_type) = Symbol(String(atom_type[:name])[3:min(3, end)])
atom_type_name(atom_type) = Symbol(String(atom_type[:name])[1:2])

@generated function atom_type_params(uff::UFF, atom_type::Symbol)
    filename_uff_params = joinpath(Pkg.dir("Alpine"), "src/potentials/uff", "uff.txt")
    header, params = readdlm(filename_uff_params, header = true)
    params = rename_param.(params)
    atom_types = Dict(
        Symbol(header[i, 1]) => Dict(
            Symbol(param) => header[i, j] for (j, param) in enumerate(params)
        ) for i in 1:size(header, 1)
    )

    :($atom_types[atom_type])
end

# van der Waals
function vdw(uff::UFF, i, j)
    x_ij = √(i[:x₁] * j[:x₁])
    d_ij = √(i[:d₁] * j[:d₁])
    ϵ = d_ij
    σ = x_ij * 0.5^(1/6)
    
    # d = d_ij, x = x_ij
    return ϵ, σ
end

# Bond stretching
function stretching(uff::UFF, i, j, n)
    λ = 0.1332
    r_bo = -λ * (i[:r₁] + j[:r₁]) * log(n)
    r_en = i[:r₁] * j[:r₁] * (√(i[:xᵢ]) - √(j[:xᵢ]))^2 / (i[:xᵢ] * i[:r₁] + j[:xᵢ] * j[:r₁])
    r_ij = i[:r₁] + j[:r₁] + r_bo - r_en
    
    g = 332.06
    k_ij = (2 * g * i[:z₁] * j[:z₁]) / (r_ij^3)
    k_ij /= 2
    
    r, k = r_ij, k_ij
    return r, k
end

# Bond bending (θ₀ [in degrees] is optional)
function bending(uff::UFF, i, j, k, n_ij, n_jk, θ₀ = j[:θ₀])
    c₂ = 1 / (4 * sind(θ₀)^2)
    c₁ = -4 * c₂ * cosd(θ₀)
    c₀ = c₂ * (2 * cosd(θ₀)^2 + 1)
    
    r_ij = stretching(i, j, n_ij)[1]
    r_jk = stretching(j, k, n_jk)[1]
    r_ik = √(r_ij^2 + r_jk^2 - 2*r_ij*r_jk*cosd(θ₀))
    β = 664.12 / (r_ij * r_jk)
    
    k_ijk = β * ((i[:z₁] * k[:z₁]) / r_ik^5) * r_ij * r_jk * (3 * r_ij * r_jk * (1 - cosd(θ₀)^2) - (r_ik^2)*cosd(θ₀))
    
    # c₀ = c₀, c₁ = c₁, c₂ = c₂, k = k_ijk
    return (k_ijk / (2 * sind(θ₀)^2)), θ₀
end

# Bond torsion
function torsion(uff::UFF, i, j, k, ℓ, bond_order_jk, num_torsions)
    hybr_i, hybr_j, hybr_k, hybr_ℓ = atom_type_number.((i, j, k, ℓ))
    name_j, name_k = atom_type_name.((j, k))
    
    hybrid = (hbyr_j, hybr_k)
    
    local n, ϕ₀, v_ijkℓ
    
    if hybrid == (:3, :3)
        n = 3
        ϕ₀ = 180
        v_ijkℓ = √(j[:vᵢ] * k[:vᵢ])
    elseif hybrid in ((:2, :3), (:3, :2), (:3, :R), (:R, :3))
        n = 6
        ϕ₀ = 0
        v_ijkℓ = 1
    elseif hybrid in ((:2, :2), (:2, :R), (:R, :2), (:R, :R))
        n = 2
        ϕ₀ = 180
        v_ijkℓ = 5 * √(j[:uⱼ] * k[:uⱼ]) * (1 + 4.18 * log(bond_order_jk)) # Official paper says 4.18; 2nd paper says 4.8
    else
        n = 0
        ϕ₀ = 0
        v_ijkℓ = 0
    end
    
    group16_elements = (:O_, :S_, :Se, :Te, :Po, :Lv)
    if name_j ∈ group16_elements || name_k ∈ group16_elements
        if hybrid == (:3, :3)
            n = 2
            ϕ₀ = 90
            v_j = name_j == :O_ ? 2 : 6.8
            v_k = name_k == :O_ ? 2 : 6.8
            v_ijkℓ = √(v_j * v_k)
        elseif hybrid in ((:2, :3), (:3, :2), (:R, :3), (:3, :R))
            n = 2
            ϕ₀ = 90
            v_ijkℓ = 5 * √(j[:uⱼ] * k[:uⱼ]) # Eq. 17 w/ bond-order = 1
            bond_order_jk != 1 && println("Bond order ≂̸ 1")
        end
        
        if (hybrid ∈ ((:2, :3), (:R, :3)) && hybr_i == (:2, :R)) || (hybrid ∈ ((:3, :2), (:3, :R)) && hybr_ℓ == (:2, :R))
            n = 3
            ϕ₀ = 180
            v_ijkℓ = 2
        end
    end
    
    v_ijkℓ /= num_torsions
    
    # k = v_ijkℓ / 2, ϕ₀ = ϕ₀, n = n
    
    # return k, ϕₒ, n
    return v_ijkℓ / 2, -round(cosd(n * ϕ₀), 0), n
end

# Bond inversion — only supports ω₀ = 0 for now
function inversion(uff::UFF, i, j, k, ℓ)
    hybr_i, hybr_j, hybr_k, hybr_ℓ = atom_type_number.((i, j, k, ℓ))
    name_i, name_j, name_k = atom_type_name.((i, j, k))
    
    # a₂ = 1 / (cosd(γ₀)^2 + 2*sind(γ₀) - 2)
    # a₁ = 4 * a₂ * sind(γ₀)
    # a₀ = a₂ * (2*cosd(γ₀)^2 - 3)
    
    if name_i == :C_ && hybr_i == (:2, :R)
        c₀ = 1
        c₁ = -1
        c₂ = 0
        
        if :O_2 in (name_j, name_k, name_ℓ)
            k_ijkℓ = 50
        else
            k_ijkℓ = 6
        end
        
        # k = k_ijkℓ, c₀ = c₀, c₁ = c₁, c₂ = c₂
        return k_ijkℓ
    else
        return 0
    end
end
