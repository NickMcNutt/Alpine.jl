include("types.jl")

import Base.LinAlg.BLAS: gemm!

function push!(func::Function, neighbor_atoms::Cluster, n::Int)
    num = neighbor_atoms.num
    new_num = num + n

    if new_num > neighbor_atoms.max_num
        error("Too many atoms in neighborhood! Expand yo buffer fool.")
    end

    neighbor_atoms.num = new_num
    r = (num + 1):new_num
    func(sub(neighbor_atoms.indices, r),
         sub(neighbor_atoms.coords, :, r),
         sub(neighbor_atoms.r, r))
end

function push!(neighbor_atoms::Cluster, n::Int)
    num = neighbor_atoms.num
    new_num = num + n

    if new_num > neighbor_atoms.max_num
        error("Too many atoms in neighborhood! Expand yo buffer fool.")
    end

    neighbor_atoms.num = new_num
    r = (num + 1):new_num
    return sub(neighbor_atoms.indices, r), sub(neighbor_atoms.r, r), sub(neighbor_atoms.coords, :, r)
end

function push!(neighbor_atoms::Cluster, atom_index::Int, coords::Vector{Float64})
    push!(neighbor_atoms, 1) do indices, n_coords, r
        indices[1] = atom_index
        copy!(n_coords, coords)
        r[1] = norm3(coords)
    end
end

function push!(neighbor_atoms::Cluster, atom_indices::AbstractVector{Int}, coords::AbstractMatrix{Float64})
    n = length(atom_indices)

    push!(neighbor_atoms, n) do indices, n_coords, r
        copy!(indices, atom_indices)
        copy!(n_coords, coords)
        map!(i -> norm3(sub(coords, :, i)), r, 1:n)
    end
end

nonzeroNeighbors(num_all_neighbors::Dict{Symbol, Int}) = filter((k, v) -> v > 0, num_all_neighbors)

function allocate_clusters(group::Group, cluster_partition_types::Dict{Symbol, Int})
    Dict(Pair{Symbol, Dict{AtomID, Cluster}}[t => Dict(Pair{AtomID, Cluster}[id => Cluster(k) for id in group]) for (t, k) in cluster_partition_types])
end

function kabsch(file_out, file_out_frames, neighbor_clusters, ref_cluster, atom_types)
    nc = neighbor_clusters
    n = length(nc)
    K = zeros(Float64, (n, n))

    num_neighbors = first(nc).num
    A = Matrix{Float64}(3, 3)
    D₁ = Matrix{Float64}(3, num_neighbors)
    D₂ = Matrix{Float64}(3, num_neighbors)
    D_min = Matrix{Float64}(3, num_neighbors)

    println(file_out, round(Int, (n-1)*num_neighbors+ 1))
    #println(file_out, round(Int, n*(n-1)*num_neighbors/2 + 1))
    println(file_out)
    println(file_out, "Li\t0.0\t0.0\t0.0\t1.0")

    #@time @inbounds for i in 1:n
    #@time @inbounds for i in 1:1
    #for (i1, nc1) in nc
    #i1 = ref_i
    nc1 = ref_cluster
    for q in 1:1
        C = nc1.coords
        for nc2 in nc
            dref = 1.0e10

            E = C'
            for p in permutations(1:num_neighbors)
                copy!(D₁, sub(nc2.coords, :, p))
                F = D₁'

                #gemm!('N', 'T', 1.0, C, D₁, 0.0, A)
                #F = svdfact!(A)
                #gemm!('N', 'N', 1.0, F[:U], F[:Vt], 0.0, A)
                #gemm!('N', 'N', 1.0, A, D₁, 0.0, D₂)
                #gemm!('N', 'N', 1.0, F[:U], R, 0.0, A)
                #gemm!('N', 'N', 1.0, A, F[:Vt], 0.0, A)
                #gemm!('N', 'N', 1.0, A, D₁, 0.0, D₁)

                #r1 = rmsd(C, D₁)
                #r2 = rmsd(C, D₂)


                A = E' * F
                AF = svdfact(A)
                V = AF[:U]
                WT = AF[:Vt]
                R1 = WT' * diagm([1.0, 1.0, 1.0]) * V'
                v1 = rmsd(E, F*R1)

                if v1 < dref
                    dref = v1
                    copy!(D_min, (F*R1)')
                end
            end

            for s = 1:num_neighbors
                println(file_out, atom_types[nc2.indices[s]], '\t', D_min[1, s], '\t', D_min[2, s], '\t', D_min[3, s], '\t', s/num_neighbors)
            end

            println(file_out_frames, num_neighbors + 1)
            println(file_out_frames)
            println(file_out_frames, "Li\t0.0\t0.0\t0.0\t1.0")
            for s = 1:num_neighbors
                println(file_out_frames, atom_types[nc2.indices[s]], '\t', D_min[1, s], '\t', D_min[2, s], '\t', D_min[3, s])
            end

        end
        break
    end

    return K
end


function kabsch2(file_out, file_out_frames, neighbor_clusters, ref_cluster, atom_types)
    nc = neighbor_clusters
    n = length(nc)
    K = zeros(Float64, (n, n))

    num_neighbors = first(nc).num
    A = Matrix{Float64}(3, 3)
    D₁ = Matrix{Float64}(3, num_neighbors)
    D₂ = Matrix{Float64}(3, num_neighbors)
    D_min = Matrix{Float64}(3, num_neighbors)
    D2_min = Matrix{Float64}(3, num_neighbors)

    println(file_out, round(Int, 6))
    #println(file_out, round(Int, n*(n-1)*num_neighbors/2 + 1))
    println(file_out)
    println(file_out, "Li\t0.0\t0.0\t0.0\t1.0")

    #@time @inbounds for i in 1:n
    #@time @inbounds for i in 1:1
    #for (i1, nc1) in nc
    #i1 = ref_i
    nc1 = ref_cluster
    for q in 1:1
        C = nc1.coords
        for nc2 in nc
            dref = 1.0e10

            E = C'
            for p in permutations(1:num_neighbors)
                copy!(D₁, sub(nc2.coords, :, p))
                F = D₁'

                #gemm!('N', 'T', 1.0, C, D₁, 0.0, A)
                #F = svdfact!(A)
                #gemm!('N', 'N', 1.0, F[:U], F[:Vt], 0.0, A)
                #gemm!('N', 'N', 1.0, A, D₁, 0.0, D₂)
                #gemm!('N', 'N', 1.0, F[:U], R, 0.0, A)
                #gemm!('N', 'N', 1.0, A, F[:Vt], 0.0, A)
                #gemm!('N', 'N', 1.0, A, D₁, 0.0, D₁)

                #r1 = rmsd(C, D₁)
                #r2 = rmsd(C, D₂)


                A = E' * F
                AF = svdfact(A)
                V = AF[:U]
                WT = AF[:Vt]
                R1 = WT' * diagm([1.0, 1.0, 1.0]) * V'
                v1 = rmsd(E, F*R1)

                if v1 < dref
                    dref = v1
                    copy!(D2_min, (F*R1)')

                    D2_min = (D_min + D2_min) / 2.0
                    v2 = rmsd(D2_min, F*R1)

                    if v2 < dref
                        dref = v2
                        copy!(D_min, D2_min)
                    end
                end
            end

            #for s = 1:num_neighbors
                #println(file_out, atom_types[nc2.indices[s]], '\t', D_min[1, s], '\t', D_min[2, s], '\t', D_min[3, s], '\t', s/num_neighbors)
            #end

            println(file_out_frames, num_neighbors + 1)
            println(file_out_frames)
            println(file_out_frames, "Li\t0.0\t0.0\t0.0\t1.0")
            for s = 1:num_neighbors
                println(file_out_frames, atom_types[nc2.indices[s]], '\t', D_min[1, s], '\t', D_min[2, s], '\t', D_min[3, s])
            end

        end
        nc2 = first(nc)
        for s = 1:num_neighbors
            println(file_out, atom_types[nc2.indices[s]], '\t', D_min[1, s], '\t', D_min[2, s], '\t', D_min[3, s], '\t', s/num_neighbors)
        end
        break
    end

    return K
end

function inner_product_spectrum(neighbor_clusters)
    nc = neighbor_clusters

    samples = zeros(Float64, length(nc) * (first(nc).num)^2)
    t = 1
    for c in nc
        n = c.num
        for i in 1:n
            for j in 1:n
                samples[t] = dot(sub(c.coords, :, i), sub(c.coords, :, j))
                t += 1
            end
        end
    end

    return samples
end


function gram(file_out, file_out_frames, neighbor_clusters::Clusters, atom_types)
    nc = neighbor_clusters
    n = length(nc)
    K = zeros(Float64, (n, n))

    num_neighbors = first(values(neighbor_clusters)).num
    println(num_neighbors)
    #num_neighbors = nc[1, :Li].num
    A = Matrix{Float64}(3, 3)
    D₁ = Matrix{Float64}(3, num_neighbors)
    D₂ = Matrix{Float64}(3, num_neighbors)
    D_min = Matrix{Float64}(3, num_neighbors)
    R = diagm(Float64[1.0, 1.0, -1.0]) 

    println(file_out, round(Int, (n-1)*num_neighbors+ 1))
    #println(file_out, round(Int, n*(n-1)*num_neighbors/2 + 1))
    println(file_out)
    println(file_out, "Li\t0.0\t0.0\t0.0\t1.0")

    #@time @inbounds for i in 1:n
    #@time @inbounds for i in 1:1
    for q in 1:1 # Nonsense
        i = 1 # Reference configuration
        C = nc[i, :all].coords
        #for j in i+1:n
        for j in 1:n
            println("$i  $j")

            D₁ = nc[j, :all].coords
            d_ref_min, M_min = d_ref(C' * C, D₁' * D₁)
            D_min = full(chol(M_min + 1e-13*eye(M_min)))[1:3, :]

            println(file_out_frames, num_neighbors + 1)
            println(file_out_frames)
            println(file_out_frames, "Li\t0.0\t0.0\t1.0")
            for s = 1:num_neighbors
                println(file_out_frames, atom_types[nc[j, :all].indices[s]], '\t', D_min[1, s], '\t', D_min[2, s], '\t', D_min[3, s])
            end

            for s = 1:num_neighbors
                println(file_out, atom_types[nc[j, :all].indices[s]], '\t', D_min[1, s], '\t', D_min[2, s], '\t', D_min[3, s], '\t', s/num_neighbors)
            end

            K[i, j] = K[j, i] = d_ref_min
        end
    end

    return K
end

function mprintln(M)
    map(x -> round(x, 3), M) |> println
end

# Computes the polar decomposition U, H = pd(H) with U unitary and H is Gramian
function pd(A)
    W, Σ, Vt = svd(A)
    U = W*Vt
    H = Vt'*diagm(Σ)*Vt
    return U, H
end

# Finds P that minimizes the Frobenius norm of (A - PB)
function optimalperm_quick(A, B)
    m = size(A, 1)
    n = size(A, 2)

    σ² = 0.1^2
    G = Matrix{Float64}(m, m)

    for i in 1:m, j in 1:m
        r_sq = zero(Float64)
        for k in 1:n
            r_sq += (A[i, k] - B[j, k])^2
        end
        G[i, j] = exp(-r_sq / (2σ²))
    end

    T, V, Ut = svd(G)
    U = Ut'
    P = T * U
    Q = map(round, abs(P))
    return Q
end

# Finds P that minimizes the Frobenius norm of (A - PB)
function optimalperm(A, B)
    m = size(A, 1)
    min_rmsd, min_p = Inf, Vector{Int}(m)
    for p in permutations(1:m)
        r = rmsd(A - sub(B, p, :))
        if r < min_rmsd
            min_rmsd = r
            copy!(min_p, p)
        end
    end
    p = eye(Float64, m)[min_p, :]

    return p
end

function gram3!(clusters_rotated, clusters_joint, clusters_mimicked, ref_id)
    key, clusters = first(clusters_joint)
    A = clusters[ref_id].coords

    for (id, cluster) in clusters
        dref = Inf

        E = cluster.coords'
        D = similar(cluster.coords)
        D_min = zeros(D)
        R_opt = zeros(3, 3)
        for p in permutations(1:cluster.num)
            copy!(D, sub(cluster.coords, :, p))

            S = A * D'
            AF = svdfact(S)
            V = AF[:U]
            Wt = AF[:Vt]
            R = Wt' * V'
            d = rmsd(cluster.coords, R'*D)

            if d < dref
                dref = d
                copy!(R_opt, R')
                copy!(D_min, R_opt*D)
            end
        end

        copy!(clusters_rotated[key][id].coords, D_min)

        for (banana, clusters2) in clusters_mimicked
            gemm!('N', 'N', 1.0, R_opt, clusters2[id].coords, 0.0, clusters_rotated[key2][id].coords)
        end
    end
end

function gram3!_saved(clusters_rotated, clusters_joint, clusters_mimicked, ref_id)
    key, clusters = first(clusters_joint)
    A = clusters[ref_id].coords'

    for (id, cluster) in clusters
        B = cluster.coords'
        B_min = zeros(B)
        
        ε = 0.0
        k = 1
        k_max = 1000
        V = eye(Float64, 3)
        P = optimalperm(A, B)
        #ρ_prev = rmsd(A, P*B)
        ρ_prev = Inf
        ρ = ρ_prev
        while ρ > ε
            k += 1
            P = optimalperm(A, B*V)
            ρ = rmsd(A, P*B*V)

            if abs(ρ_prev - ρ) <= ε || k > k_max
                break
            end

            U, H = pd((P*B)' * A)
            V = U
            ρ_prev = ρ
            ρ = rmsd(A, P*B*V)
        end

        #println("Iter: $k")
        #println(ρ_prev)
        #println(ρ)
        #println()

        copy!(B_min, P*B*V)
        copy!(clusters_rotated[key][id].coords, B_min')

        for (key2, clusters2) in clusters_mimicked
            gemm!('T', 'N', 1.0, V, clusters2[id].coords, 0.0, clusters_rotated[key2][id].coords)
        end
    end
end

function gram2(file_out, file_out_frames, neighbor_clusters, ref_cluster, atom_types)
    num_clusters = length(neighbor_clusters)
    println(length(neighbor_clusters))

    A = ref_cluster.coords'
    n = ref_cluster.num
    n = size(A, 1)

    println(file_out, round(Int, (num_clusters)*n + 1))
    println(file_out)
    println(file_out, "Li\t0.0\t0.0\t0.0\t1.0")

    B_min = zeros(Float64, n, 3)

    for nc in neighbor_clusters
        B = nc.coords'

        ρ_min = Inf
        copy!(B_min, B)
        for p in permutations(1:n)
            PB = B[p, :]
            U, H = pd(PB' * A)
            ρ = rmsd(A, PB*U)

            if ρ < ρ_min
                ρ_min = ρ
                copy!(B_min, PB*U)
            end
        end

        println("ρ: $ρ_min")
        for s in 1:n
            println(file_out, "Li", '\t', B_min[s, 1], '\t', B_min[s, 2], '\t', B_min[s, 3], '\t', s/n)
        end
    end

    f = open("ref.xyz", "w")
    println(f, n + 1)
    println(f)
    println(f, "Li\t0.0\t0.0\t0.0\t1.0")
    for s in 1:n
        println(f, "Li", '\t', A[s, 1], '\t', A[s, 2], '\t', A[s, 3], '\t', s/n)
    end
    close(f)

end

# Constructs a Gram matrix and stores it in M
function gramMatrix(coords)
    M = Array{Float64}(size(coords, 2), size(coords, 2))
    gramMatrix!(M, coords)
end

function gramMatrix!(M, coords)
    for i = 1:size(coords, 2)
        for j = 1:size(coords, 2)
            M[i, j] = dot(coords[:, i], coords[:, j])
        end
    end
end

function gramMatrix!(Σ, k, coords)
    #Σ[:, :, k] = coords' * coords + 1e-13 * eye(Float64, size(coords, 2))
    @inbounds for i = 1:size(coords, 2)
        for j = 1:size(coords, 2)
            Σ[i, j, k] = dot(coords[:, i], coords[:, j])
        end
    end

    @inbounds for i = 1:size(coords, 2)
        Σ[i, i, k] = Σ[i, i, k] + 1e-13
    end

    return Σ
end

function d_ref(Σ₁::AbstractCoords, Σ₂::AbstractCoords)
    n = size(Σ₁, 1)
    i = eye(Σ₁)
    P = copy(i)
    
    d_min = Inf
    Σ_min = copy(P)
    for p in permutations(1:n)
        P = i[p, :]
        d = norm(Σ₁ - P * Σ₂ * P')
        if d < d_min
            d_min = d
            Σ_min = P * Σ₂ * P'
        end
    end

    return (d_min, Σ_min)
end
    
function d_ref1(coords1, coords2)
    num_coords = size(coords1, 2)

    gram1 = zeros(Float64, (num_coords, num_coords))
    gram2 = zeros(Float64, (num_coords, num_coords))

    gramMatrix!(gram1, coords1)
    indices_order = [1:num_coords]
    for i = 1:factorial(num_coords)
        gramMatrix!(gram2, coords2[:, nthperm(indices_order, i)])
        chol 
        println(norm(gram1 - gram2))
    end
end
