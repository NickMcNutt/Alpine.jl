type Graph{T}
    V::Vector{T}
    E::Vector{Tuple{T, T}}
end

vertices(G::Graph) = G.V
edges(G::Graph) = G.E

function acyclic_path_table{T}(G::Graph{T}, max_path_length::Int = 20)
    paths = [Dict{Vector{T}, Vector{T}}()]
    
    for E in edges(G)
        push!(get!(paths[1], [E[1]], []), E[2])
        push!(get!(paths[1], [E[2]], []), E[1])
    end
    
    for i in 1:max_path_length
        push!(paths, Dict{Vector{T}, Vector{T}}())

        for (stem, leaves) in paths[i]
            for leaf in leaves
                new_stem = [stem ; leaf]
                paths[i + 1][new_stem] = T[]
                prefix, suffix = new_stem[1], new_stem[2:end]
                for p in paths[i][suffix]
                    p != prefix && push!(paths[i + 1][new_stem], p)
                end
            end
        end

        isempty(paths[i + 1]) && break
    end
    
    return paths
end

reversal_invariant(v) = hash(v) ⊻ hash(reverse(v))

unique_paths{T}(paths::Vector{Dict{Vector{T}, Vector{T}}}) = (p -> unique(reversal_invariant, p)).(collect.(keys.(paths)))

acyclic_paths{T}(G::Graph{T}, max_path_length::Int = 20) = unique_paths(acyclic_path_table(G, max_path_length))
