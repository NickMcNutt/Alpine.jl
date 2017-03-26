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

function generate_bonds{T}(atoms::Atoms, bond_types::Dict{Tuple{Symbol, Symbol}, Tuple{T, T}})
    vertices = Set{Atom}()
    edges = Vector{Tuple{Atom, Atom}}()
    
    coords = atoms[:coords]
    types = atoms[:type]
    
    n = size(coords, 2)::Int
    
    for i in 2:n, j in 1:i-1
        type_i, type_j = types[i]::Symbol, types[j]::Symbol
        
        dist_min, dist_max = get(bond_types, (type_i, type_j), get(bond_types, (type_j, type_i), (Inf, -Inf)))
        if dist_min != Inf
            if dist_min^2 <= distancesq(Inf, coords, i, j) <= dist_max^2
                push!(vertices, atoms[i], atoms[j])
                push!(edges, (atoms[i], atoms[j]))
            end
        end
    end
    
    return Graph(collect(vertices), edges)
end

function unique_upto_reversal{T}(v::Vector{Vector{T}})
    unique = Set{Vector{T}}()
    for a in v
        if a ∉ unique && reverse(a) ∉ unique
            push!(unique, a)
        end
    end
    collect(unique)
end

function unique_neighborhoods(path_table::Vector{Dict{Vector{Atom}, Vector{Atom}}})
    neighborhoods = Dict{Symbol, Dict{Vector{Symbol}, Vector{Atom}}}()
    for (atom, neighbors) in path_table[1]
        atom_type = atom[1][:type]
        neighborhood = sort(Atoms(neighbors)[:type])
        push!(get!(get!(neighborhoods, atom_type, Dict{Vector{Symbol}, Vector{Atom}}()), neighborhood, Symbol[]), atom[1])
    end
    neighborhoods
end

function neighborhood_type_ids(neighborhoods)
    Dict(
        atom_type => Dict(
            first(neighborhood_type) => i
        for (i, neighborhood_type) in enumerate(neighborhood_types))
    for (atom_type, neighborhood_types) in neighborhoods)
end
