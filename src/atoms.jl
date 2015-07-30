using DataFrames

include("types.jl")

function getIndicesByAtomType(atom_type::Symbol, atom_types::AbstractVector{Symbol})
    indices = Group()

    for i = 1:length(atom_types)
        if atom_types[i] == atom_type
            push!(indices, i)
        end
    end

    return indices
end

function getIndices(f::Function, atom_types::AbstractVector{Symbol}, atom_coords::AbstractMatrix{Float64}, atom_misc::AbstractMatrix{Any})
    indices = Int64[]
    
    for i = 1:length(atom_types)
        if f(atom_types[i], atom_coords[:, i], atom_misc[:, i])
            push!(indices, i)
        end
    end

    return indices
end

function createGroupsFromAtomTypes(atom_type_names::Set{Symbol}, atom_types::AbstractVector{Symbol})
    Dict{Symbol, Group}([t => getIndicesByAtomType(t, atom_types) for t in atom_type_names])
end

function getGroupCoords(group_indices::Dict{Symbol, Vector{Int}}, atom_coords::AbstractMatrix{Float64})
    Dict{Symbol, Matrix{Float64}}([t => atom_coords[:, group_indices[t]] for t in keys(group_indices)])
end

filterGroup(func::Function, group::Group) = filter(func, group)

filterGroup(func::Function, group::Group, categories::AbstractVector) = filter(x -> categories[func(x)], group)

function groupBy(func::Function, group::Group, atom_property, categories)
    groups = Dict([c => Group() for c in categories])

    for id in group
        for c in categories
            if func(atom_property[id], c)
                push!(groups[c], id)
                break
            end
        end
    end

    return groups
end

function coalesce(group::Group; props...)
    atom_ids = collect(group)
    new_props = map(props) do prop
        dims = ndims(prop[2])
        dims == 1 && return (prop[1], prop[2][atom_ids])
        dims == 2 && return (prop[1], prop[2][:, atom_ids])
        return nothing
    end
    push!(new_props, (:id, atom_ids))
    Dict(new_props)
end
