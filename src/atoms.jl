function group(f::Function, indices::AbstractVector{Int}, atom_property, categories)
    groups = Dict([c => Indices() for c in categories])

    for i in indices
        for c in categories
            if f(atom_property[i], c)
                push!(groups[c], i)
                break
            end
        end
    end

    return groups
end

function coalesce(indices::Indices; props...)
    new_props = map(props) do prop
        dims = ndims(prop[2])
        dims == 1 && return (prop[1], prop[2][indices])
        dims == 2 && return (prop[1], prop[2][:, indices])
        return nothing
    end
    push!(new_props, (:indices => indices))
    Dict(new_props)
end
