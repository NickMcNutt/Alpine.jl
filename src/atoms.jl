import Base: length, start, next, done, eachindex, get, getindex, setindex!, show, display, keys, haskey, merge, similar, hash, ==, isequal, endof, merge

# Types
FrameProps = Dict{Symbol, Any}

immutable Frame <: Associative{Symbol, Any}
    props::FrameProps
end

AtomProps = Dict{Symbol, Array}

immutable Atoms <: Associative{Symbol, Array}
    indices::Indices
    props::AtomProps
end

length(atoms::Atoms) = length(atoms.indices)

immutable Atom <: Associative{Symbol, Array}
    index::Int
    props::AtomProps
end

length(atom::Atom) = 1


# Constructors
Frame(props::Pair...) = Frame(FrameProps(props))

Atoms(num_atoms::Int, props::AtomProps) = Atoms(Indices(1:num_atoms), props)
Atoms(num_atoms::Int, props::Pair...) = Atoms(num_atoms, AtomProps(props))
Atoms{T <: AbstractVector{Int}}(indices::T, props::AtomProps) = Atoms(Indices(indices), props)
Atoms(atoms1::Atoms, atoms2::Atoms) = Atoms(atoms1.indices, atoms2.props)
Atom(atom1::Atom, atoms2::Atoms) = Atom(atom1.index, atoms2.props)
function Atoms(atoms_list::Vector{Atom})
    indices = Indices()
    for item in atoms_list
        push!(indices, item.index)
    end

    # Assumes that all atoms_list.props is same props structure.  Find fast way to check this for sure.
    Atoms(indices, first(atoms_list).props)
end

# Frame accessors
length(f::Frame) = length(f.props)
haskey(f::Frame, key) = haskey(f.props, key)
start(f::Frame) = start(f.props)
next(f::Frame, args...) = next(f.props, args...)
done(f::Frame, args...) = done(f.props, args...)
get(f::Frame, args...) = get(f.props, args...)
setindex!(f::Frame, args...) = setindex!(f.props, args...)

# Atoms accessors
keys(atoms::Atoms) = keys(atoms.props)
haskey(atoms::Atoms, key) = haskey(atoms.props, key)
start(atoms::Atoms) = start(atoms.indices)
function next(atoms::Atoms, state)
    index, state = next(atoms.indices, state)
    Atom(index, atoms.props), state
end
done(atoms::Atoms, state) = done(atoms.indices, state)
eachindex(atoms::Atoms) = 1:length(atoms)
endof(atoms::Atoms) = endof(atoms.indices)
function getindex(atoms::Atoms, index::Int)
    @inbounds i = atoms.indices[index]
    Atom(i, atoms.props)
end
getindex{T <: AbstractVector{Int}}(atoms::Atoms, indices::T) = Atoms(atoms.indices[indices], atoms.props)

function getindex(atoms::Atoms, prop::Symbol)
    i = atoms.indices
    p = atoms.props[prop]

    if ndims(p) == 1
        @inbounds v = view(p, i)
        return v
    elseif ndims(p) == 2
        @inbounds v = view(p, :, i)
        return v
    elseif ndims(p) == 3
        @inbounds v = view(p, :, :, i)
    else
        return 0
    end
end

getindex(atoms::Atoms, atom::Atom) = Atom(atom, atoms)

function setindex!(atoms::Atoms, a::Vector, prop::Symbol)
    atoms.props[prop] = a
end

function merge(atoms::Atoms...)
    num_atoms = sum(length, atoms)
    prop_types = intersect(map(keys, atoms)...)
    props = Dict(
        prop => cat(
            ndims(first(atoms).props[prop]),
            map(a -> a[prop], atoms)...
        ) for prop in prop_types
    )

    Atoms(num_atoms, props)
end

function add_property!{T}(atoms::Atoms, prop::Symbol, default_value::T)
    num_atoms = last(size(first(values(atoms.props))))
    atoms.props[prop] = fill(default_value, num_atoms)
end

# Atom accessors
function getindex(atom::Atom, prop::Symbol)
    i = atom.index
    p = atom.props[prop]

    if ndims(p) == 1
        @inbounds v = p[i]
        return v
    elseif ndims(p) == 2
        @inbounds v = p[:, i]
        return v
    elseif ndims(p) == 3
        @inbounds v = p[:, :, i]
    else
        return 0
    end
end

function setindex!{T}(atom::Atom, value::T, prop::Symbol)
    @inbounds atom.props[prop][atom.index] = value
end

hash(atom::Atom) = hash(object_id(atom.props)) ⊻ hash(atom.index)
hash(atom::Atom, h::UInt) = hash(atom) ⊻ hash(h)
==(atom1::Atom, atom2::Atom) = hash(atom1) == hash(atom2)
isequal(atom1::Atom, atom2::Atom) = atom1 == atom2

# Display

function show(io::IO, frame::Frame)
    props = frame.props
    show_special = (:num_atoms, :num_bonds, :timestep, :atoms, :bonds)
    haskey(props, :atoms) && print(io, length(props[:atoms]), " atoms  ")
    haskey(props, :bonds) && print(io, size(props[:bonds], 2), " bonds  ")
    haskey(props, :timestep) && print(io, "t = ", props[:timestep], "  ")
    for (k, v) in props
        k ∈ show_special || print(io, k, " = ", v, "  ") 
    end
end

display(frame::Frame) = show(frame)

function show(io::IO, atoms::Atoms)
    print(io, "$(length(atoms)) atoms")
    print(io, " of types $(tuple(unique(atoms[:type])...))")
end

function show(io::IO, atom::Atom)
    for (k, v) in atom.props
        println(io, k, ": ", atom[k])
    end
end

display(atoms::Atoms) = show(atoms)
display(atom::Atom) = show(atom)
