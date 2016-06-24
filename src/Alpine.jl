module Alpine

using Distances, VPTrees

include("types.jl")

export
    # types.jl
    Indices,
    Clusters,

    # distances.jl
    sd,
    msd,
    rmsd,
    mag,
    wrapcoords!,
    unwrapcoords!,
    distancesq,
    distances!,

    # atoms.jl
    group,
    coalesce,

    # files.jl
    Frame,
    AtomData,
    read_atomfile,
    read_atomfile_structure,
    get_atom_indices,
    types,
    coords,
    charges,
    energies,
    num_atoms,
    timestep,
    box_dims,

    # analysis.jl
    density,
    get_atoms,

    # rdf.jl
    rdf_slow,
    rdf,

    # neighborhood_similarity.jl
    allocate_clusters

include("distances.jl")
include("atoms.jl")
include("files.jl")
include("analysis.jl")
include("rdf.jl")
include("neighborhood_similarity.jl")

end
