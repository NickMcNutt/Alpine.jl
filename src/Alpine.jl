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

    # trj.jl
    readtrj,

    # analysis.jl
    density,
    get_atoms,
    rdf,

    # neighborhood_similarity.jl
    allocate_clusters

include("distances.jl")
include("atoms.jl")
include("trj.jl")
include("analysis.jl")
include("neighborhood_similarity.jl")

end
