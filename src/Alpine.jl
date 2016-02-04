module Alpine

using Distances, VantagePointTrees

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

    # geometry.jl
    rotation,
    randomrotation,

    # atoms.jl
    group,
    coalesce,

    # trj.jl
    readtrj,

    # analysis.jl
    density,
    rdf,

    # neighborhood_similarity.jl
    allocate_clusters

include("distances.jl")
include("geometry.jl")
include("atoms.jl")
include("trj.jl")
include("analysis.jl")
include("neighborhood_similarity.jl")

end
