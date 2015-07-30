module Alpine

include("types.jl")

export
    # geometry.jl
    Indices,
    Coords,
    AbstractIndices,
    AbstractCoords,

    normN,
    norm3,
    normFrobenius,
    rmsd,
    kNN!,
    kNN,
    centerOfMass!,
    centerOfMass,
    rodriguesRotationMatrix,
    rotateCoords!,
    rotateCoordsVector!,
    transposeCoords!,
    unwrapCoords!,
    wrapCoord!,
    wrapCoords!,
    alignWithAxes!,
    rdf,
    @distanceSqPeriodic,
    distanceSqPeriodic,
    distancePeriodic,
    distancesSqPeriodic!,
    distancesPeriodic!,

    # atoms.jl
    getIndicesByAtomType,
    getIndices,
    createGroupsFromAtomTypes,
    getGroupCoords,
    filterGroup,
    groupBy,
    coalesce,

    # files.jl
    AtomFile,
    AtomStream,
    
    start,
    next,
    done,
    length,
    endof,
    first,
    last,
    getindex,
    close,
    open,
    readFrame,
    readPos,
    newFilename,

    # neighborhood_similarity.jl
    Group,
    Neighbors,
    NeighborClusters,

    push!,
    nonzeroNeighbors,
    allocate_clusters,
    kabsch,
    kabsch2,
    inner_product_spectrum,
    gram,
    pd,
    optimalperm,
    gram2,
    gram3!,
    gramMatrix,
    gramMatrix!,
    d_ref

include("geometry.jl")
include("atoms.jl")
include("files.jl")
include("neighborhood_similarity.jl")

end
