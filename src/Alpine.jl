module Alpine
using Devectorize

export
    # geometry.jl
    centerOfMass,
    rodriguesRotationMatrix,
    rotateCoords!,
    rotateCoordsVector!,
    transposeCoords!,
    unwrapCoords!,
    wrapCoord!,
    wrapCoords!,
    distanceSqPeriodic,
    alignWithAxes!,
    rdf,

    # atoms.jl
    whichEnergyGroup,
    atomIndicesToCoords,
    atomIndicesToCoords!,
    atomIndexToCoords,
    atomIndexToCoords!,

    # files.jl
    XYZFile,
    XYZStream,
    
    convert,
    start,
    next,
    done,
    length,
    getindex,
    open,
    readFrame,
    readPos,
    newFilename,

    # neighborhood_similarity.jl
    gramMatrix!

include("geometry.jl")
include("atoms.jl")
include("files.jl")
include("neighborhood_similarity.jl")

end
