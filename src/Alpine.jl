module AtomTools
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
    weylMatrix!,
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
    newFilename

include("geometry.jl")
include("atoms.jl")
include("files.jl")

end
