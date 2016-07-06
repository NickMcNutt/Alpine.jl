module Alpine

using VPTrees, CompressedIndices

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
    Frame,
    Atoms,
    Atom,
    add_property!,

    # atoms_analysis.jl
    group,
    box_dims,
    unwrap!,

    # files.jl
    AtomFileType,
    AtomFile,
    FrameHeader,
    read_atoms,
    read_frame_headers,
    write_atoms,
    LAMMPSTRJ, TRJ,
    LAMMPSDAT, DAT,
    XYZ,

    # analysis.jl
    density,
    get_atoms,

    # rdf.jl
    rdf_slow,
    rdf,
    rdf_axis

include("distances.jl")
include("atoms.jl")
include("atoms_analysis.jl")

include("files.jl")
include("files/lammpstrj.jl")
include("files/lammpsdat.jl")
include("files/xyz.jl")

include("analysis.jl")
include("rdf.jl")

end
