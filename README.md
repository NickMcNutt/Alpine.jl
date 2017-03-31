# Alpine.jl

[![Build Status](https://travis-ci.org/NickMcNutt/Alpine.jl.svg?branch=master)](https://travis-ci.org/NickMcNutt/Alpine.jl)
[![Coverage Status](https://coveralls.io/repos/NickMcNutt/Alpine.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/NickMcNutt/Alpine.jl?branch=master)
[![codecov.io](http://codecov.io/github/NickMcNutt/Alpine.jl/coverage.svg?branch=master)](http://codecov.io/github/NickMcNutt/Alpine.jl?branch=master)

## About

Alpine is a Julia library for analyzing and processing the results of molecular dynamics simulations.

### Overview

* Load and save a variety of molecular dynamics data/trajectory files (XYZ, LAMMPSTRJ, LAMMPSDAT, etc.)

* Work with extremely large (1 TB+) data files by indexing frames and loading only the relevant data

* Key/value property system supports arbitrary properties defined on atoms and simulation boxes

* Quickly perform *k* nearest neighbor and range searches within large simulation boxes using Vantage Point trees (which can handle periodic box boundaries)

* Very fast computation of radial distribution functions (RDFs) with support for arbitrary binning functions

* Ability to generate 3D density distributions from atom position data

* Licensed under the [MIT License](https://opensource.org/licenses/MIT)

### How to install

In Julia:
```julia
Pkg.clone("https://github.com/NickMcNutt/Alpine.jl")
```

## Usage

### Reading and writing molecular dynamics data/trajectory files

```julia
AtomFile(filename::AbstractString)

AtomFile(filetype::AtomFileType, filename::AbstractString)
```

Create a reference to a molecular dynamics data file.  The file type is determined from the file extension.  The second form of this function may be used to directly specify the file type, which can be one of:

* XYZ
 * An XYZ file. The first line of this file must be an integer specifying the number of atoms.  The second line of the file can be anything. Multiple frames can be contained in the XYZ file.  Each frame must immediately follow the previous one.
* LAMMPSTRJ
 * A LAMMPS trajectory file
* LAMMPSDAT
 * A LAMMPS data file

 Note that Alpine **DOES NOT CHECK** that the specified file format is correct, or even that the data file is error free. This is because a data integrity check would drastically slow down the reading of large files. To this end, ensure that your data is correct, otherwise the results will be unpredictable.


```julia
read_frame_headers(atom_file::AtomFile)
```

Collect header information from all of the simulation frames in a molecular dynamics trajectory file. This function indexes the location of each frame within the data file so that atom data can later be quickly read without searching the entire file. This function also reads any properties defined on the simulation box in each frame (e.g., current timestep, box dimensions, etc.).

```julia
read_atoms(atom_file::AtomFile; [frames::AbstractVector{Int}])

read_atoms(filename::AbstractString, [end_frames::Int])
```

Read atom data from an AtomFile.

The first form of this function takes an AtomFile and an optional keyword argument `frames` that specifies the indices of the frames that are to be read from a trajectory file.

The second form of this function is a shortcut that creates an AtomFile from a given `filename` and then reads the atom data from that file. The optional argument `end_frames` specifies that only the last `end_frames` frames in the trajectory file should be read.

```julia
write_atoms(atom_file::AtomFile, frame::Frame)

write_atoms(atom_file::AtomFile, frames::Vector{Frame})

write_atoms(filename::AbstractString, frame::Frame)

write_atoms(filename::AbstractString, frames::Vector{Frame})
```

Write one or more `Frame`s to an atom file.

### Analysis

Alpine includes a collection of tools for analyzing the results of MD simulations.

For computing [radial distribution functions](https://en.wikipedia.org/wiki/Radial_distribution_function) (RDFs), the following functions are provided:

```julia
rdf{T}(atoms::Atoms, xbw::T, ybw::T, zbw::T, radius::T, Δr::T)
```

Compute the RDF for a collection of `atoms` in a periodic rectangular prism box with dimensions (`xbw`, `ybw`, `zbw`) out to radius `radius` with bin size `Δr`.

```julia
rdf{T}(f::Function, nf::Int, atoms::Atoms, xbw::T, ybw::T, zbw::T, radius::T, Δr::T)
```

Compute the multiple RDFs simultaneous for a collection of `atoms` where each RDF consists of atoms that satisfy a particular condition.
