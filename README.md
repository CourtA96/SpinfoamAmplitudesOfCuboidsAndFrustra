# CuboidsAndFrustra

Julia code for numerically evaluating the amplitudes of spinfoams composed of cuboids and frustra. Please refer to our paper: https://arxiv.org/pdf/2201.09902.pdf for information about how the intertwiners are labeled and so on.

## Basic Usage

The code is broken into separate directories for the frustra and cuboid cases. Each directory contains four files:

 - createIntertwiners.jl
    - Module that creates the intertwiners.
 - ContractionMod.jl
    - Module that contracts the intertwiners.
 - createTest.jl
    - Gives an example of how to create intertwiners.
 - conTest.jl
    - Gives and example of how to contract intertwiners.

### `createIntertwiners.jl`

To illustrate how to use this module, take the case where all spins are one half. To use this module in the frustra case:

```
julia> include("frustra\\createIntertwiners.jl")
julia> createIntertwiners.createInts([1,1,1])
```

To use this module in the cuboid case:

```
julia> include("frustra\\createIntertwiners.jl")
julia> createIntertwiners.createInts([1,1,1,1,1,1])
```

The `createInts(Jvector)` computes the intertwiners in a spinfoam defined by the spins in the vector `Jvector` which is a 1D array of `Int64`. So if the spinfoam is defined by spins *j*=[0.5,0.5,0.5], then `Jvector = floor(2`*j*`)`. So in the above examples, the spins are all equal to one half.

The `createInts` function saves the intertwiners to `.txt` files named `IntNj[x,y,z].txt` where `N` is the intertwiner's label and `[x,y,z]` are the spins that define the intertwiner. For example, intertwiner 1 in the spin one half case would be saved as `Int1j[1,1,1].txt`. To save time, `createInts` automatically checks the directory to see if the intertwiner is already in the directory, if it is, the intertwiner is not computed again.

== Note: Since both the frustra and cuboid cases obey the same naming convention, the intertwiners for each case must be computed and stored in separate directories to avoid confusion. ==

### `ContractionMod.jl`

To use this module in the frustra case (again, using the spin one half case as an example):

```
julia> include("ContractionMod.jl")
julia> Jvector = [1,1,1]
julia> int1,int2,int3,int4,int5,int6,int7,int8 = ContractionMod.loadInts(Jvector)
julia> contract = ContractionMod.Contraction(int1,int2,int3,int4,int5,int6,int7,int8)
```

In the cuboid case, everything is the same as above except `Jvector = [1,1,1,1,1,1]`.

`Jvector` is the vector of spins that define the spinfoam, as in `createIntertwiners.jl`. `loadInts(jVector)` reads in the intertwiners from the `.txt` files, and `Contraction(int1,int2,int3,int4,int5,int6,int7,int8)` performs the contraction.

# Authors of the Code:

 - Courtney Allen (University of Guelph, Canada)
    - Contact at: callen15@uoguelph.ca
 - Sebastian Steinhaus (Friedrich-Schiller-Universitat Jena, Germany)
