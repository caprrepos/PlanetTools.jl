# PlanetTools.jl

```@docs
PlanetTools
```

## Installation

PlanetTools is registered in the Julia General registry:

```julia
using Pkg
Pkg.add("PlanetTools")
```

The development version is installed with
`Pkg.add(url="https://github.com/caprrepos/PlanetTools.jl")`.

PlanetTools requires Julia 1.10 or later. The plotting methods `compose_palette` and
`finite_contour!` additionally require Makie or one of its backends, e.g.
`Pkg.add("CairoMakie")`.
