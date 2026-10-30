# PlanetTools.jl

```@docs
PlanetTools
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/caprrepos/PlanetTools.jl")
# or, when it has become part of the General Registry:
Pkg.add("PlanetTools")
```

PlanetTools requires Julia 1.10 or later. The plotting methods `compose_palette` and
`finite_contour!` additionally require Makie or one of its backends, e.g.
`Pkg.add("CairoMakie")`.
