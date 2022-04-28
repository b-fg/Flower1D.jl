# Flower1D.jl

Flux reconstruction fluid flow solver in Julia.

<!-- ![](assets/viscous_burgers.png?raw=true "Viscous burgers equation") -->
<p align="center">
	<img align="center" width="600" src="assets/viscous_burgers.png?raw=true">
</p>

## Overview

**Flower1D.jl** is a pure Julia implementation of a [flux reconstruction](https://www.researchgate.net/publication/309715679) (high-order) solver of 1D PDEs: Linear advection, Burgers, viscous Burgers, and Euler equations.
Currently, it is used for educational and academic purposes, but future developments will be focused in 2D and 3D versions for unstructured grids.

## Numerical methods
The mesh can be arbitrarily discretized into `N` elements, and each element is represented with Lagrange polynomials of degree `P`, hence using `P + 1` solution points within each element. Furthermore, a continuos flux is obtained with [VCJH](https://doi.org/10.1007/s10915-010-9420-z) correction schemes, governed by the `C` parameter, which can recover different high-order formulations.
The solution points within each element can follow an `"equidistant"` distribution or a `"legendre"` distribution.
The time-integration methods available are Euler (1st order) and Runge-Kutta (2nd and 3rd order).

## Installation, usage, and development

To install the package (plus dependencies) use
```
julia
(v1.x) pkg> add "git@github.com:b-fg/Flower1D.jl.git"
```

You can now download the examples and run them by previously loading the package
```
using Flower1D
```

If you wish to download the source code and play with it: clone the repo, install the package, and mark it for local development
```
git clone git@github.com:b-fg/Flower1D.jl
cd Flower1D.jl
julia
(v1.x) pkg> add .
(v1.x) pkg> dev --local Flower1D
```
Now you can run any example using
```
cd Flower1D.jl
julia --project
julia > include("examples/burgers_advection_diffusion/main.jl")

```
Modifications in the source code will be re-compiled with `Revise.jl`.

Note that examples for all the available PDEs and flux reconstruction functionalities can be found in `examples/`.
Documentation will soon be available in `docs/`.

## Contributions and future work
Contributions are most welcomed!
Please submit a pull request and we will work together to include your changes in the package 😄.

Future features development will be related to:
- Extension to structured 2D and 3D geometries.
- Extension to unstructured grids.
- Artificial viscosity for stabilisation.
- Compressible Navier-Stokes equations.
- Parallel computing using [`Distributed.jl`](https://docs.julialang.org/en/v1/manual/distributed-computing/).
- Adaptive h/p refinement.
- VTK input/output.

## Cite the repo!

If you find this repository useful, please cite it:

```
@misc{Flower1D.jl,
  author = {B. Font},
  title = {Flower1D.jl: Flux reconstruction fluid flow solver for 1D {PDE}s},
  year = {2022},
  publisher = {GitHub},
  journal = {GitHub repository},
  url = {https://github.com/b-fg/Flower1D.jl},
}
```