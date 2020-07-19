# TropicalTensors

Solving tropical tensor network with Yao. It contains

* Spinglass solvers for three predefined lattices,
  square lattice, chimera lattice and second-nearest neighor coupled square lattice,
* CUDA programming and forward mode automatic differentiation,
* A reversible programming automatic differention.
* A visualization toolkit for spinglass solvers.

[![Build Status](https://travis-ci.com/TensorBFS/TropicalTensors.jl.svg?branch=master)](https://travis-ci.com/TensorBFS/TropicalTensors.jl)
[![Codecov](https://codecov.io/gh/TensorBFS/TropicalTensors.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TensorBFS/TropicalTensors.jl)

## To Start

Type `]` in Julia REPL to enter pkg mode and
```julia pkg
pkg> add https://github.com/JuliaReverse/TropicalYao.jl#master
pkg> dev git@github.com:TensorBFS/TropicalTensors.jl.git
```
The last line is required only when you use GPU for computing.

Then run the tests
```bash
$ cd ~/.julia/dev/TropicalTensors
$ julia test/runtests.jl
```
