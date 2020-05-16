# TropicalTensors

Tensors of Tropical numbers.

[![Build Status](https://travis-ci.com/TensorBFS/TropicalTensors.jl.svg?branch=master)](https://travis-ci.com/TensorBFS/TropicalTensors.jl)
[![Codecov](https://codecov.io/gh/TensorBFS/TropicalTensors.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TensorBFS/TropicalTensors.jl)

## To Start

Type `]` in Julia REPL to enter pkg mode and
```julia pkg
pkg> add https://github.com/GiggleLiu/NiLangCore.jl#log-numbers add https://github.com/GiggleLiu/NiLang.jl#matchcore https://github.com/GiggleLiu/Viznet.jl#master https://github.com/JuliaReverse/TropicalYao.jl#master https://github.com/JuliaReverse/NiLogLikeNumbers.jl#master https://github.com/GiggleLiu/TropicalNumbers.jl#master
pkg> dev git@github.com:TensorBFS/TropicalTensors.jl.git
```
The last line is required only when you use GPU for computing.

Then run the tests
```bash
$ cd ~/.julia/dev/TropicalTensors
$ julia test/runtests.jl
```
