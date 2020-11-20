# Note

## How to run the 6x6x6 cubic benchmark with multi-threading

```julia
pkg> add YaoArrayRegister#multi-threading-ngates

pkg> dev https://github.com/TensorBFS/TropicalTensors.jl.git
```

```bash
$ cd .julia/dev/TropicalTensors
$ git checkout random-tensor
$ JULIA_NUM_THREADS=8 julia project/bigmem.jl 6
```
