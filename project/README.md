# Note

## How to run the 6x6x6 cubic benchmark with multi-threading

```julia
pkg> add https://github.com/TensorBFS/TropicalTensors.jl.git#random-tensor
```

```bash
JULIA_NUM_THREADS=8 julia bigmem.jl 6
```
