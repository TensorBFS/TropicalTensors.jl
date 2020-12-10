# Note

```julia
pkg> dev https://github.com/TensorBFS/TropicalTensors.jl.git
```

```bash
$ cd .julia/dev/TropicalTensors
```

1. Computing the maximum energy of a 5x5x5 cubic lattice spin glass on CPU by quantum simulation.
```
$ JULIA_NUM_THREADS=8 julia project/cubic_spinglass.jl 5
```
The parameter `5` is the size of the cubic lattice.

2. Sample the degeneracy of spinglass on a cylinder on GPU by quantum simulation.
```
$ julia project/counting_cylinder.jl 0 18
```
The first parameter `0` is the GPU device ID, `18` is the lattice size in both `x` and `y` (periodic) direcitons.

3. Sample the degeneracy of 3-state Potts model on GPU by tensor network contraction.
```
$ julia project/counting_potts.jl 0 15
```
The parameters are GPU device ID, and square lattice size.
The default number of samples is 100.

4. Sample the residual entropy for the Ising model and the 2-SAT model defined on a random 3-regular graph on GPU by tensor network contraction.
```
$ julia project/counting_ising_and_2sat.jl 0 ising 140

$ julia project/counting_ising_and_2sat.jl 0 2sat 140
```
The parameters are GPU device ID, model name and graph size.
The graph size should be one of 60, 80, 120, 140, 160, 180 and 200.
The default number of samples is 100.