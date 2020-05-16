using TropicalTensors

sg = rand_spinglass(Int32, SquareLattice(18, 18); jt=Randpm(), ht=Zero(), seed=2)
res = opt_config(sg)
res |> display
