using TropicalTensors

sg = rand_spinglass(Int32, SquareLattice(7, 7); jt=Randpm(), ht=Zero(), seed=2)
res = opt_config(sg)
res |> display
