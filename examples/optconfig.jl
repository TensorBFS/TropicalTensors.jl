using TropicalTensors

sg = rand_spinglass(Int32, SquareLattice(7, 7); jt=Randpm(), ht=Zero(), seed=2)
eng, g = opt_config(sg)
using Compose
set_default_graphic_size(12cm, 12cm)
vizgrad_J(sg, g)
