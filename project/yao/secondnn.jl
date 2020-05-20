using TropicalTensors
using Viznet, Compose

sg = rand_spinglass(Int32, rand_maskedsquare(12, 8, 0.5); jt=Randpm(), ht=Zero(), seed=2)
res = opt_config_largemem(sg)

lt = sg.lattice
g =showlattice(lt)


lt = MaskedSquareLattice([0 1; 1 0; 0 1])
sg = rand_spinglass(Int64, lt; jt=[-1,1], ht=Zero(), seed=7)
solve(sg).n == 2
#=
dash = compose(context(), bondstyle(:dashed), stroke("black"))
g2 = canvas() do
    for i=2.5:2:size(lt, 1)
        dash >> lt[(i, 0); (i, size(lt, 2)+1)]
    end
    for j=2.5:2:size(lt, 2)
        dash >> lt[(0, j); (size(lt, 1)+1, j)]
    end
end
compose(context(), g, g2)
=#
