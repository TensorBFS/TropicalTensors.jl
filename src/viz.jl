export viz

function Viznet.viz(sg::Spinglass; r=0.3*unit(sg.lattice))
    lt = sg.lattice
    nb0 = nodestyle(:default, fill("white"), stroke("black"); r=r)
    nb1 = nodestyle(:default, fill("skyblue"), stroke("black"); r=r)
    nb2 = nodestyle(:default, fill("orange"), stroke("black"); r=r)
    eb1 = bondstyle(:default, linewidth(4mm*unit(lt)), stroke("skyblue"))
    eb2 = bondstyle(:default, linewidth(4mm*unit(lt)), stroke("orange"))
    tb = textstyle(:default, fontsize(150pt*unit(lt)))
    cdots = canvas() do
        for (i, v) in enumerate(sgvertices(lt))
            (sg.hs[i]>0 ? nb1 : (sg.hs[i]<0 ? nb2 : nb0)) >> lt[v]
            tb >> (lt[v], "$v")
        end
        for ((i,j),v) in zip(sgbonds(lt), sg.Js)
            if v > 0
                eb1 >> lt[i;j]
            elseif v < 0
                eb2 >> lt[i;j]
            end
        end
    end
    compose(context(), cdots)
end
function Base.show(io::IO, mime::MIME"text/html", sg::Spinglass)
    Base.show(io, mime, viz(sg))
end

Viznet.viz(res::SpinglassOptConfig; r=0.3*unit(res.sg.lattice)) = vizgrad_J(res.sg, res.grad_Js, res.grad_hs; r=r)
function Base.show(io::IO, mime::MIME"text/html", sg::SpinglassOptConfig)
    Base.show(io, mime, viz(sg))
end

export vizgrad_J
function vizgrad_J(sg::Spinglass, grad_Js::AbstractVector, grad_hs::AbstractVector; r=0.3*unit(sg.lattice))
    lt = sg.lattice
    grid = assign_Js_hs(lt, grad_Js, grad_hs)
    nb1 = nodestyle(:default, fill("white"), stroke("black"); r=r)
    nb2 = nodestyle(:default, fill("black"), stroke("black"); r=r)
    eb1 = bondstyle(:default, linewidth(4mm*unit(lt)), stroke("skyblue"))
    eb2 = bondstyle(:default, linewidth(4mm*unit(lt)), stroke("orange"))
    cdots = canvas() do
        for i in sgvertices(lt)
            if grid[i] > 0
                nb1 >> lt[i]
            elseif grid[i] < 0
                nb2 >> lt[i]
            else
                error("index $i not set!")
            end
        end
        for ((i,j),v) in zip(sgbonds(lt), sg.Js)
            if v > 0
                eb1 >> lt[i;j]
            else
                eb2 >> lt[i;j]
            end
        end
    end
    compose(context(), cdots)
end

# visualize lattice in spin-glass definition order
export viz_lattice

function viz_lattice(lt; line_style=bondstyle(:default, stroke("black"), linewidth(4mm*unit(lt))),
        node_style=nodestyle(:default, stroke("black"), fill("white"), linewidth(2mm*unit(lt)); r=0.3*unit(lt)),
        text_style=textstyle(:default, fontsize(150pt*unit(lt))))
    Viznet.empty_cache!()
    for node in sgvertices(lt)
        node_style >> lt[node]
        text_style >> (lt[node], "$node")
    end
    for bond in sgbonds(lt)
        line_style >> lt[bond[1]; bond[2]]
    end
    Viznet.flush!()
end

