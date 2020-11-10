export viz_sg, viz_lattice

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

function viz_sg(sg::Spinglass; r=0.3*unit(sg.lattice))
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

function Base.show(io::IO, mime::MIME"text/html", lt::Viznet.AbstractLattice)
    Base.show(io, mime, viz_lattice(lt))
end

function Base.show(io::IO, mime::MIME"text/html", sg::Spinglass)
    Base.show(io, mime, viz_sg(sg))
end

export viz_optconfig
function Base.show(io::IO, mime::MIME"text/html", sg::SpinglassOptConfig)
    Base.show(io, mime, viz_optconfig(sg))
end

viz_optconfig(res::SpinglassOptConfig; r=0.3*unit(res.sg.lattice)) = vizgrad_J(res.sg, res.grad_Js, res.grad_hs; r=r)

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

export viz_tnet

function viz_tnet(tnet::TensorNetwork; r=0.25/sqrt(length(tnet.tensors)+1), show_edgeindex=false, node_fontsize=100pt/sqrt(length(tnet.tensors)+1), edge_fontsize=200pt/sqrt(length(tnet.tensors)+1))
    nt = length(tnet.tensors)
    nb = nodestyle(:default, fill("white"), stroke("black"), linewidth(2mm/sqrt(length(tnet.tensors)+1)); r=r)
    eb = bondstyle(:default, linewidth(4mm/sqrt(length(tnet.tensors)+1)), stroke("skyblue"))
    tb1 = textstyle(:default, fontsize(node_fontsize))
    tb2 = textstyle(:default, fontsize(edge_fontsize))
    compose(context(r, r, 1-2r, 1-2r), canvas() do
        for (t, meta) in zip(tnet.tensors, tnet.metas)
            nb >> meta.loc
            if !isempty(meta.name)
                tb1 >> (meta.loc, meta.name)
            end
        end
        for i=1:nt
            for j=i+1:nt
                li = tnet.tensors[i].labels
                lj = tnet.tensors[j].labels
                loci, locj = tnet.metas[i].loc, tnet.metas[j].loc
                common_labels = li ∩ lj
                if !isempty(common_labels)
                    eb >> (loci, locj)
                    show_edgeindex && tb2 >> ((loci .+ locj) ./ 2, join(common_labels, ", "))
                end
            end
        end
    end)
end

function Base.show(io::IO, mime::MIME"text/html", tnet::TensorNetwork)
    show(io, mime, viz_tnet(tnet))
end
