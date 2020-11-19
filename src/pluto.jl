using .Pluto

@info "You can use the notebooks now by typing, e.g. `TropicalTensors.notebook(\"spinglass\")`"

export notebook

project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

"""
    notebook(which; dev=false, kwargs...)

Open a notebook, the first argument can be
* "spinglass": solving spinglass model with Yao.
* "randomgraph": solving randomgraph model with tensor network contraction.
"""
function notebook(which; dev=false, kwargs...)
    src = project_relative_path("notebooks", "$which.jl")
    if dev
        dest = src
    else
        dest = tempname()
        cp(src, dest)
        chmod(dest, 0o664)
    end
    Pluto.run(; notebook=dest, project=project_relative_path(), kwargs...)
end
