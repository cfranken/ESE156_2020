using IJulia
using Pkg.Artifacts
""" Shorthand for @artifact_str """
artifact(name) = joinpath(@artifact_str(name), name)
jupyterlab(dir=joinpath(dirname(@__FILE__), ".."))
