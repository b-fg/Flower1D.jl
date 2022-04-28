push!(LOAD_PATH,"../src/")
using Flower1D, Documenter
makedocs(
         sitename = "Flower1D.jl",
         modules  = [Flower1D],
         pages=["Home" => "index.md"])
deploydocs(;
    repo="github.com/b-fg/Flower1D.jl",
)
