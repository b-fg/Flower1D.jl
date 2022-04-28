using Glob, DelimitedFiles, LaTeXStrings

"""
    quadPointsAndWeights(elementType = 1, p = 2, quadType = "legendre"; T = Float64)

Reads the integration point location and the Gauss-quadrature weight according to the (integer) polynomial order `p`,
and the (string) `quadratureType`. The only supported `elementType` is 1.
The `quadratureType` takes can take the value:
- "legendre"
- "equidistant"
- "lobatto"
Note that the Gauss-Lobatto flux reconstruction is still not implemented.
"""
function quadPointsAndWeights(elementType = 1, p = 2, quadType = "legendre"; T = Float64)
    projectPath(parts...) = normpath(joinpath(@__DIR__, parts...))
    src = projectPath("quadratures/")
    for (root, dirs, _) in walkdir(src)
        for dir in dirs
            if occursin(string(elementType), dir)
                file = glob("*" * quadType * "*" * string(p) * "*", joinpath(root, dir))[1]
                data = readdlm(file, '\t', T, skipstart = 1)
                return data[:, 1] , data[:, 2]
            end
        end
    end
end

"""
    out(str, value)

Pretty print given a `str` variable name and the `value` fo the variable.
"""
function out(str, value)
    println(str)
    display(value)
end

"""
    function pifrac(x)

Plotting with a custom scale [0 π/2 π 3π/2 2π]

https://stackoverflow.com/questions/58962826/plotting-with-a-custom-scale-0-%CF%80-2-%CF%80-3%CF%80-2-2%CF%80-in-julia
"""
function pifrac(x)
    isintegralmultipleof(f, x) = f != 0 && x / f ≈ round(x / f)
    if isintegralmultipleof(pi/2, x)
        n = Int(round(2 * x / π))
        str = iseven(n) ? "\$ $(n ÷ 2)\\pi \$" : "\$ $(n)\\pi/2 \$"
        str = replace(str, " 0\\pi" => " 0")
        str = replace(str, " 1\\pi" => " \\pi")
        return str
    end
    return latexstring(x)
end