using Documenter, Astrometry

makedocs(;
    modules=[Astrometry],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Example Usage" => "exampleusage.md",
        "APIs" => "internalapis.md",
    ],
    repo="https://github.com/pbarrett/Astrometry.jl/blob/{commit}{path}#L{line}",
    sitename="Astrometry.jl",
    authors="Paul Barrett and contributors",
)

deploydocs(;
    repo="github.com/pbarrett/Astrometry.jl",
    push_preview=true
)
