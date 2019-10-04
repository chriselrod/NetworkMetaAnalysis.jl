using Documenter, NetworkMetaAnalysis

makedocs(;
    modules=[NetworkMetaAnalysis],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/chriselrod/NetworkMetaAnalysis.jl/blob/{commit}{path}#L{line}",
    sitename="NetworkMetaAnalysis.jl",
    authors="Chris Elrod",
    assets=String[],
)

deploydocs(;
    repo="github.com/chriselrod/NetworkMetaAnalysis.jl",
)
