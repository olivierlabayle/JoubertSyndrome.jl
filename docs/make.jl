using JoubertSyndrome
using Documenter

DocMeta.setdocmeta!(JoubertSyndrome, :DocTestSetup, :(using JoubertSyndrome); recursive=true)

makedocs(;
    modules=[JoubertSyndrome],
    authors="Olivier Labayle",
    repo="https://github.com/olivierlabayle/JoubertSyndrome.jl/blob/{commit}{path}#{line}",
    sitename="JoubertSyndrome.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olivierlabayle.github.io/JoubertSyndrome.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olivierlabayle/JoubertSyndrome.jl",
)
