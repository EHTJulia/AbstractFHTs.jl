using AbstractFastHartleyTransforms
using Documenter

DocMeta.setdocmeta!(AbstractFastHartleyTransforms, :DocTestSetup, :(using AbstractFastHartleyTransforms); recursive=true)

makedocs(;
    modules=[AbstractFastHartleyTransforms],
    authors="Kazunori Akiyama",
    sitename="AbstractFastHartleyTransforms.jl",
    format=Documenter.HTML(;
        canonical="https://EHTJulia.github.io/AbstractFastHartleyTransforms.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EHTJulia/AbstractFastHartleyTransforms.jl",
    devbranch="main",
)
