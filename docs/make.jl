using AbstractFHTs
using Documenter

DocMeta.setdocmeta!(FastHartleyTransforms, :DocTestSetup, :(using FastHartleyTransforms); recursive=true)

makedocs(;
    modules=[AbstractFHTs],
    authors="Kazunori Akiyama",
    sitename="AbstractFHTs.jl",
    format=Documenter.HTML(;
        canonical="https://EHTJulia.github.io/AbstractFHTs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EHTJulia/AbstractFHTs.jl",
    devbranch="main",
)
