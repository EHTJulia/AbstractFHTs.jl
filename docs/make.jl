using FastHartleyTransforms
using Documenter

DocMeta.setdocmeta!(FastHartleyTransforms, :DocTestSetup, :(using FastHartleyTransforms); recursive=true)

makedocs(;
    modules=[FastHartleyTransforms],
    authors="Kazunori Akiyama",
    sitename="FastHartleyTransforms.jl",
    format=Documenter.HTML(;
        canonical="https://EHTJulia.github.io/FastHartleyTransforms.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EHTJulia/FastHartleyTransforms.jl",
    devbranch="main",
)
