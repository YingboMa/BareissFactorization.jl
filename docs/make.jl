using BareissFactorization
using Documenter

DocMeta.setdocmeta!(BareissFactorization, :DocTestSetup, :(using BareissFactorization); recursive=true)

makedocs(;
    modules=[BareissFactorization],
    authors="Yingbo Ma <mayingbo5@gmail.com> and contributors",
    repo="https://github.com/YingboMa/BareissFactorization.jl/blob/{commit}{path}#{line}",
    sitename="BareissFactorization.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://YingboMa.github.io/BareissFactorization.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YingboMa/BareissFactorization.jl",
    devbranch="master",
)
