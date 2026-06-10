using Documenter
using PVlib

DocMeta.setdocmeta!(PVlib, :DocTestSetup, :(using PVlib); recursive = true)

makedocs(
    modules = [PVlib],
    sitename = "PVlib.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = ["Home" => "index.md", "API" => "api.md"],
)

deploydocs(repo = "github.com/jtgrasb/PVlib.jl", devbranch = "main")
