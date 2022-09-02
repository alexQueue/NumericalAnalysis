push!(LOAD_PATH, "src")
using Documenter 

makedocs(
        sitename = "ProjectNumericalAnalysis",
        format = Documenter.HTML(prettyurls = false),
        pages = [
            "Home" => "index.md",
            "Introduction" => "introduction.md",
            "Mathematical Model" => "math_model.md",
            "Implementation" =>
                ["space_discretization.md", "time_discretization.md"], 
            "Results" => "results.md",
            "Conclusion" => "conclusion.md",
        ],
)

deploydocs(
    repo = "github.com/alexQueue/NumericalAnalysis",
    versions = nothing,
)
