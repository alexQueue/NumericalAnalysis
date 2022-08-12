using Documenter 

push!(LOAD_PATH, "src")

makedocs(
        sitename = "ProjectNumericalAnalysis",
        pages = [
            "Home" => "index.md",
            "Introduction" => "introduction.md",
            "Mathematical Model" => "math_model.md"
            "Implementation" =>
                ["space_discretization.md", "time_discretization.md"], 
            "Results" => "results.md"
            "Conclusion" => "conclusion.md"
        ],
)

deploydocs(
    repo = "github.com/alexQueue/NumericalAnalysis",
    versions = nothing,
)
