using DrWatson
@quickactivate "RailwayOptSitePos"

# Here you may include files from the scripts/source directory
# Running the model with the Simulated Data
include(scriptsdir("model_solvit.jl"))

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Have fun with your new project!

You can help us improve DrWatson by opening
issues on GitHub, submitting feature requests,
or even opening your own Pull Requests!
"""
)
