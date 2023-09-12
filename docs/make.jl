using Documenter
using RfbScFVM

makedocs(
    sitename = "RfbScFVM",
    format = Documenter.HTML(),
    modules = [RfbScFVM]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
