
using Pkg
Pkg.activate(".")
using RfbScFVM

function run_simulation()
    @assert length(ARGS)==1 && typeof(ARGS[1])===String
    RfbScFVM.run(ARGS[1])
end

run_simulation()
