# RfbScFVM.jl: Performance Prediction of Flow Battery Cells

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

The RfbScFVM.jl package provides functionalities to predict the performance of single flow battery cells. The implemented model considers a simplified two-dimensional domain of the flow battery cell extending both in the through-plane direction of the electrode-membrane-assembly and the forced flow direction. The balance laws of species concentrations, electrostatic potentials in the electrode and electrolyte, pressure, and optionally temperature are discretized over the specified domain using the Voronoi-Finite-Volume scheme implemented in [VoronoiFVM](https://github.com/j-fu/VoronoiFVM.jl).

The main model features include:

- two-dimensional spatial discretization with the Voronoi-FV scheme
- prediction of polarization curves and power density vs current density
- simple handling of model inputs and outputs with JSON configuration files

**Important Notes**: 

- The current version of RfbScFVM.jl supports **Julia 1.9**.
- RfbScFVM.jl is a work in progress, which means that model interface changes are likely to occur.

### Installing Julia

Julia binaries can be installed directly from the official website [https://julialang.org](https://julialang.org/). Alternatively, Julia can be installed via the [Juliaup](https://github.com/JuliaLang/juliaup) cross-platform installer.

To check that the correct version of Julia is installed run `julia -e "println(VERSION)"`, or even simpler `julia --version`.

### Installing RfbScFVM.jl

Currently, the package can be installed locally by cloning this repository:
```bash
git clone https://github.com/Isomorph-Electrochemical-Cells/RfbScFVM.git
```

### Running RfbScFVM from the CLI

A simulation can be executed by running the script `run_simulation.jl` and passing a single argument with the path to the input configuration file. E.g., to run one simulation with the input configuration file `input/mv_temptma_polarization_soc50.json` we could execute
```bash
$ julia scripts/run_simulation.jl "input/mv_temptma_polarization_soc50.json"
```
Due to the just-in-time compilation of Julia, the runtime of the above command also includes possible (re-)compilation of code and package loading. To avoid this to occur between subsequent calls to Julia from the command line, we can make use of the [DaemonMode.jl](https://github.com/dmolina/DaemonMode.jl) package. This package uses a server/client model, where an instance of Julia runs in a background process (server), which allows subsequent calls to Julia from a client to re-use the last state of the Julia.

The `DaemonMode.jl` package can be installed directly from the command line with `julia -e "using Pkg; Pkg.add(\"DaemonMode\")"`.

As suggested on [DaemonMode.jl](https://github.com/dmolina/DaemonMode.jl), an alias for running Julia scripts as a client can be used as follows:
```bash
alias juliaclient='julia --startup-file=no -e "using DaemonMode; runargs()"'
```
#### Workflow


To start a Julia background process run:
```bash
$ julia -t auto --startup-file=no -e 'using DaemonMode; serve(async=true)' &
```
Here we set the number of threads automatically and allow for parallel execution of clients. See [DaemonMode.jl](https://github.com/dmolina/DaemonMode.jl) for a discussion of the individual options.

Subsequently, each simulation can be executed by running:
```bash
$ juliaclient scripts/run_simulation.jl "input/some_input_file.json"
```


## Acknowledgements

This work is part of the SONAR project supported by the European Union / European Commission funding program Horizon 2020.