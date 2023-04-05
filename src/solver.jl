function solve_steady_state_problem(system, inival; verbose=false)
    control = VoronoiFVM.SolverControl()
    control.verbose = verbose

    # Stationary solution of the problem
    solution = VoronoiFVM.solve(system, control=control, inival=inival,
                                time=0.0, tstep=Inf, abstol=1e-8, reltol=1e-8, maxiter=10)
    return solution
end

function solve_transient_problem(system, inival, times; verbose=false)
    # overall time span of the problem
    ΔT = times[end] - times[1]
    # minimum time span between two time stamps
    Δt = minimum(times[2:end] - times[1:(end-1)])

    control = VoronoiFVM.SolverControl()
    control.verbose = verbose
    control.Δt = 0.01*Δt
    control.Δt_min = 0.01*Δt
    control.Δt_max = 0.1*ΔT
    solution = VoronoiFVM.solve(system,
        inival=inival, abstol=1e-8, reltol=1e-8, times=tspan, control=control)
end
