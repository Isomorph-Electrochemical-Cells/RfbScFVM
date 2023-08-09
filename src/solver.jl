function solve_steady_state_problem(system, inival, applied_voltage_values; verbose=false)
    control = VoronoiFVM.SolverControl(verbose=true,
                                       handle_exceptions=true)

    voltage_start = applied_voltage_values[1]
    voltage_stop = applied_voltage_values[end]

    norm_embed = (applied_voltage_values .- voltage_start)/(voltage_stop - voltage_start)

    @assert maximum(norm_embed)<=one(eltype(norm_embed)) &&
            minimum(norm_embed)>=zero(eltype(norm_embed))

    norm_step = maximum(norm_embed[2:end]-norm_embed[1:(end-1)])

    maxiters = system.physics.data.discr.maxiters
    abstol = system.physics.data.discr.abstol
    reltol = system.physics.data.discr.reltol
    # Stationary solution of the problem
    solution = VoronoiFVM.solve(system, control=control, inival=inival, embed=norm_embed,
                                Δp=norm_step, Δp_max=norm_step, Δp_min=norm_step/10,
                                maxiters=maxiters, abstol=abstol, reltol=reltol, Δu_opt=1,
                                pre=(sol, embed)->pre_step(system, sol, embed,
                                                           voltage_start, voltage_stop),
                                post=(sol, oldsol, embed, Δembed)->post_step(system, sol))
    return solution
end

function pre_step(system, sol, norm_embed, voltage_start, voltage_stop)
    cell_voltage = voltage_start + (voltage_stop-voltage_start) * norm_embed

    #@info "applied cell voltage: " cell_voltage

    system.physics.data.boundary.ϕₛ_neg[] = -cell_voltage/2
    system.physics.data.boundary.ϕₛ_pos[] = cell_voltage/2

end

function post_step(system, solution)
    data = system.physics.data
    current_density_values = integrate_reaction_terms(system, solution, data)
    current_density_values /= (data.geom.ly_cell * data.scaling_params.ϵL0)

    #check that charge is conserved by asserting that the current density
    # is the same in both half-cells
    @assert isapprox(current_density_values[1], -current_density_values[2], atol=1e-6)

    push!(data.results.current_density, current_density_values[2])
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
    solution = VoronoiFVM.solve(system, inival=inival, abstol=1e-8, reltol=1e-8,
                                times=tspan, control=control)
    return solution
end
