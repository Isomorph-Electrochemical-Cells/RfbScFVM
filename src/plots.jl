function plot_polarization_curve(dict_result)

    current_density_unit = dict_result["current_density"]["unit"]
    current_density_values = dict_result["current_density"]["value"]
    voltage_unit = dict_result["cell_voltage"]["unit"]
    voltage_values = dict_result["cell_voltage"]["value"]
    power_density_unit = dict_result["power_density"]["unit"]
    power_density_values = dict_result["power_density"]["value"]

    CairoMakie.activate!()
    fontsize_theme = Theme(fontsize = 22)
    set_theme!(fontsize_theme)

    fig = Makie.Figure()
    fig_axis_voltage = Makie.Axis(fig[1,1];
                                  title="Polarization Curve",
                                  xlabel="Current Density [" * current_density_unit * "]",
                                  ylabel="Voltage [" * voltage_unit * "]",
                                  yticklabelcolor = :blue)
    fig_axis_voltage.xreversed = true
    fig_axis_power_density = Makie.Axis(fig[1,1];
                                        yticklabelcolor = :red,
                                        ylabel="Power Density [" * power_density_unit * "]",
                                        yaxisposition = :right)
    fig_axis_power_density.xreversed = true
    hidespines!(fig_axis_power_density)
    hidexdecorations!(fig_axis_power_density)

    lines_voltage = lines!(fig_axis_voltage,
                          current_density_values, voltage_values,
                          label=L"Δϕ", color=:blue)
    scatter_voltage = scatter!(fig_axis_voltage,
                               current_density_values, voltage_values,
                               label=L"Δϕ", color=:blue)

    power_density_values = -(current_density_values .* voltage_values)

    lines_power_density = lines!(fig_axis_power_density,
                                 current_density_values, power_density_values,
                                 label=L"PD", color=:red)
    scatter_power_density = scatter!(fig_axis_power_density,
                                     current_density_values, power_density_values,
                                     label=L"PD", color=:red)

    Legend(fig[1,2], [[lines_voltage, scatter_voltage],
                      [lines_power_density, scatter_power_density]],
            ["Voltage", "Power Density"])

    return fig
end


function reaction_terms!(f, u, node, data)
    #f .= 0.0
    region = node.region

    if occursin("1d", data.discr.spatial_discr)
        if region == data.dom[:el_neg] || region == data.dom[:el_pos]
            temp = data.study.non_isothermal ? u[data.idx[region].temp] : 1.0

            if region == data.dom[:el_neg]
                reactions = data.el_neg.reactions
            else
                reactions = data.el_pos.reactions
            end

            for r in reactions
                f[data.idx[region].ϕₛ] += integrated_volumetric_current_density(r, u,
                temp, region, data)
            end
        end
    else
        if region == data.dom[:el_neg]
            Δϕ_neg = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
            temp = data.study.non_isothermal ? u[data.idx[region].temp] : 1.0
            for r in data.el_neg.reactions
                η_l = Δϕ_neg - ϕ_eq(r, idx -> u[data.idx[region].c[idx]], temp)
                idx_c2u = data.idx[region].c
                iᵥ = volumetric_current_density(r,
                                            η_l,
                                            u[idx_c2u[r.idx_ox]], u[idx_c2u[r.idx_red]],
                                            data.boundary.v_out_neg,
                                            temp,
                                            region, data)
                f[data.idx[region].ϕₛ] += iᵥ
            end
        elseif region == data.dom[:el_pos]
            Δϕ_pos = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
            temp = data.study.non_isothermal ? u[data.idx[region].temp] : 1.0
            for r in data.el_pos.reactions
                η_r = Δϕ_pos - ϕ_eq(r, idx -> u[data.idx[region].c[idx]], temp)
                idx_c2u = data.idx[region].c
                iᵥ = volumetric_current_density(r,
                                            η_r,
                                            u[idx_c2u[r.idx_ox]], u[idx_c2u[r.idx_red]],
                                            data.boundary.v_out_pos,
                                            temp,
                                            region, data)
                f[data.idx[region].ϕₛ] += iᵥ
            end
        end
    end
end

function integrate_reaction_terms(sys, solution, data)
    integrals = integrate(sys,
                        (f, u, node, data) -> reaction_terms!(f, u, node, data),
                        solution)
    return [integrals[data.idx[data.dom[:el_neg]].ϕₛ, 2],
            integrals[data.idx[data.dom[:el_pos]].ϕₛ, 3]]
end


function plot_all_fields_1d(solution, grid, subgrids, data)
    CairoMakie.activate!()

    figures = Dict{String, Makie.Figure}()

    L0_unit = u"cm"
    L0_value = ustrip(data.scales.L0 |> L0_unit)

    coord_neg = vec(subgrids.subgrid_neg[ExtendableGrids.Coordinates]) * L0_value
    coord_pos = vec(subgrids.subgrid_pos[ExtendableGrids.Coordinates]) * L0_value
    coord_el_neg = vec(subgrids.subgrid_el_neg[ExtendableGrids.Coordinates]) * L0_value
    coord_el_pos = vec(subgrids.subgrid_el_pos[ExtendableGrids.Coordinates]) * L0_value
    coord_sep = vec(subgrids.subgrid_sep[ExtendableGrids.Coordinates]) * L0_value

    lx_sep = data.geom.lx_sep

    coord_sep .+= lx_sep/2 * L0_value
    coord_el_pos .+= lx_sep * L0_value
    coord_pos .+= lx_sep * L0_value

    P0_unit = u"kPa"
    P0_value = ustrip(data.scales.P0 |> P0_unit)

    xlim_cell = (0, cell_thickness(data.geom) * L0_value)
    ylim_cell = (0, data.geom.ly_cell * L0_value)
    axis_kwargs = Dict(:xlabel => "Through-plane direction x [$L0_unit]")

    var = data.idx
    dom_ids = data.dom

    p_neg_view = view(solution[var[dom_ids[:el_neg]].p,:], subgrids.subgrid_el_neg)
    p_pos_view = view(solution[var[dom_ids[:el_pos]].p,:], subgrids.subgrid_el_pos)

    axis_kwargs[:title] = "Pressure"
    axis_kwargs[:ylabel] = "Pressure [$P0_unit]"

    fig = line_plot_1d([[p_neg_view; p_pos_view] * P0_value],
                       [[coord_el_neg; coord_el_pos]];
                        labels=[latexstring("P")],
                        axis_kwargs=axis_kwargs)
    figures["p"] = fig


    nodes_el_neg = vec(subgrids.subgrid_el_neg[ExtendableGrids.Coordinates])
    dx_el_neg = nodes_el_neg[2:end]-nodes_el_neg[1:(end-1)]
    nodes_el_pos = vec(subgrids.subgrid_el_pos[ExtendableGrids.Coordinates]) .+ lx_sep
    dx_el_pos = nodes_el_pos[2:end]-nodes_el_pos[1:(end-1)]

    nodes_el_neg_center = similar(dx_el_neg)
    nodes_el_pos_center = similar(dx_el_pos)
    vx_neg = similar(dx_el_neg)
    vx_pos = similar(dx_el_pos)

    for idx in eachindex(dx_el_neg)
        dx = dx_el_neg[idx]
        nodes_el_neg_center[idx] = (nodes_el_neg[idx] + nodes_el_neg[idx+1]) / 2
        vx_neg[idx] = -(p_neg_view[idx+1] - p_neg_view[idx])/dx * data.el_neg.kₕ / data.el_neg.μ
    end

    for idx in eachindex(dx_el_pos)
        dx = dx_el_pos[idx]
        nodes_el_pos_center[idx] = (nodes_el_pos[idx] + nodes_el_pos[idx+1]) / 2
        vx_pos[idx] = -(p_pos_view[idx+1] - p_pos_view[idx])/dx * data.el_pos.kₕ / data.el_pos.μ
    end

    ϕₛ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₛ,:], subgrids.subgrid_neg)
    ϕₛ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₛ,:], subgrids.subgrid_pos)

    ϕₗ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₗ,:], subgrids.subgrid_el_neg)
    ϕₗ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₗ,:], subgrids.subgrid_el_pos)

    node_sep = [nodes_el_neg[end] + lx_sep/2]

    vxₕ_sep = -((p_pos_view[1] - p_neg_view[end]) / lx_sep) * data.sep.kh / data.sep.μ
    ϕₗ_neg = ϕₗ_neg_view[end]
    ϕₗ_pos = ϕₗ_pos_view[1]
    vxₑ_sep =  data.sep.kϕ / data.sep.μ * data.sep.cf * data.sep.zf * (ϕₗ_pos - ϕₗ_neg) / lx_sep / data.scaling_params.PE0
    vx_sep = [vxₕ_sep + vxₑ_sep]

    vel0_unit = u"mm/s"
    vel0_value = ustrip(data.scales.VEL0 |> vel0_unit)

    axis_kwargs[:title] = "Velocity"
    axis_kwargs[:ylabel] = "Velocity [$vel0_unit]"

    fig = line_plot_1d([[vx_neg; vx_sep; vx_pos] * vel0_value],
                       [[nodes_el_neg_center; node_sep; nodes_el_pos_center] * L0_value];
                        labels=[latexstring("v_x")],
                        axis_kwargs=axis_kwargs)
    figures["velocity_x"] = fig

    V0_unit = u"V"
    V0_value = ustrip(data.scales.V0 |> V0_unit)


    axis_kwargs[:title] = "Electrostatic Potentials "
    axis_kwargs[:ylabel] = "Electrostatic Potential [$V0_unit]"

    fig = line_plot_1d([[ϕₛ_neg_view; ϕₛ_pos_view]] * V0_value,
                       [[coord_neg; coord_pos]];
                       labels=[latexstring("\\phi_s")],
                       axis_kwargs=axis_kwargs)
    figures["phi_s"] = fig

    idx_c_neg = var[dom_ids[:el_neg]].c
    idx_c_pos = var[dom_ids[:el_pos]].c
    species_names = data.electrolyte.species.row
    num_free_species = count_free_species(data.electrolyte.species)
    num_all_species = num_free_species+1
    num_el_nodes_neg = length(coord_el_neg)
    num_el_nodes_pos = length(coord_el_pos)
    c_all_species_neg = Matrix{RealType}(undef, num_all_species, num_el_nodes_neg)
    c_all_species_pos = Matrix{RealType}(undef, num_all_species, num_el_nodes_pos)

    for idx in eachindex(idx_c_neg)
        c_all_species_neg[idx, :] = view(solution[idx_c_neg[idx],:], subgrids.subgrid_el_neg)
    end
    for idx in eachindex(idx_c_pos)
        c_all_species_pos[idx, :] = view(solution[idx_c_pos[idx],:], subgrids.subgrid_el_pos)
    end

    for idx in 1:num_el_nodes_neg
        c_all_species_neg[end, idx] = concentration_eliminated_species(c_all_species_neg[1:(end-1), idx], data)
    end
    for idx in 1:num_el_nodes_pos
        c_all_species_pos[end, idx] = concentration_eliminated_species(c_all_species_pos[1:(end-1), idx], data)
    end

    z = data.electrolyte.species[col=Key("charge")]
    cf = data.sep.cf
    zf = data.sep.zf


    ϕₗ_neg = ϕₗ_neg_view[end]
    ϕₗ_pos = ϕₗ_pos_view[1]

    if data.study.non_isothermal
        temp_neg = view(solution[var[dom_ids[:el_neg]].temp, :], subgrids.subgrid_el_neg)
        temp_pos = view(solution[var[dom_ids[:el_pos]].temp, :], subgrids.subgrid_el_pos)
        temp_sep = solution[var[dom_ids[:sep]].temp, :]
        temp_sep = temp_sep[findfirst(x->!isnan(x), temp_sep)]
    else
        temp_neg = fill(data.boundary.temp_amb, length(ϕₗ_neg_view))
        temp_pos = fill(data.boundary.temp_amb, length(ϕₗ_pos_view))
        temp_sep = data.boundary.temp_amb
    end

    ϕₗ_sep_neg = donnan_equlibrium_potential(ϕₗ_neg, c_all_species_neg[:, end], z, cf, zf, temp_neg[end])
    ϕₗ_sep_pos = donnan_equlibrium_potential(ϕₗ_pos, c_all_species_pos[:, 1], z, cf, zf, temp_pos[1])


    c_all_species_sep_neg = @MVector zeros(num_all_species)
    donnan_equilibrium_concentrations!(c_all_species_sep_neg, ϕₗ_neg, ϕₗ_sep_neg,
                                       c_all_species_neg[:, end], z, temp_sep)
    c_all_species_sep_pos = @MVector zeros(num_all_species)
    donnan_equilibrium_concentrations!(c_all_species_sep_pos, ϕₗ_pos, ϕₗ_sep_pos,
                                       c_all_species_pos[:, 1], z, temp_sep)

    axis_kwargs[:title] = "Electrostatic Potentials "
    axis_kwargs[:ylabel] = "Electrostatic Potential [$V0_unit]"

    fig = line_plot_1d([[ϕₗ_neg_view; ϕₗ_sep_neg; ϕₗ_sep_pos; ϕₗ_pos_view]] * V0_value,
                        [[coord_el_neg; coord_el_neg[end]; coord_el_pos[1]; coord_el_pos]];
                       labels=[latexstring("\\phi_l")],
                       axis_kwargs=axis_kwargs)
    figures["phi_l"] = fig

    C0_unit_str = "mol/l"
    C0_unit = uparse(C0_unit_str)
    C0_value = ustrip(data.scales.C0 |> C0_unit)

    for idx_var in 1:num_all_species
        c_l_view = c_all_species_neg[idx_var, :] * C0_value
        c_r_view = c_all_species_pos[idx_var, :] * C0_value
        c_sep_view = [c_all_species_sep_neg[idx_var]; c_all_species_sep_pos[idx_var]] * C0_value

        axis_kwargs[:title] = "Molar Species Concentrations"
        axis_kwargs[:ylabel] = "Molar Concentration [$C0_unit]"

        label = latexstring("c_{\\mathrm{$(species_names[idx_var])}}")
        fig = line_plot_1d([[c_l_view; c_sep_view; c_r_view]],
                            [[coord_el_neg; coord_el_neg[end];
                              coord_el_pos[1]; coord_el_pos]];
                           labels=[label], axis_kwargs=axis_kwargs)
        figures["c_" * species_names[idx_var]] = fig
    end

    c_views = Vector{Vector{Float64}}(undef, num_all_species)
    axis_kwargs[:title] = "Molar Species Concentrations"
    axis_kwargs[:ylabel] = "Molar Concentration [$C0_unit]"
    labels = Vector{LaTeXStrings.LaTeXString}(undef, num_all_species)

    for idx_var in 1:num_all_species
        c_views[idx_var] = [c_all_species_neg[idx_var, :];
                            c_all_species_sep_neg[idx_var];
                            c_all_species_sep_pos[idx_var];
                            c_all_species_pos[idx_var, :]]
        labels[idx_var] = latexstring("c_{\\mathrm{$(species_names[idx_var])}}")
    end
    c_views .*= C0_value

    coords = repeat([[coord_el_neg; coord_el_neg[end];
                      coord_el_pos[1]; coord_el_pos]], outer=num_all_species)

    fig = line_plot_1d(c_views, coords; labels=labels, axis_kwargs=axis_kwargs)
    figures["c"] = fig

    sigma0_unit_str = "mS/cm"
    sigma0_unit = uparse(sigma0_unit_str)
    sigma0_value = ustrip(data.scales.σ0 |> sigma0_unit)

    σₗ_views_neg = Vector{RealType}(undef, num_el_nodes_neg)
    σₗ_views_pos = Vector{RealType}(undef, num_el_nodes_pos)
    for idx in 1:num_el_nodes_neg
        σₗ_views_neg[idx] = σl_eff(c_all_species_neg[:,idx], 1.0, temp_neg[idx], data, data.dom[:el_neg]) # FIXME!!!!
    end
    for idx in 1:num_el_nodes_pos
        σₗ_views_pos[idx] = σl_eff(c_all_species_pos[:,idx], 1.0, temp_pos[idx], data, data.dom[:el_pos]) # FIXME!!!!
    end
    σₗ_sep_neg = σl_eff(c_all_species_sep_neg, 1.0, temp_sep, data, data.dom[:sep]) # FIXME!!!!
    σₗ_sep_pos = σl_eff(c_all_species_sep_pos, 1.0, temp_sep, data, data.dom[:sep]) # FIXME!!!!

    σₗ_views = [σₗ_views_neg; σₗ_sep_neg; σₗ_sep_pos; σₗ_views_pos]
    coords = [[coord_el_neg; coord_el_neg[end]; coord_el_pos[1]; coord_el_pos]]

    axis_kwargs[:title] = "Electrolyte Conductivity"
    axis_kwargs[:ylabel] = "Electrical Conductivity [$sigma0_unit]"
    label = latexstring("\\sigma_{l}")
    fig = line_plot_1d([σₗ_views] * sigma0_value, coords;
                        labels=[label], axis_kwargs=axis_kwargs)
    figures["sigma_l"] = fig

    ϕₛ_cc_neg_view = view(solution[var[dom_ids[:cc_neg]].ϕₛ,:], subgrids.subgrid_cc_neg)
    ϕₛ_cc_pos_view = view(solution[var[dom_ids[:cc_pos]].ϕₛ,:], subgrids.subgrid_cc_pos)
    ϕₛ_el_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₛ,:], subgrids.subgrid_el_neg)
    ϕₛ_el_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₛ,:], subgrids.subgrid_el_pos)
    ϕₗ_el_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₗ,:], subgrids.subgrid_el_neg)
    ϕₗ_el_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₗ,:], subgrids.subgrid_el_pos)

    num_system_var = size(solution)[1]
    solution_el_pos_view = zeros(num_system_var, length(ϕₛ_el_pos_view))
    [solution_el_pos_view[idx, :] = view(solution[idx,:], subgrids.subgrid_el_pos)
                                    for idx in 1:num_system_var]
    solution_el_neg_view = zeros(num_system_var, length(ϕₛ_el_neg_view))
    [solution_el_neg_view[idx, :] = view(solution[idx,:], subgrids.subgrid_el_neg)
                                    for idx in 1:num_system_var]

    Δϕ_neg = ϕₛ_el_neg_view - ϕₗ_el_neg_view
    iᵥ_neg = zeros(size(Δϕ_neg))
    η_neg = similar(Δϕ_neg)
    r_neg = data.el_neg.reactions
    for idx_node in eachindex(Δϕ_neg)
        η_neg[idx_node] = Δϕ_neg[idx_node] -
                    ϕ_eq(r_neg[1], idx -> c_all_species_neg[idx, idx_node], temp_neg[idx_node]) # FIXME: Handle multiple reactions
        for r in r_neg
            iᵥ_neg[idx_node] += integrated_volumetric_current_density(r, solution_el_neg_view[:, idx_node],
                                                temp_neg[idx_node], data.dom[:el_neg], data)
        end
    end
    Δϕ_pos = ϕₛ_el_pos_view - ϕₗ_el_pos_view
    iᵥ_pos = zeros(size(Δϕ_pos))
    η_pos = similar(Δϕ_pos)
    r_pos = data.el_pos.reactions
    for idx_node in eachindex(Δϕ_pos)
        η_pos[idx_node] = Δϕ_pos[idx_node] -
                        ϕ_eq(r_pos[1], idx -> c_all_species_pos[idx, idx_node], temp_pos[idx_node]) # FIXME: Handle multiple reactions
        for r in r_pos
            iᵥ_pos[idx_node] += integrated_volumetric_current_density(r, solution_el_pos_view[:, idx_node],
                                                temp_pos[idx_node], data.dom[:el_pos], data)
        end

    end

    axis_kwargs[:title] = "Overpotential"
    axis_kwargs[:ylabel] = "Electrical Overpotential [$V0_unit]"
    label = latexstring("\\eta_l")
    fig = line_plot_1d([[η_neg; η_pos]] * V0_value, [[coord_el_neg; coord_el_pos]];
                        labels=[label], axis_kwargs=axis_kwargs)
    figures["eta_l"] = fig

    iv0_unit = u"mA/mm^3"
    iv0_value = ustrip(data.scales.iv0 |> iv0_unit)

    axis_kwargs[:title] = "Volumetric Current Density"
    axis_kwargs[:ylabel] = "Volumetric Current Density [$iv0_unit]"
    label = latexstring("i_v")
    fig = line_plot_1d([[iᵥ_neg; iᵥ_pos]] * iv0_value, [[coord_el_neg; coord_el_pos]];
                        labels=[label], axis_kwargs=axis_kwargs)
    figures["iv"] = fig

    if data.study.non_isothermal
        T0_unit = u"K"
        T0_unit_str = "K"
        T0_value = ustrip(data.scales.TEMP0 |> T0_unit)

        temp_l_view = view(solution[var[dom_ids[:el_neg]].temp,:], subgrids.subgrid_neg) * T0_value
        temp_r_view = view(solution[var[dom_ids[:el_pos]].temp,:], subgrids.subgrid_pos) * T0_value

        temp_i_sol = solution[var[dom_ids[:sep]].temp,:] * T0_value
        temp_i_sol = temp_i_sol[temp_i_sol .> 0]

        axis_kwargs[:title] = "Temperature"
        axis_kwargs[:ylabel] = "Temperature [$T0_unit]"
        label = latexstring("T\\;[\\mathrm{$(T0_unit_str)}]")

        fig = line_plot_1d([[temp_l_view; temp_i_sol; temp_r_view]],
                           [[coord_neg; coord_sep; coord_pos]];
                            labels=[label], axis_kwargs=axis_kwargs)
        figures["temp"] = fig

        s_η_neg = iᵥ_neg .* η_neg
        s_η_pos = iᵥ_pos .* η_pos

        r_neg = data.el_neg.reactions[1]
        s_r_neg = iᵥ_neg .* temp_neg * r_neg.Δs / r_neg.ν_el # FIXME: Allow multiple reactions
        r_pos = data.el_pos.reactions[1]
        s_r_pos = iᵥ_pos .* temp_pos * r_pos.Δs / r_pos.ν_el # FIXME: Allow multiple reactions

        nodes_cc_neg = vec(subgrids.subgrid_cc_neg[ExtendableGrids.Coordinates])
        dx_cc_neg = nodes_cc_neg[2:end]-nodes_cc_neg[1:(end-1)]
        nodes_cc_pos = vec(subgrids.subgrid_cc_pos[ExtendableGrids.Coordinates]) .+ data.geom.lx_sep
        dx_cc_pos = nodes_cc_pos[2:end]-nodes_cc_pos[1:(end-1)]
        nodes_el_neg = vec(subgrids.subgrid_el_neg[ExtendableGrids.Coordinates])
        dx_el_neg = nodes_el_neg[2:end]-nodes_el_neg[1:(end-1)]
        nodes_el_pos = vec(subgrids.subgrid_el_pos[ExtendableGrids.Coordinates]) .+ data.geom.lx_sep
        dx_el_pos = nodes_el_pos[2:end]-nodes_el_pos[1:(end-1)]

        s_j_cc_neg = similar(dx_cc_neg)
        s_j_cc_pos = similar(dx_cc_pos)
        s_j_el_neg = similar(dx_el_neg)
        s_j_el_pos = similar(dx_el_pos)

        nodes_cc_neg_center = similar(dx_cc_neg)
        nodes_cc_pos_center = similar(dx_cc_pos)
        nodes_el_neg_center = similar(dx_el_neg)
        nodes_el_pos_center = similar(dx_el_pos)

        σₛ_cc_neg = data.cc_neg.σₑ
        for idx in eachindex(dx_cc_neg)
            dx = dx_cc_neg[idx]
            s_j_cc_neg[idx] = σₛ_cc_neg * ((ϕₛ_cc_neg_view[idx+1] - ϕₛ_cc_neg_view[idx]) / dx).^2
            nodes_cc_neg_center[idx] = (nodes_cc_neg[idx] + nodes_cc_neg[idx+1]) / 2
        end

        σₛ_cc_pos = data.cc_pos.σₑ
        for idx in eachindex(dx_cc_neg)
            dx = dx_cc_pos[idx]
            s_j_cc_pos[idx] = σₛ_cc_pos * ((ϕₛ_cc_pos_view[idx+1] - ϕₛ_cc_pos_view[idx]) / dx).^2
            nodes_cc_pos_center[idx] = (nodes_cc_pos[idx] + nodes_cc_pos[idx+1]) / 2
        end

        for idx in eachindex(dx_el_neg)
            σₗ_center = (σₗ_views_neg[idx] + σₗ_views_neg[idx+1]) / 2
            dx = dx_el_neg[idx]
            s_j_el_neg[idx] = σₗ_center * ((ϕₗ_el_neg_view[idx+1] - ϕₗ_el_neg_view[idx]) / dx).^2
            nodes_el_neg_center[idx] = (nodes_el_neg[idx] + nodes_el_neg[idx+1]) / 2
        end

        for idx in eachindex(dx_el_pos)
            σₗ_center = (σₗ_views_pos[idx] + σₗ_views_pos[idx+1]) / 2
            dx = dx_el_pos[idx]
            s_j_el_pos[idx] = σₗ_center * ((ϕₗ_el_pos_view[idx+1] - ϕₗ_el_pos_view[idx]) / dx).^2
            nodes_el_pos_center[idx] = (nodes_el_pos[idx] + nodes_el_pos[idx+1]) / 2
        end

        s_j_sep_neg = σₗ_sep_neg * ((ϕₗ_sep_pos - ϕₗ_sep_neg)/data.geom.lx_sep)^2
        s_j_sep_pos = σₗ_sep_pos * ((ϕₗ_sep_pos - ϕₗ_sep_neg)/data.geom.lx_sep)^2

        h0_unit_str = "kW/m^3" # this is equivalent to "mJ/(cm^3*s)"
        h0_unit = uparse(h0_unit_str)
        V0_iv0_value = ustrip(data.scales.V0 * data.scales.iv0  |> h0_unit)

        axis_kwargs[:title] = "Heat Production / Consumption"
        axis_kwargs[:ylabel] = "Energy source [$h0_unit]"
        labels = [latexstring("S^{(\\eta)}"), latexstring("S^{(r)}"), latexstring("S^{(j)}")]

        nodes_cc_neg_center *= L0_value
        nodes_el_neg_center *= L0_value
        nodes_el_pos_center *= L0_value
        nodes_cc_pos_center *= L0_value

        ϵL0 = data.scaling_params.ϵL0
        heat_production_values = [[s_η_neg; s_η_pos],
                                  [s_r_neg; s_r_pos],
                                 ([s_j_cc_neg; s_j_el_neg; s_j_sep_neg;
                                   s_j_sep_pos; s_j_el_pos; s_j_cc_pos] * ϵL0^2) ]
        heat_production_values *= V0_iv0_value
        fig = line_plot_1d(heat_production_values,
                            [[coord_el_neg; coord_el_pos],
                            [coord_el_neg; coord_el_pos],
                            [nodes_cc_neg_center; nodes_el_neg_center; coord_el_neg[end];
                            coord_el_pos[1]; nodes_el_pos_center; nodes_cc_pos_center]];
                            labels=labels, axis_kwargs=axis_kwargs)
        figures["heat_production"] = fig
    end
    return figures
end

function line_plot_1d(sol, coords; labels=[""], axis_kwargs=axis_kwargs)
    fig = Makie.Figure()
    min_sol = Float32(minimum(Iterators.flatten(sol)))
    max_sol = Float32(maximum(Iterators.flatten(sol)))

    δ = 1e-4
    limits = nothing
    if isapprox(min_sol, max_sol, rtol=1e-4)
        limits = (min_sol * (1 - δ) - δ, max_sol * (1 + δ) + δ)
    end
    fig_axis = Makie.Axis(fig[1,1]; collect(axis_kwargs)..., limits=(nothing, limits))

    vec_l = []
    vec_s = []
    for idx in eachindex(coords)
        l = lines!(fig_axis, coords[idx], sol[idx], label=labels[idx])
        s = scatter!(fig_axis, coords[idx], sol[idx], label=labels[idx])

        push!(vec_l, l)
        push!(vec_s, s)
    end

    Legend(fig[1,2], vec_l, labels)

    return fig
end

function plot_all_fields_2d(solution, grid, subgrids, data)
    CairoMakie.activate!()

    dict_figures = Dict{String, Makie.Figure}()

    L0_unit = u"cm"
    L0_value = ustrip(data.scales.L0 |> L0_unit)

    coord_neg = subgrids.subgrid_neg[ExtendableGrids.Coordinates] * L0_value
    coord_pos = subgrids.subgrid_pos[ExtendableGrids.Coordinates] * L0_value
    coord_el_neg = subgrids.subgrid_el_neg[ExtendableGrids.Coordinates] * L0_value
    coord_el_pos = subgrids.subgrid_el_pos[ExtendableGrids.Coordinates] * L0_value
    coord_sep = subgrids.subgrid_sep[ExtendableGrids.Coordinates] * L0_value

    coord_sep .+= [data.geom.lx_sep/2, 0.0] * L0_value
    coord_el_pos .+= [data.geom.lx_sep, 0.0] * L0_value
    coord_pos .+= [data.geom.lx_sep, 0.0] * L0_value

    xlim_cell = (0, cell_thickness(data.geom) * L0_value)
    ylim_cell = (0, data.geom.ly_cell * L0_value)
    axis_kwargs = Dict(:xlabel => "Through-plane direction x [$L0_unit]",
                       :ylabel => "In-plane direction y [$L0_unit]",
                       :limits => ((nothing, nothing), ylim_cell))
    axis_kwargs_cell = Dict(:xlabel => "Through-plane direction x [$L0_unit]",
                       :ylabel => "In-plane direction y [$L0_unit]",
                       :limits => (xlim_cell, ylim_cell))

    P0_unit = u"kPa"
    P0_value = ustrip(data.scales.P0 |> P0_unit)
    var = data.idx
    dom_ids = data.dom

    p_neg_view = view(solution[var[dom_ids[:el_neg]].p,:], subgrids.subgrid_el_neg) * P0_value
    p_pos_view = view(solution[var[dom_ids[:el_pos]].p,:], subgrids.subgrid_el_pos) * P0_value

    axis_kwargs_cell[:title] = L"\text{Pressure}"
    fig = contour_plot_2d([p_neg_view, p_pos_view],
                          [coord_el_neg, coord_el_pos];
                          label=latexstring("P\\;[\\mathrm{$(P0_unit)}]"),
                          axis_kwargs=axis_kwargs_cell)
    dict_figures["p"] = fig

    V0_unit = u"V"
    V0_value = ustrip(data.scales.V0 |> V0_unit)
    ϕₛ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₛ,:], subgrids.subgrid_neg) * V0_value
    ϕₛ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₛ,:], subgrids.subgrid_pos) * V0_value
    axis_kwargs_cell[:title] = L"\text{Electrostatic Potential} \; \phi_s"
    fig = contour_plot_2d([ϕₛ_neg_view, ϕₛ_pos_view],
                          [coord_neg, coord_pos];
                          label=latexstring("\\phi_s\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs_cell)
    dict_figures["ϕₛ"] = fig

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_s \; \text{in Negative Electrode}"
    fig = contour_plot_2d([ϕₛ_neg_view],
                          [coord_neg];
                          label=latexstring("\\phi_s\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["phi_s_neg"] = fig

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_s \; \text{in Positive Electrode}"
    fig = contour_plot_2d([ϕₛ_pos_view],
                          [coord_pos];
                          label=latexstring("\\phi_s\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["phi_s_pos"] = fig

    ϕₗ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₗ,:], subgrids.subgrid_el_neg) * V0_value
    ϕₗ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₗ,:], subgrids.subgrid_el_pos) * V0_value
    axis_kwargs_cell[:title] = L"\text{Electrostatic Potential} \; \phi_l"
    fig = contour_plot_2d([ϕₗ_neg_view, ϕₗ_pos_view],
                          [coord_el_neg, coord_el_pos];
                          label=latexstring("\\phi_l\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs_cell)
    dict_figures["phi_l"] = fig

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_l \; \text{in Negative Electrode}"
    fig = contour_plot_2d([ϕₗ_neg_view],
                          [coord_el_neg];
                          label=latexstring("\\phi_l\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["phi_l_neg"] = fig

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_l \; \text{in Positive Electrode}"
    fig = contour_plot_2d([ϕₗ_pos_view],
                          [coord_el_pos];
                          label=latexstring("\\phi_l\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["phi_l_pos"] = fig

    C0_unit_str = "mol/l"
    C0_unit = uparse(C0_unit_str)
    C0_value = ustrip(data.scales.C0 |> C0_unit)

    idx_c_neg = var[dom_ids[:el_neg]].c
    idx_c_pos = var[dom_ids[:el_pos]].c
    species_names = data.electrolyte.species.row[1:(end-1)]
    for idx_var in eachindex(species_names)
        c_l_view = view(solution[idx_c_neg[idx_var],:], subgrids.subgrid_el_neg) * C0_value
        c_r_view = view(solution[idx_c_pos[idx_var],:], subgrids.subgrid_el_pos) * C0_value
        title = L"\text{Molar Concentration} \; c_{\mathrm{%$(species_names[idx_var])}}"
        axis_kwargs_cell[:title] = title;
        label = latexstring("c_{\\mathrm{$(species_names[idx_var])}}\\;[\\mathrm{$C0_unit_str}]")
        fig = contour_plot_2d([c_l_view, c_r_view],
                              [coord_el_neg, coord_el_pos];
                               label=label,
                               axis_kwargs=axis_kwargs_cell)
        plot_name = "c_" * species_names[idx_var]
        dict_figures[plot_name] = fig

        title = L"\text{Molar Concentration} \; c_{\mathrm{%$(species_names[idx_var])}} \;
                  \text{in Negative Electrode}"
        axis_kwargs[:title] = title
        fig = contour_plot_2d([c_l_view],
                              [coord_el_neg];
                              label=label,
                              axis_kwargs=axis_kwargs)
        plot_name = "c_" * species_names[idx_var] * "_neg"
        dict_figures[plot_name] = fig

        title = L"\text{Molar Concentration} \; c_{\mathrm{%$(species_names[idx_var])}} \;
                  \text{in Positive Electrode}"
        axis_kwargs[:title] = title
        fig = contour_plot_2d([c_r_view],
                              [coord_el_pos];
                              label=label,
                              axis_kwargs=axis_kwargs)
        plot_name = "c_" * species_names[idx_var] * "_pos"
        dict_figures[plot_name] = fig
    end

    ϕₛ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₛ,:], subgrids.subgrid_el_neg)
    ϕₗ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₗ,:], subgrids.subgrid_el_neg)
    ϕₛ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₛ,:], subgrids.subgrid_el_pos)
    ϕₗ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₗ,:], subgrids.subgrid_el_pos)

    c_neg_views = [view(solution[idx,:], subgrids.subgrid_el_neg)
                        for idx in var[dom_ids[:el_neg]].c]
    c_pos_views = [view(solution[idx,:], subgrids.subgrid_el_pos)
                        for idx in var[dom_ids[:el_pos]].c]

    Δϕ_neg = ϕₛ_neg_view - ϕₗ_neg_view
    Δϕ_pos = ϕₛ_pos_view - ϕₗ_pos_view
    η_l = similar(Δϕ_neg)
    η_r = similar(Δϕ_pos)
    iᵥ_l = similar(Δϕ_neg)
    iᵥ_r = similar(Δϕ_neg)

    T0_unit = u"K"
    T0_unit_str = "K"
    T0_value = ustrip(data.scales.TEMP0 |> T0_unit)

    if data.study.non_isothermal
        temp_l_view = view(solution[var[dom_ids[:el_neg]].temp,:], subgrids.subgrid_el_neg)
        temp_r_view = view(solution[var[dom_ids[:el_pos]].temp,:], subgrids.subgrid_el_pos)
    else
        temp_l_view = ones(size(Δϕ_neg))
        temp_r_view = ones(size(Δϕ_pos))
    end

    for idx in eachindex(η_l)
        c_neg = [c_neg_views[c_idx][idx] for c_idx in eachindex(c_neg_views)]
        r = data.el_neg.reactions[1] # FIXME

        c_neg_ox = c_neg[r.idx_ox]
        c_neg_red = c_neg[r.idx_red]

        η_l[idx] = Δϕ_neg[idx] - ϕ_eq(data.el_neg.reactions,
                                      idx -> c_neg[idx],
                                      temp_l_view[idx])
        iᵥ_l[idx] = volumetric_current_density(r,
                                                η_l[idx],
                                                c_neg_ox, c_neg_red,
                                                data.boundary.v_out_neg,
                                                temp_l_view[idx],
                                                dom_ids[:el_neg],
                                                data)
    end
    for idx in eachindex(η_r)
        c_pos = [c_pos_views[c_idx][idx] for c_idx in eachindex(c_pos_views)]
        r = data.el_pos.reactions[1] # FIXME
        c_pos_ox = c_pos[r.idx_ox]
        c_pos_red = c_pos[r.idx_red]

        η_r[idx] = Δϕ_pos[idx] - ϕ_eq(data.el_pos.reactions,
                                      idx -> c_pos[idx],
                                      temp_r_view[idx])
        iᵥ_r[idx] = volumetric_current_density(r,
                                                η_r[idx],
                                                c_pos_ox, c_pos_red,
                                                data.boundary.v_out_pos,
                                                temp_r_view[idx],
                                                dom_ids[:el_pos],
                                                data)
    end

    η_l .*= V0_value
    η_r .*= V0_value

    iv0_unit = u"A/cm^3"
    iv0_value = ustrip(data.scales.iv0 |> iv0_unit)
    iᵥ_l .*= iv0_value
    iᵥ_r .*= iv0_value

    # TODO: Reactivate after implementation of η and iᵥ
    title = L"\text{Overpotential}"
    axis_kwargs[:title] = title;
    fig = contour_plot_2d([η_l], [coord_el_neg];
                          label=latexstring("\\eta\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["eta_neg"] = fig
    fig = contour_plot_2d([η_r], [coord_el_pos];
                          label=latexstring("\\eta\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["eta_pos"] = fig

    title = L"\text{Volumetric Reaction Source}"
    axis_kwargs[:title] = title;
    fig = contour_plot_2d([iᵥ_l], [coord_el_neg];
                          label=latexstring("i_v\\;[\\mathrm{$(iv0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["iv_neg"] = fig
    fig = contour_plot_2d([iᵥ_r], [coord_el_pos];
                          label=latexstring("i_v\\;[\\mathrm{$(iv0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    dict_figures["iv_pos"] = fig
    # ####

    if data.study.non_isothermal
        temp_l_view = view(solution[var[dom_ids[:el_neg]].temp,:], subgrids.subgrid_neg) * T0_value
        temp_r_view = view(solution[var[dom_ids[:el_pos]].temp,:], subgrids.subgrid_pos) * T0_value

        temp_i_sol = solution[var[dom_ids[:sep]].temp,:] * T0_value
        temp_i_sol = temp_i_sol[temp_i_sol .!== NaN]

        axis_kwargs_cell[:title] = L"\text{Temperature}"

        ϵ = 1e-8
        coord_sep_l = coord_sep .- [data.geom.lx_sep/2 * (1 - ϵ), 0.0] * L0_value
        coord_sep_r = coord_sep .+ [data.geom.lx_sep/2 * (1 - ϵ), 0.0] * L0_value

        fig = contour_plot_2d([vcat(temp_l_view,
                                    temp_i_sol, temp_i_sol,
                                    temp_r_view)],
                              [hcat(coord_neg,
                                    coord_sep_l, coord_sep_r,
                                    coord_pos)];
                              label=latexstring("T\\;[\\mathrm{$(T0_unit_str)}]"),
                              axis_kwargs=axis_kwargs_cell)
        dict_figures["temp"] = fig
    end
    return dict_figures
end

function contour_plot_2d(solutions, coords;
                        label="",
                        levels=30,
                        axis_kwargs=axis_kwargs)
    fig = Makie.Figure()
    fig_axis = Makie.Axis(fig[1,1];
                          backgroundcolor=:gray,
                          xgridvisible = false,
                          ygridvisible = false,
                          axis_kwargs...)
    if axis_kwargs[:limits][1][1] === nothing
        tightlimits!(fig_axis, Left())
    end
    if axis_kwargs[:limits][1][2] === nothing
        tightlimits!(fig_axis, Right())
    end
    if axis_kwargs[:limits][1][1] === nothing || axis_kwargs[:limits][1][2] === nothing
        tightlimits!(fig_axis, Bottom())
        tightlimits!(fig_axis, Top())
    end

    data_range = (minimum(minimum.(solutions)), maximum(maximum.(solutions)))

    # in case of constant data we need to slightly enlarge the colorrange
    # to avoid an issue in Makie
    δ_tol = 1e-4
    if isapprox(data_range[1], data_range[2]; rtol=δ_tol, atol=δ_tol)
        data_range = (data_range[1] * (1 - δ_tol) - δ_tol,
                      data_range[2] * (1 + δ_tol) + δ_tol)
    end
    # choose the levels according to the overall data limits
    levels = range(data_range[1], stop=data_range[2], length=levels)

    for idx in eachindex(coords)
        contourf!(fig_axis,
                vec(coords[idx][1,:]), vec(coords[idx][2,:]), vec(solutions[idx]);
                levels=levels,
                colorrange=data_range,
                extendlow=:auto, extendhigh=:auto,
                colormap=:viridis)
    end
    ticks = collect(range(data_range..., length=5))

    ticks[1] = round(ticks[1], sigdigits=6, RoundUp)
    ticks[end] = round(ticks[end], sigdigits=6, RoundDown)
    ticks[2:(end-1)] .= round.(ticks[2:(end-1)], sigdigits=6, RoundNearest)

    Colorbar(fig[:, end+1],
            colorrange=data_range,
            label=label,
            labelrotation=0,
            ticks=ticks)
    return fig
end

function plot_grid(grid; plotter=CairoMakie, size=size)
    vis = GridVisualizer(Plotter=plotter, layout=(1,1); size=(800,800), show=false, reveal=false)
    gridplot!(vis, grid; show=false, reveal=false)
    return vis
end
