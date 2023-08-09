function plot_polarization_curve(current_density_values, voltage_values)
    CairoMakie.activate!()
    fontsize_theme = Theme(fontsize = 22)
    set_theme!(fontsize_theme)

    fig = Makie.Figure()
    fig_axis_voltage = Makie.Axis(fig[1,1];
                                  title="Polarization Curve",
                                  xlabel="Current Density [mA/cm^2]",
                                  ylabel="Voltage [V]",
                                  yticklabelcolor = :blue)
    fig_axis_voltage.xreversed = true
    fig_axis_power_density = Makie.Axis(fig[1,1];
                                        yticklabelcolor = :red,
                                        ylabel="Power Density [mW/cm^2]",
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

    power_density_values = abs.(current_density_values .* voltage_values)

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
    f .= 0.0
    region = node.region
    if region == data.dom[:el_neg]
        r = data.el_neg.reactions
        Δϕ_neg = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
        temp = data.study.non_isothermal ? u[data.idx[region].temp] : 1.0
        η_l = Δϕ_neg - ϕ_eq(r, u, data.idx[region].c, temp)
        idx_c2u = data.idx[region].c
        iᵥ = volumetric_current_density(r,
                                       η_l,
                                       u[idx_c2u[r.idx_ox]], u[idx_c2u[r.idx_red]],
                                       data.boundary.v_out_neg,
                                       temp,
                                       region, data)
        f[data.idx[region].ϕₛ] = iᵥ
    elseif region == data.dom[:el_pos]
        r = data.el_pos.reactions
        Δϕ_pos = u[data.idx[region].ϕₛ] - u[data.idx[region].ϕₗ]
        temp = data.study.non_isothermal ? u[data.idx[region].temp] : 1.0
        η_r = Δϕ_pos - ϕ_eq(r, u, data.idx[region].c, temp)

        idx_c2u = data.idx[region].c
        iᵥ = volumetric_current_density(r,
                                       η_r,
                                       u[idx_c2u[r.idx_ox]], u[idx_c2u[r.idx_red]],
                                       data.boundary.v_out_pos,
                                       temp,
                                       region, data)

        f[data.idx[region].ϕₛ] = iᵥ
    end
end

function integrate_reaction_terms(sys, solution, data)
    integrals = integrate(sys,
                        (f, u, node, data) -> reaction_terms!(f, u, node, data),
                        solution)
    return [integrals[data.idx[data.dom[:el_neg]].ϕₛ, 2],
            integrals[data.idx[data.dom[:el_pos]].ϕₛ, 3]]
end


# TODO: FIX USAGE OF KINETIC / THERMODYNAMIC FUNCTIONS!
function plot_all_fields_2d(solution, grid, subgrids, prms, dir)
    CairoMakie.activate!()

    L0_unit = u"cm"
    L0_value = ustrip(prms.scales.L0 |> L0_unit)

    coord_neg = subgrids.subgrid_neg[ExtendableGrids.Coordinates] * L0_value
    coord_pos = subgrids.subgrid_pos[ExtendableGrids.Coordinates] * L0_value
    coord_el_neg = subgrids.subgrid_el_neg[ExtendableGrids.Coordinates] * L0_value
    coord_el_pos = subgrids.subgrid_el_pos[ExtendableGrids.Coordinates] * L0_value
    coord_sep = subgrids.subgrid_sep[ExtendableGrids.Coordinates] * L0_value

    coord_sep .+= [prms.geom.lx_sep/2, 0.0] * L0_value
    coord_el_pos .+= [prms.geom.lx_sep, 0.0] * L0_value
    coord_pos .+= [prms.geom.lx_sep, 0.0] * L0_value

    dir_path = pwd()
    output_dir_path = joinpath(dir_path, prms.study.output_folder, dir, "2d")
    mkpath(output_dir_path)

    xlim_cell = (0, cell_thickness(prms.geom) * L0_value)
    ylim_cell = (0, prms.geom.ly_cell * L0_value)
    axis_kwargs = Dict(:xlabel => "Through-plane direction x [$L0_unit]",
                       :ylabel => "In-plane direction y [$L0_unit]",
                       :limits => ((nothing, nothing), ylim_cell))
    axis_kwargs_cell = Dict(:xlabel => "Through-plane direction x [$L0_unit]",
                       :ylabel => "In-plane direction y [$L0_unit]",
                       :limits => (xlim_cell, ylim_cell))

    P0_unit = u"kPa"
    P0_value = ustrip(prms.scales.P0 |> P0_unit)
    var = prms.idx
    dom_ids = prms.dom

    p_neg_view = view(solution[var[dom_ids[:el_neg]].p,:], subgrids.subgrid_el_neg) * P0_value
    p_pos_view = view(solution[var[dom_ids[:el_pos]].p,:], subgrids.subgrid_el_pos) * P0_value

    axis_kwargs_cell[:title] = L"\text{Pressure}"
    fig = contour_plot_2d([p_neg_view, p_pos_view],
                          [coord_el_neg, coord_el_pos];
                          label=latexstring("P\\;[\\mathrm{$(P0_unit)}]"),
                          axis_kwargs=axis_kwargs_cell)
    Makie.save(joinpath(output_dir_path, "p.pdf"), fig)

    V0_unit = u"V"
    V0_value = ustrip(prms.scales.V0 |> V0_unit)
    phi_s_l_view = view(solution[var[dom_ids[:el_neg]].ϕₛ,:], subgrids.subgrid_neg) * V0_value
    phi_s_r_view = view(solution[var[dom_ids[:el_pos]].ϕₛ,:], subgrids.subgrid_pos) * V0_value
    axis_kwargs_cell[:title] = L"\text{Electrostatic Potential} \; \phi_s"
    fig = contour_plot_2d([phi_s_l_view, phi_s_r_view],
                          [coord_neg, coord_pos];
                          label=latexstring("\\phi_s\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs_cell)
    Makie.save(joinpath(output_dir_path, "phi_s.pdf"), fig)

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_s \; \text{in Negative Electrode}"
    fig = contour_plot_2d([phi_s_l_view],
                          [coord_neg];
                          label=latexstring("\\phi_s\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "phi_s_neg.pdf"), fig)

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_s \; \text{in Positive Electrode}"
    fig = contour_plot_2d([phi_s_r_view],
                          [coord_pos];
                          label=latexstring("\\phi_s\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "phi_s_pos.pdf"), fig)

    phi_l_l_view = view(solution[var[dom_ids[:el_neg]].ϕₗ,:], subgrids.subgrid_el_neg) * V0_value
    phi_l_r_view = view(solution[var[dom_ids[:el_pos]].ϕₗ,:], subgrids.subgrid_el_pos) * V0_value
    axis_kwargs_cell[:title] = L"\text{Electrostatic Potential} \; \phi_l"
    fig = contour_plot_2d([phi_l_l_view, phi_l_r_view],
                          [coord_el_neg, coord_el_pos];
                          label=latexstring("\\phi_l\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs_cell)
    Makie.save(joinpath(output_dir_path, "phi_l.pdf"), fig)

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_l \; \text{in Negative Electrode}"
    fig = contour_plot_2d([phi_l_l_view],
                          [coord_el_neg];
                          label=latexstring("\\phi_l\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "phi_l_neg.pdf"), fig)

    axis_kwargs[:title] = L"\text{Electrostatic Potential} \; \phi_l \; \text{in Positive Electrode}"
    fig = contour_plot_2d([phi_l_r_view],
                          [coord_el_pos];
                          label=latexstring("\\phi_l\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "phi_l_pos.pdf"), fig)


    C0_unit_str = "mol/l"
    C0_unit = uparse(C0_unit_str)
    C0_value = ustrip(prms.scales.C0 |> C0_unit)

    idx_c_neg = var[dom_ids[:el_neg]].c
    idx_c_pos = var[dom_ids[:el_pos]].c
    species_names = prms.electrolyte.species.row
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
        file_name = "c_" * species_names[idx_var] * ".pdf"
        Makie.save(joinpath(output_dir_path, file_name), fig)

        title = L"\text{Molar Concentration} \; c_{\mathrm{%$(species_names[idx_var])}} \;
                  \text{in Negative Electrode}"
        axis_kwargs[:title] = title
        fig = contour_plot_2d([c_l_view],
                              [coord_el_neg];
                              label=label,
                              axis_kwargs=axis_kwargs)
        file_name = "c_" * species_names[idx_var] * "_neg.pdf"
        Makie.save(joinpath(output_dir_path, file_name), fig)

        title = L"\text{Molar Concentration} \; c_{\mathrm{%$(species_names[idx_var])}} \;
                  \text{in Positive Electrode}"
        axis_kwargs[:title] = title
        fig = contour_plot_2d([c_r_view],
                              [coord_el_pos];
                              label=label,
                              axis_kwargs=axis_kwargs)
        file_name = "c_" * species_names[idx_var] * "_pos.pdf"
        Makie.save(joinpath(output_dir_path, file_name), fig)
    end

    ####
    ϕₛ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₛ,:], subgrids.subgrid_el_neg)
    ϕₗ_neg_view = view(solution[var[dom_ids[:el_neg]].ϕₗ,:], subgrids.subgrid_el_neg)
    ϕₛ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₛ,:], subgrids.subgrid_el_pos)
    ϕₗ_pos_view = view(solution[var[dom_ids[:el_pos]].ϕₗ,:], subgrids.subgrid_el_pos)

    #FIXME: Determine indices from reactions!
    c_neg_views = [view(solution[idx,:], subgrids.subgrid_el_neg)
                        for idx in var[dom_ids[:el_neg]].c]
    c_pos_views = [view(solution[idx,:], subgrids.subgrid_el_pos)
                        for idx in var[dom_ids[:el_pos]].c]

    # c_ox_neg_l_view = view(solution[var[dom_ids[:el_neg]].c[1],:], subgrids.subgrid_el_neg)
    # c_red_neg_l_view = view(solution[var[dom_ids[:el_neg]].c[2],:], subgrids.subgrid_el_neg)
    # c_ox_pos_r_view = view(solution[var[dom_ids[:el_pos]].c[3],:], subgrids.subgrid_el_pos)
    # c_red_pos_r_view = view(solution[var[dom_ids[:el_pos]].c[4],:], subgrids.subgrid_el_pos)

    Δϕ_neg = ϕₛ_neg_view - ϕₗ_neg_view
    Δϕ_pos = ϕₛ_pos_view - ϕₗ_pos_view
    η_l = similar(Δϕ_neg)
    η_r = similar(Δϕ_pos)
    iᵥ_l = similar(Δϕ_neg)
    iᵥ_r = similar(Δϕ_neg)

    T0_unit = u"K"
    T0_unit_str = "K"
    T0_value = ustrip(prms.scales.TEMP0 |> T0_unit)

    if prms.study.non_isothermal
        temp_l_view = view(solution[var[dom_ids[:el_neg]].temp,:], subgrids.subgrid_el_neg)
        temp_r_view = view(solution[var[dom_ids[:el_pos]].temp,:], subgrids.subgrid_el_pos)
    else
        temp_l_view = ones(size(Δϕ_neg))
        temp_r_view = ones(size(Δϕ_pos))
    end

    # TODO: Evaluate these in a separate function based on the nondimensional quantities!
    for idx in eachindex(η_l)
        c_neg = [c_neg_views[c_idx][idx] for c_idx in eachindex(c_neg_views)]
        r = prms.el_neg.reactions

        c_neg_ox = c_neg[r.idx_ox]
        c_neg_red = c_neg[r.idx_red]

        η_l[idx] = Δϕ_neg[idx] - ϕ_eq(prms.el_neg.reactions, c_neg, 1:length(c_neg),
                                      temp_l_view[idx])
        iᵥ_l[idx] = volumetric_current_density(prms.el_neg.reactions,
                                                η_l[idx],
                                                c_neg_ox, c_neg_red,
                                                prms.boundary.v_out_neg, # TODO: use local velocity
                                                temp_l_view[idx],
                                                dom_ids[:el_neg],
                                                prms)
    end
    for idx in eachindex(η_r)
        c_pos = [c_pos_views[c_idx][idx] for c_idx in eachindex(c_pos_views)]
        r = prms.el_pos.reactions
        c_pos_ox = c_pos[r.idx_ox]
        c_pos_red = c_pos[r.idx_red]

        η_r[idx] = Δϕ_pos[idx] - ϕ_eq(prms.el_pos.reactions, c_pos, 1:length(c_pos),
                                      temp_r_view[idx])
        iᵥ_r[idx] = volumetric_current_density(prms.el_neg.reactions,
                                                η_r[idx],
                                                c_pos_ox, c_pos_red,
                                                prms.boundary.v_out_pos, # TODO: use local velocity
                                                temp_r_view[idx],
                                                dom_ids[:el_pos],
                                                prms)
    end

    η_l .*= V0_value
    η_r .*= V0_value

    iv0_unit = u"A/cm^3"
    iv0_value = ustrip(prms.scales.iv0 |> iv0_unit)
    iᵥ_l .*= iv0_value
    iᵥ_r .*= iv0_value

    # TODO: Reactivate after implementation of η and iᵥ
    title = L"\text{Overpotential}"
    axis_kwargs[:title] = title;
    fig = contour_plot_2d([η_l], [coord_el_neg];
                          label=latexstring("\\eta\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "eta_neg.pdf"), fig)
    fig = contour_plot_2d([η_r], [coord_el_pos];
                          label=latexstring("\\eta\\;[\\mathrm{$(V0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "eta_pos.pdf"), fig)

    title = L"\text{Volumetric Reaction Source}"
    axis_kwargs[:title] = title;
    fig = contour_plot_2d([iᵥ_l], [coord_el_neg];
                          label=latexstring("i_v\\;[\\mathrm{$(iv0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "i_v_neg.pdf"), fig)
    fig = contour_plot_2d([iᵥ_r], [coord_el_pos];
                          label=latexstring("i_v\\;[\\mathrm{$(iv0_unit)}]"),
                          axis_kwargs=axis_kwargs)
    Makie.save(joinpath(output_dir_path, "i_v_pos.pdf"), fig)
    # ####

    if prms.study.non_isothermal
        temp_l_view = view(solution[var[dom_ids[:el_neg]].temp,:], subgrids.subgrid_neg) * T0_value
        temp_r_view = view(solution[var[dom_ids[:el_pos]].temp,:], subgrids.subgrid_pos) * T0_value

        temp_i_sol = solution[var[dom_ids[:sep]].temp,:] * T0_value
        temp_i_sol = temp_i_sol[temp_i_sol .!== NaN]

        axis_kwargs_cell[:title] = L"\text{Temperature}"

        ϵ = 1e-8
        coord_sep_l = coord_sep .- [prms.geom.lx_sep/2 * (1 - ϵ), 0.0] * L0_value
        coord_sep_r = coord_sep .+ [prms.geom.lx_sep/2 * (1 - ϵ), 0.0] * L0_value

        fig = contour_plot_2d([vcat(temp_l_view,
                                    temp_i_sol, temp_i_sol,
                                    temp_r_view)],
                              [hcat(coord_neg,
                                    coord_sep_l, coord_sep_r,
                                    coord_pos)];
                              label=latexstring("T\\;[\\mathrm{$(T0_unit_str)}]"),
                              axis_kwargs=axis_kwargs_cell)
        Makie.save(joinpath(output_dir_path, "temp.pdf"), fig)
    end

end

function contour_plot_2d(solutions, coords;
                        label="",
                        levels=30,
                        axis_kwargs=axis_kwargs)
    fig = Figure()
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
    δ = 1e-3
    if isapprox(data_range[2]-data_range[1], 0; atol=δ)
        data_range = (data_range[1] - δ, data_range[2] + δ)
    end
    # choose the levels according to the overall data limits
    levels = range(data_range[1],stop=data_range[2],length=levels)

    for idx in eachindex(coords)
        contourf!(fig_axis,
                vec(coords[idx][1,:]), vec(coords[idx][2,:]), vec(solutions[idx]);
                levels=levels,
                colorrange=data_range,
                extendlow=:auto, extendhigh=:auto,
                colormap=:viridis)
    end
    ticks = collect(range(data_range..., length=5))
    ticks[1] = round(ticks[1], sigdigits=5, RoundUp)
    ticks[end] = round(ticks[end], sigdigits=5, RoundDown)
    ticks[2:(end-1)] .= round.(ticks[2:(end-1)], sigdigits=5, RoundNearest)

    Colorbar(fig[:, end+1],
            colorrange=data_range,
            label=label,
            labelrotation=0,
            ticks=ticks)
    return fig
end

function plot_all_fields_1d(combined_solutions_1d)
    CairoMakie.activate!()

    path_dir = joinpath(pwd(), "output/steady_state/through_plane_plots")
    mkpath(path_dir)

    fig_p = plot_pressure_1d(combined_solutions_1d)
    Makie.save(joinpath(path_dir, "p.pdf"), fig_p)

    fig_phi = plot_phi_1d(combined_solutions_1d)
    Makie.save(joinpath(path_dir, "phi.pdf"), fig_phi)

    fig_c = plot_c_1d(combined_solutions_1d)
    Makie.save(joinpath(path_dir, "c.pdf"), fig_c)

    if haskey(combined_solutions_1d, :temp)
        fig_temp = plot_temp_1d(combined_solutions_1d)
        Makie.save(joinpath(path_dir, "temp.pdf"), fig_temp)
    end
end

# function combined_solution_fields_1d(solutions_on_subgrids, prms)
#     dict_solutions = Dict{Symbol, Tuple{Vector{Float64},Vector{Float64}}}()

#     (coords_neg, sol_neg) = copy.(solutions_on_subgrids[:subgrid_neg])
#     (coords_pos, sol_pos) = copy.(solutions_on_subgrids[:subgrid_pos])
#     (coords_el_neg, sol_el_neg) = copy.(solutions_on_subgrids[:subgrid_el_neg])
#     (coords_el_pos, sol_el_pos) = copy.(solutions_on_subgrids[:subgrid_el_pos])
#     (coords_sep, sol_sep) = copy.(solutions_on_subgrids[:subgrid_sep])

#     # Account for membrane thickness
#     coords_sep .+= prms.geom.lx_sep / 2
#     coords_el_pos .+= prms.geom.lx_sep
#     coords_pos .+= prms.geom.lx_sep

#     # Coordinates
#     coords = vcat(coords_neg, coords_pos)
#     coords_el_neg_pos = vcat(coords_el_neg, coords_el_pos)
#     coords_neg_sep_pos = vcat(coords_neg, coords_sep, coords_pos)

#     # Pressure field
#     sol_p = vcat(sol_el_neg[Int(p_l),:], sol_el_pos[Int(p_r),:])
#     dict_solutions[:p] = (coords_el_neg_pos, sol_p)

#     # Electrostatic fields
#     sol_phi_s = vcat(sol_neg[Int(ϕₛ_l),:], sol_pos[Int(ϕₛ_r),:])
#     dict_solutions[:phi_s] = (coords, sol_phi_s)

#     sol_phi_l = vcat(sol_el_neg[Int(ϕₗ_l),:], sol_el_pos[Int(ϕₗ_r),:])
#     dict_solutions[:phi_l] = (coords_el_neg_pos, sol_phi_l)

#     # Species concentrations
#     sol_c_ox_neg = vcat(sol_el_neg[Int(c_ox_neg_l),:], sol_el_pos[Int(c_ox_neg_r),:])
#     dict_solutions[:c_ox_neg] = (coords_el_neg_pos, sol_c_ox_neg)

#     sol_c_red_neg = vcat(sol_el_neg[Int(c_red_neg_l),:], sol_el_pos[Int(c_red_neg_r),:])
#     dict_solutions[:c_red_neg] = (coords_el_neg_pos, sol_c_red_neg)

#     sol_c_ox_pos = vcat(sol_el_neg[Int(c_ox_pos_l),:], sol_el_pos[Int(c_ox_pos_r),:])
#     dict_solutions[:c_ox_pos] = (coords_el_neg_pos, sol_c_ox_pos)

#     sol_c_red_pos = vcat(sol_el_neg[Int(c_red_pos_l),:], sol_el_pos[Int(c_red_pos_r),:])
#     dict_solutions[:c_red_pos] = (coords_el_neg_pos, sol_c_red_pos)

#     if prms.study.non_isothermal
#         # Temperature
#         sol_temp = vcat(sol_neg[Int(temp_l),:], sol_sep[Int(temp_i),:], sol_pos[Int(temp_r),:])
#         dict_solutions[:temp] = (coords_neg_sep_pos, sol_temp)
#     end

#     return dict_solutions
# end

# function plot_pressure_1d(dict_combined_solutions)
#     fig = Makie.Figure()
#     fig_axis = Makie.Axis(fig[1,1];
#                     title="Pressure",
#                     xlabel=L"x",
#                     ylabel=L"P")

#     (coords, sol) = dict_combined_solutions[:p]

#     if(maximum(sol) - minimum(sol) < 1e-6)
#         ylims!(fig_axis, (minimum(sol) - 1e-6, minimum(sol) + 1e-6))
#     end

#     lines_p = lines!(fig_axis, coords, sol,
#                      label=L"P", color=:blue)
#     scatter_p = scatter!(fig_axis, coords, sol,
#                          label=L"P", color=:blue)

#     Legend(fig[1,2], [[lines_p, scatter_p]], [L"P"])
#     return fig
# end

# function plot_phi_1d(dict_combined_solutions)
#     fig = Makie.Figure()
#     fig_axis = Makie.Axis(fig[1,1];
#                     title="Electrostatic Potentials",
#                     xlabel=L"x",
#                     ylabel=L"\phi")

#     (coords_phi_s, sol_phi_s) = dict_combined_solutions[:phi_s]
#     (coords_phi_l, sol_phi_l) = dict_combined_solutions[:phi_l]

#     lines_ϕₛ = lines!(fig_axis, coords_phi_s, sol_phi_s,
#                      label=L"\phi_s", color=:green)
#     scatter_ϕₛ = scatter!(fig_axis, coords_phi_s, sol_phi_s,
#                           label=L"\phi_s", color=:green)

#     lines_ϕₗ = lines!(fig_axis, coords_phi_l, sol_phi_l;
#                      label=L"\phi_l", color=:blue)
#     scatter_ϕₗ = scatter!(fig_axis, coords_phi_l, sol_phi_l;
#                          label=L"\phi_l", color=:blue)

#     Legend(fig[1,2], [[lines_ϕₛ, scatter_ϕₛ], [lines_ϕₗ, scatter_ϕₗ]], [L"\phi_s", L"\phi_l"])
#     return fig
# end

# function plot_c_1d(dict_combined_solutions)
#     fig = Makie.Figure()
#     fig_axis = Makie.Axis(fig[1,1];
#                     title="Species Molar Concentration",
#                     xlabel=L"x",
#                     ylabel=L"\text{Concentration}")

#     (coords_c_ox_neg, sol_c_ox_neg) = dict_combined_solutions[:c_ox_neg]
#     (coords_c_red_neg, sol_c_red_neg) = dict_combined_solutions[:c_red_neg]
#     (coords_c_ox_pos, sol_c_ox_pos) = dict_combined_solutions[:c_ox_pos]
#     (coords_c_red_pos, sol_c_red_pos) = dict_combined_solutions[:c_red_pos]

#     lines_c_ox_neg = lines!(fig_axis, coords_c_ox_neg, sol_c_ox_neg,
#                             color=:blue)
#     scatter_c_ox_neg = scatter!(fig_axis, coords_c_ox_neg, sol_c_ox_neg,
#                             color=:blue)

#     lines_c_red_neg = lines!(fig_axis, coords_c_red_neg, sol_c_red_neg,
#                             color=:blue3, linestyle=:dot)
#     scatter_c_red_neg = scatter!(fig_axis, coords_c_red_neg, sol_c_red_neg,
#                             color=:blue3)

#     lines_c_ox_pos = lines!(fig_axis, coords_c_ox_pos, sol_c_ox_pos,
#                             color=:red)
#     scatter_c_ox_pos = scatter!(fig_axis, coords_c_ox_pos, sol_c_ox_pos,
#                             color=:red)

#     lines_c_red_pos = lines!(fig_axis, coords_c_red_pos, sol_c_red_pos,
#                             color=:red3, linestyle=:dot)
#     scatter_c_red_pos = scatter!(fig_axis, coords_c_red_pos, sol_c_red_pos,
#                             color=:red3)

#     Legend(fig[1,2], [[lines_c_ox_neg, scatter_c_ox_neg],
#                         [lines_c_red_neg, scatter_c_red_neg],
#                         [lines_c_ox_pos, scatter_c_ox_pos],
#                         [lines_c_red_pos, scatter_c_red_pos]],
#                         [L"c_{\mathrm{ox,n}}", L"c_{\mathrm{red,n}}",
#                         L"c_{\mathrm{ox,p}}", L"c_{\mathrm{red,p}}"])

#     return fig
# end


# function plot_temp_1d(dict_combined_solutions)
#     fig = Makie.Figure()
#     fig_axis = Makie.Axis(fig[1,1];
#                     title="Temperature",
#                     xlabel=L"x",
#                     ylabel=L"\text{Temperature}")


#     (coords_temp, sol_temp) = dict_combined_solutions[:temp]

#     if(maximum(sol_temp) - minimum(sol_temp) < 1e-6)
#         ylims!(fig_axis, (minimum(sol_temp) - 1e-6, minimum(sol_temp) + 1e-6))
#     end

#     lines_temp = lines!(fig_axis, coords_temp, sol_temp,
#                         label=L"T", color=:blue)
#     scatter_temp = scatter!(fig_axis, coords_temp, sol_temp,
#                             label=L"T", color=:blue)

#     Legend(fig[1,2], [[lines_temp, scatter_temp]], [L"T"])
#     return fig
# end


# function project_solutions(solution, sys, grid, subgrids, data, dim=1)
#     # TODO: Generalize this implementation to allow for projections of arbitrary subdomains
#     # Currently, only the projection on the x-axis (dim=1) is supported!

#     num_variables = length(instances(Variables))

#     coords_x = grid.components[ExtendableGrids.Coordinates][dim, :]
#     sp = sortperm(coords_x)

#     coords_x_sorted = coords_x[sp]
#     solution_sorted = solution[:, sp]
#     replace!(solution_sorted, NaN=>0.0)

#     coords_x_sorted_unique = unique(coords_x_sorted)
#     solution_x = zeros(size(solution_sorted)[1], length(coords_x_sorted_unique))
#     count_x = zeros(length(coords_x_sorted_unique))

#     x_last = -Inf
#     idx_solution_x = 0
#     for idx in 1:size(solution_sorted)[2]
#         x_current = coords_x_sorted[idx]
#         if x_current > x_last + eps(x_current)
#             x_last = x_current
#             idx_solution_x += 1
#         end
#         count_x[idx_solution_x] += 1
#         solution_x[:, idx_solution_x] .+= solution_sorted[:, idx]
#     end
#     @assert idx_solution_x == size(solution_x)[2]
#     @assert isapprox(sum(count_x), length(coords_x))
#     @assert isapprox(sum(count_x), size(solution_sorted)[2])

#     for idx_row in 1:size(solution_x)[1]
#         solution_x[idx_row, :] ./= count_x
#     end

#     solutions_on_subgrids = Dict{Symbol, Tuple{Vector{Float64}, Matrix{Float64}}}()

#     for key_subgrid in keys(subgrids)
#         subgrid = subgrids[key_subgrid]
#         coords_x_subgrid = vec(subgrid.components[ExtendableGrids.Coordinates][dim, :])
#         coords_x_subgrid_sorted_unique = unique(sort(coords_x_subgrid))
#         idx_coords_x_subgrid = 1
#         solution_x_subgrid = zeros(size(solution_sorted)[1], length(coords_x_subgrid_sorted_unique))
#         for idx_coords_x in eachindex(coords_x_sorted)
#             if idx_coords_x_subgrid <= length(coords_x_subgrid_sorted_unique) &&
#                 isapprox(coords_x[idx_coords_x],
#                 coords_x_subgrid_sorted_unique[idx_coords_x_subgrid])
#                 solution_x_subgrid[:, idx_coords_x_subgrid] = solution_x[:, idx_coords_x]
#                 idx_coords_x_subgrid += 1
#             end
#         end

#         solutions_on_subgrids[key_subgrid] = (coords_x_subgrid_sorted_unique,
#                                               solution_x_subgrid)
#     end

#     return (coords_x_sorted_unique, solution_x, solutions_on_subgrids)
# end
