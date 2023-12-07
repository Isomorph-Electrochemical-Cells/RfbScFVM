### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 77f26554-1a7c-11ee-2353-5b8000f42f71
begin
	try
		using Revise
	catch
		@info "Not using Revise"
	end

    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	import CairoMakie
	using RfbScFVM
	using PlutoVista
	using PlutoUI
	import PlutoUI: combine
	using GridVisualize
	using Unitful
	using HypertextLiteral: @htl, @htl_str
	using LaTeXStrings
end


# ╔═╡ 8f1b5627-7630-4eaa-a955-c2fa47614c37
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 7e655b36-29e0-46c7-a8fc-b64291bb1c05
md"""# RFB Cell Performance Model"""

# ╔═╡ ba1699bf-ac99-47c8-9b2a-693278b205e1
begin
	plotter = CairoMakie
	nothing
end

# ╔═╡ e9537640-f686-44c3-b055-5019776b561c
md"# Select Parameter File"

# ╔═╡ 76569a6d-84bc-4573-b0ea-6cf1411b2490
@bind input_param_file FilePicker([MIME("application/json")])

# ╔═╡ 9f632fd1-bd9d-442c-b637-97a338ab1d77
@bind relaad_file Button("Reload File")

# ╔═╡ 46ee0afc-d471-4d6e-87c2-542afeacfd60
begin
	relaad_file
	path_input_folder =  joinpath("..", "input")
	default_param_file_name = joinpath(path_input_folder, "example_mv_temptma_polarization_isothermal_1d.json")
	file_name = input_param_file === nothing ? default_param_file_name : joinpath(path_input_folder, input_param_file["name"])
	base_params = read_config_file(file_name)
	base_data = preprocess_parameters(base_params)
	nothing
end

# ╔═╡ cb2dc619-a5fc-4f76-93da-e7636824297f
md"# Create and visualize Mesh"

# ╔═╡ 9214d9ee-3983-4272-8b15-0c7033ba8109
md"# Adapt Model Parameters"

# ╔═╡ 017bdcd1-369f-43a4-ab0c-92a29e26b370
begin
	@bind selected_non_isothermal Select([false => "Isothermal", true => "Non-Isothermal"], default=base_data.study.non_isothermal)
end

# ╔═╡ cb5f9dc9-8d2c-4ec1-8a3c-8a6b0dc358a5
begin
	@bind selected_spatial_discr Select(["fvm_1d" => "1D FV Discretization", "fvm_2d" => "2D FV Discretization"]; default=base_data.discr.spatial_discr)
end

# ╔═╡ e743db7d-f57c-4dff-afe5-c7042937df76
begin
	function param_inputs(params::Vector; title="", ranges, defaults::Vector)
		return confirm(PlutoUI.combine() do Child

			inputs = [
				(@htl("""<li> $(params[idx][2]): </br>
				 $(Child(params[idx][1], Slider(ranges[idx], show_value=true, default=defaults[idx])))
				</li>"""))

				for idx in eachindex(params)
			]

			@htl("""
			<h2> $title </h2>
			<ul>
			$(inputs)
			</ul>
			""")
		end)
	end

	param_velocity = ["velocity_y" => "Velocity (y-direction) [mm/s]"]
	param_species_neg = "c_" .* base_data.boundary.species_neg.row .* "_neg" .=> "Concentration " .* base_data.boundary.species_neg.row .* "_neg" .* " [mol/l]"
	param_species_pos = "c_" .* base_data.boundary.species_pos.row .* "_pos" .=> "Concentration " .* base_data.boundary.species_pos.row .* "_pos" .* " [mol/l]"

	params_selectable = Pair{String, String}[]

	append!(params_selectable, param_velocity)
	append!(params_selectable, param_species_neg)
	append!(params_selectable, param_species_pos)

	defaults=[base_data.boundary.v_out_neg]
	append!(defaults, collect(base_data.boundary.species_neg))
	append!(defaults, collect(base_data.boundary.species_pos))
	range_velocity = 1.0:0.1:20.0
	range_species_concentration = 0.0:0.1:5.0
	ranges = [range_velocity]
	append!(ranges, collect(Iterators.repeated(range_species_concentration, length(param_species_neg)+length(param_species_pos))))

	@bind slider_values param_inputs(params_selectable, title="Operating Conditions", ranges=ranges, defaults=defaults)
end

# ╔═╡ bf21ec29-8bed-47e0-8663-778423d7692c
md"# Run Polarization Simulation"

# ╔═╡ 5384a63b-56ec-4253-a1c2-fe6c77718cd0
	@bind run_simulation Button("Run Simulation")

# ╔═╡ 1de2a57e-cf12-4224-983f-3dedf33cff4a
begin
	slider_values

	params = copy(base_params)

	params["discretization_parameters"]["spatial_discretization"] = selected_spatial_discr

	dict_species_inlet_neg = params["boundary_conditions"]["species_inlet_neg"]
	for idx in eachindex(dict_species_inlet_neg)
		dict_species = dict_species_inlet_neg[idx]

		idx_slider_value = findfirst(x->x==Symbol("c_" * base_data.boundary.species_neg.row[idx] * "_neg"), keys(slider_values))

		dict_species["concentration"] = unitful_to_dict(slider_values[idx_slider_value], "mol/l")
	end

	dict_species_inlet_pos = params["boundary_conditions"]["species_inlet_pos"]
	for idx in eachindex(dict_species_inlet_neg)
		dict_species = dict_species_inlet_pos[idx]

		idx_slider_value = findfirst(x->x==Symbol("c_" * base_data.boundary.species_pos.row[idx] * "_pos"), keys(slider_values))

		dict_species["concentration"] = unitful_to_dict(slider_values[idx_slider_value], "mol/l")
	end

	params["boundary_conditions"]["velocity_outlet_neg"] = unitful_to_dict(slider_values.velocity_y, "mm/s")

	params["boundary_conditions"]["velocity_outlet_pos"] = unitful_to_dict(slider_values.velocity_y, "mm/s")

	params["study_parameters"]["non_isothermal"] = selected_non_isothermal

	nothing
end

# ╔═╡ 4bbd4f38-32be-446b-9308-4e0ed79d3c6f
data = preprocess_parameters(params); # pre-process updated parameters

# ╔═╡ 16727286-e669-4132-9577-b8e53a11e009
# ╠═╡ show_logs = false
begin
	if data.discr.spatial_discr == "fvm_1d"
    	(grid, subgrids) = create_grid_1d(data.geom, data.mesh, data.dom, data.bnd)
		grid_plot_size = (800,400)
	elseif data.discr.spatial_discr == "fvm_2d"
		(grid, subgrids) = create_grid_2d(data.geom, data.mesh, data.dom, data.bnd)
		grid_plot_size = (1600,1600)
	else
		@error "unknown spatial discretization"
	end

    vis=GridVisualize.GridVisualizer(Plotter=plotter, layout=(1,1); size=grid_plot_size)
    gridplot!(vis, grid; show=true)
end

# ╔═╡ 76292fd9-0322-4479-81f0-d70d283b605c
begin
	CairoMakie.activate!()
	run_simulation
    ocv = run_ocv(data)
    sys = system(grid, data)
    polarization_results = run_polarization(sys, grid, subgrids, data, ocv.Δϕₛ)
	nothing
end

# ╔═╡ c01fa748-b21c-4cc4-8ad4-25313ef742c9
dict_polarization = postprocess_results(polarization_results, data); # post-process polarization data

# ╔═╡ 64a24b44-172f-431d-a2a4-adbd824e903d
md"# Polarization Curve"

# ╔═╡ b1e5a7ec-cfca-4e3f-8db4-f59a39e8b71f
begin
        fig_polarization = plot_polarization_curve(dict_polarization["polarization"])
end

# ╔═╡ 35bb661b-fe42-46d2-85e7-d395f8e2238a
md"# Spatially Resolved Fields"

# ╔═╡ e5503a88-6f37-484b-9ccf-38a50d71e4df
begin
	dim_polarization_results = dimensional(polarization_results, data)
	vec_Δϕₛ = round.(typeof(1.0u"V"), dim_polarization_results.vec_Δϕₛ,digits=5)
	vec_i = round.(typeof(1.0u"mA/cm^2"), dim_polarization_results.vec_i,digits=5)
	polarization_indices = collect(1:length(vec_Δϕₛ))
	@bind selected_cell_voltage Select(polarization_indices .=> zip(vec_Δϕₛ, vec_i))
end

# ╔═╡ 83b3908a-fb7e-40ed-9db0-81bbdf4f0f24
begin
	solutions = polarization_results.vec_solution;
	if data.discr.spatial_discr == "fvm_1d"
		lst_dict_figures = [plot_all_fields_1d(solution, grid, subgrids, data) for solution in solutions]
	else
        lst_solution_proj = [project_solutions(solution, sys, grid, subgrids, data)[2] for solution in solutions]
        grid_proj_1d, subgrids_proj_1d = create_grid_1d(data.geom, RfbScFVM.Mesh1D(data.mesh), data.dom, data.bnd)
		lst_dict_figures = [plot_all_fields_1d(solution_proj, grid_proj_1d, subgrids_proj_1d, data) for solution_proj in lst_solution_proj]
	end
	num_figures = length(lst_dict_figures)
	nothing
end

# ╔═╡ 2f8d7899-2fe7-4a77-b86f-b4201f62353e
begin
	dict_figures = sort(lst_dict_figures[selected_cell_voltage])
	keys_figures = collect(keys(dict_figures))
	@bind selected_figure (Select(keys_figures, default="c"))
end

# ╔═╡ 92e9cace-97c0-4ad7-9a9e-0d67cc4a7e06
dict_figures[selected_figure]

# ╔═╡ c5a82d16-4632-4d72-9d95-1b53d130d6a0
begin
	if data.discr.spatial_discr == "fvm_2d"
			lst_dict_figures_2d = [plot_all_fields_2d(solution, grid, subgrids, data) for solution in solutions]

		begin
			dict_figures_2d = sort(lst_dict_figures_2d[selected_cell_voltage])
			keys_figures_2d = collect(keys(dict_figures_2d))
			@bind selected_figure_2d (Select(keys_figures_2d))
		end
	end
end

# ╔═╡ d1a80cff-627c-4f10-9695-017104e2d967
if data.discr.spatial_discr == "fvm_2d"
	dict_figures_2d[selected_figure_2d]
end

# ╔═╡ 0755ef5b-0154-4619-bdc9-ad62e3b92dda
md"""# Appendix: Load Packages"""

# ╔═╡ Cell order:
# ╟─8f1b5627-7630-4eaa-a955-c2fa47614c37
# ╟─7e655b36-29e0-46c7-a8fc-b64291bb1c05
# ╟─ba1699bf-ac99-47c8-9b2a-693278b205e1
# ╟─e9537640-f686-44c3-b055-5019776b561c
# ╟─76569a6d-84bc-4573-b0ea-6cf1411b2490
# ╟─9f632fd1-bd9d-442c-b637-97a338ab1d77
# ╟─46ee0afc-d471-4d6e-87c2-542afeacfd60
# ╟─cb2dc619-a5fc-4f76-93da-e7636824297f
# ╟─16727286-e669-4132-9577-b8e53a11e009
# ╟─9214d9ee-3983-4272-8b15-0c7033ba8109
# ╟─017bdcd1-369f-43a4-ab0c-92a29e26b370
# ╟─cb5f9dc9-8d2c-4ec1-8a3c-8a6b0dc358a5
# ╟─e743db7d-f57c-4dff-afe5-c7042937df76
# ╟─bf21ec29-8bed-47e0-8663-778423d7692c
# ╟─5384a63b-56ec-4253-a1c2-fe6c77718cd0
# ╟─76292fd9-0322-4479-81f0-d70d283b605c
# ╟─1de2a57e-cf12-4224-983f-3dedf33cff4a
# ╟─4bbd4f38-32be-446b-9308-4e0ed79d3c6f
# ╟─c01fa748-b21c-4cc4-8ad4-25313ef742c9
# ╟─64a24b44-172f-431d-a2a4-adbd824e903d
# ╟─b1e5a7ec-cfca-4e3f-8db4-f59a39e8b71f
# ╟─35bb661b-fe42-46d2-85e7-d395f8e2238a
# ╟─e5503a88-6f37-484b-9ccf-38a50d71e4df
# ╟─2f8d7899-2fe7-4a77-b86f-b4201f62353e
# ╟─83b3908a-fb7e-40ed-9db0-81bbdf4f0f24
# ╟─92e9cace-97c0-4ad7-9a9e-0d67cc4a7e06
# ╟─c5a82d16-4632-4d72-9d95-1b53d130d6a0
# ╟─d1a80cff-627c-4f10-9695-017104e2d967
# ╟─0755ef5b-0154-4619-bdc9-ad62e3b92dda
# ╟─77f26554-1a7c-11ee-2353-5b8000f42f71
