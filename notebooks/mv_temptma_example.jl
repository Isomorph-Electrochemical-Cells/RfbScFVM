### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 77f26554-1a7c-11ee-2353-5b8000f42f71
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using RfbScFVM
	using PlutoVista
	using VoronoiFVM
	using GLMakie, CairoMakie
end


# ╔═╡ e9537640-f686-44c3-b055-5019776b561c
md"# Define Cell Geometry and Mesh"

# ╔═╡ 46ee0afc-d471-4d6e-87c2-542afeacfd60
begin
	params = read_simulation_parameters("../input/mv_temptma_polarization_soc50.json")

    (grid, subgrids) = create_grid_2d(params.geom, params.mesh)
end

# ╔═╡ cb2dc619-a5fc-4f76-93da-e7636824297f
md"## Mesh Visualization"

# ╔═╡ 16727286-e669-4132-9577-b8e53a11e009
# ╠═╡ show_logs = false
gridplot(grid, Plotter=GLMakie, size=(1500,2000))

# ╔═╡ c0a84e79-66de-47aa-9aaf-61204a986c30
md"# Define Material Parameters"

# ╔═╡ 3b598016-a273-4ecb-99a7-3da6df5ef130
md"# Define Chemistry"

# ╔═╡ 94edcf27-7fa7-45af-a184-68f771beae93
md"# Define Boundary Conditions"

# ╔═╡ d7bfe9a3-0425-46c5-9f99-17998e7578fc
md"# Define Operating Conditions"

# ╔═╡ b3b23c03-7dc0-48e2-8742-2c07e381a8dc
md"# Define Study Case"

# ╔═╡ bf21ec29-8bed-47e0-8663-778423d7692c
md"# Run Simulation"

# ╔═╡ a2b92792-47ee-4e61-b31a-fb9aa53f017c
md"# Post-Processing"

# ╔═╡ 35bb661b-fe42-46d2-85e7-d395f8e2238a
md"## Spatially Resolved Fields"

# ╔═╡ 879dbe1a-926e-4079-98c7-af9d56c421d4
md"## Voltage Efficiency"

# ╔═╡ Cell order:
# ╠═77f26554-1a7c-11ee-2353-5b8000f42f71
# ╟─e9537640-f686-44c3-b055-5019776b561c
# ╠═46ee0afc-d471-4d6e-87c2-542afeacfd60
# ╟─cb2dc619-a5fc-4f76-93da-e7636824297f
# ╠═16727286-e669-4132-9577-b8e53a11e009
# ╟─c0a84e79-66de-47aa-9aaf-61204a986c30
# ╟─3b598016-a273-4ecb-99a7-3da6df5ef130
# ╟─94edcf27-7fa7-45af-a184-68f771beae93
# ╠═d7bfe9a3-0425-46c5-9f99-17998e7578fc
# ╠═b3b23c03-7dc0-48e2-8742-2c07e381a8dc
# ╠═bf21ec29-8bed-47e0-8663-778423d7692c
# ╠═a2b92792-47ee-4e61-b31a-fb9aa53f017c
# ╠═35bb661b-fe42-46d2-85e7-d395f8e2238a
# ╠═879dbe1a-926e-4079-98c7-af9d56c421d4
