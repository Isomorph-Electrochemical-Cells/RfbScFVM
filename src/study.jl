function update_physics_data(data, Δϕₛ)
    bc_old = data.boundary
    bc_new = BoundaryConditions{RealType}(species_neg=bc_old.species_neg,
                                          species_pos=bc_old.species_pos,
                                          temp_amb=bc_old.temp_amb,
                                          p_in_neg=bc_old.p_in_neg,
                                          p_in_pos=bc_old.p_in_pos,
                                          v_out_neg=bc_old.v_out_neg,
                                          v_out_pos=bc_old.v_out_pos,
                                          ϕₛ_neg=-Δϕₛ/2,
                                          ϕₛ_pos=Δϕₛ/2)

    return (geom=data.geom,
            el_neg=data.el_neg,
            el_pos=data.el_pos,
            cc_neg=data.cc_neg,
            cc_pos=data.cc_pos,
            sep=data.sep,
            boundary=bc_new,
            discr=data.discr,
            study=data.study,
            scales=data.scales,
            scaling_params=data.scaling_params,
            electrolyte=data.electrolyte,
            pressure_boundary_type = data.pressure_boundary_type)
end
