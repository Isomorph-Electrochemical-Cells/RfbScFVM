function effective_diffusivity_factor(vel_long, vel_trans, region, data)
    PE0 = data.scaling_params.PE0
    if region == Int(domain_id_el_neg)
        total_dispersion_factor = data.el_neg.Deff_const
        total_dispersion_factor += data.el_neg.Deff_linear * abs(PE0 * vel_long)
        total_dispersion_factor += data.el_neg.Deff_quadratic * (PE0 * vel_long)^2
        return total_dispersion_factor
    elseif region == Int(domain_id_el_pos)
        total_dispersion_factor = data.el_pos.Deff_const
        total_dispersion_factor += data.el_pos.Deff_linear * abs(PE0 * vel_long)
        total_dispersion_factor += data.el_pos.Deff_quadratic * (PE0 * vel_long)^2
        return total_dispersion_factor
    end
    return 0.0
end
