function domain_symbols(non_isothermal=false)
    return [:cc_neg, :el_neg, :el_pos, :cc_pos, :sep]
end

function domain_id(domain_sym)
    if domain_sym == :cc_neg
        return 1
    elseif domain_sym == :el_neg
        return 2
    elseif domain_sym == :el_pos
        return 3
    elseif domain_sym == :cc_pos
        return 4
    elseif domain_sym == :sep
        return 5
    end
    @assert false "unknown domain symbol used"
end

function boundary_id(adjacent_domains::Tuple{Symbol,Symbol})
    if adjacent_domains == (:cc_neg, :el_neg)
        return 7
    elseif adjacent_domains == (:el_neg, :el_pos)
        return 5
    elseif adjacent_domains == (:el_pos, :cc_pos)
        return 8
    end
    @assert false "unknown adjacent_domains"
end

function boundary_id(boundary::Symbol) # FIXME: Optimize
    if boundary == :cc_neg_left
        return 6
    elseif boundary == :cc_neg_el_neg
        return 7
    elseif boundary == :sep
        return 5
    elseif boundary == :el_pos_cc_pos
        return 8
    elseif boundary == :cc_pos_right
        return 9
    elseif boundary == :el_neg_inflow
        return 10
    elseif boundary == :el_neg_outflow
        return 11
    elseif boundary == :el_pos_inflow
        return 12
    elseif boundary == :el_pos_outflow
        return 13
    end
end

variable_symbols() = [:p, :ϕₛ, :ϕₗ, :c, :temp]
variable_symbols(num_species) = vcat([:p, :ϕₛ, :ϕₗ], fill(:c, num_species), :temp)


function variable_symbol_to_indices(num_species)
    Dict{Symbol, Vector{Int64}}(:p => [1], :ϕₛ => [2], :ϕₗ => [3],
         :c => collect(4:4+num_species-1), :temp => [4+num_species])
end

function variable_domain_definitions(non_isothermal=true)
    d = Dict(:p => [:el_neg, :el_pos], :ϕₛ => [:cc_neg, :el_neg, :el_pos, :cc_pos],
        :ϕₗ => [:el_neg, :el_pos], :c => [:el_neg, :el_pos])
    d[:temp] = non_isothermal ? [:cc_neg, :el_neg, :sep, :el_pos, :cc_pos] : []
    return d
end

function variable_continuity()
    # continuity of field variables across domain borders
    # Dict((:cc_neg, :el_neg) => true, (:el_neg, :sep) => false,
    #      (:sep, :el_pos) => false, (:el_pos, :cc_pos) => true)
    Dict((:cc_neg, :el_neg) => true, (:el_neg, :el_pos) => false,
         (:el_pos, :cc_pos) => true, (:cc_pos, :sep) => false)
end

function domain_definitions_table(num_species, non_isothermal=true)
    domains = domain_symbols(non_isothermal)
    num_domains = length(domains)
    variables = variable_symbols(num_species)
    num_variables = length(variables)

    enabled_variables = AxisArray(zeros(Bool, num_variables, num_domains),
                                        rows=variables, cols=domains)

    sym2idx = variable_symbol_to_indices(num_species)
    var_dom_def = variable_domain_definitions(non_isothermal)

    for var in variables
        subdomains = var_dom_def[var]
        @view(enabled_variables[rows=sym2idx[var], cols=subdomains]) .= true
    end

    return enabled_variables
end

function subdomain_variable_indices(domain_def, non_isothermal)
    (num_variables, num_domains) = size(domain_def)
    subdomain_var_idx = AxisArray(zeros(Int64, num_variables, num_domains),
                                 rows=AxisArrays.axes(domain_def)[1],
                                 cols=AxisArrays.axes(domain_def)[2])

    # vector indicating whether the variables on subdomain at index i are discontinuous
    # from the variables on subdomain at index i-1
    var_discontinuity = zeros(Bool, num_domains)
    var_discontinuity[1] = false
    var_continuity = variable_continuity()
    domains = domain_symbols(non_isothermal)
    for idx_subdomain in 1:length(domains)-1
        var_discontinuity[idx_subdomain+1] = !var_continuity[(domains[idx_subdomain],
                                                              domains[idx_subdomain+1])]
    end

    offset = 1
    for idx_row in 1:num_variables
        var_indices = cumsum(domain_def[rows=idx_row] .* var_discontinuity) .+ offset
        offset = var_indices[end]+1
        var_indices .*= domain_def[rows=idx_row] .== true
        @view(subdomain_var_idx[rows=idx_row]) .= var_indices
    end
    return subdomain_var_idx
end

function subdomain_variable_id_to_domain_ids(subdomain_var_idx_table)
    dict_var_id_to_domain_ids = Dict()
    domain_ids = AxisArrays.axes(subdomain_var_idx_table)[2]
    (num_rows, num_cols) = size(subdomain_var_idx_table)
    for idx_row in 1:num_rows
        row = @view(subdomain_var_idx_table[rows=idx_row])
        unique_var_ids = filter(x->(x>0),unique(row))
        @show unique_var_ids
        for var_id in unique_var_ids
            ids = domain_ids[row .== var_id]
            if !isempty(ids)
                dict_var_id_to_domain_ids[var_id] = domain_id.(ids)
            end
        end
    end
    return dict_var_id_to_domain_ids
end
