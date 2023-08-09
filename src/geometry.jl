@kwdef struct FlowCellGeometry{T}
    lx_cc_neg::T
    lx_el_neg::T
    lx_sep::T
    lx_el_pos::T
    lx_cc_pos::T
    ly_cell::T
end

@kwdef struct Mesh1D{T}
    hx_cc_neg::Vector{T}
    hx_el_neg::Vector{T}
    hx_el_pos::Vector{T}
    hx_cc_pos::Vector{T}
end
@kwdef struct Mesh2D{T}
    hx_cc_neg::Vector{T}
    hx_el_neg::Vector{T}
    hx_el_pos::Vector{T}
    hx_cc_pos::Vector{T}
    hy_cell::Vector{T}
end

function mesh_1d(geom::FlowCellGeometry{T};
    hx_cc_neg=[geom.lx_cc_neg/4, geom.lx_cc_neg/4],
    hx_cc_pos=[geom.lx_cc_pos/4, geom.lx_cc_pos/4],
    hx_el_neg=[geom.lx_el_neg/8, geom.lx_el_neg/4, geom.lx_el_neg/32],
    hx_el_pos=[geom.lx_el_pos/32, geom.lx_el_pos/4, geom.lx_el_pos/8]) where T<:AbstractFloat

    Mesh1D(hx_cc_neg=hx_cc_neg, hx_cc_pos=hx_cc_pos,
           hx_el_neg=hx_el_neg, hx_el_pos=hx_el_pos)
end

function mesh_2d(geom::FlowCellGeometry{T};
    hx_cc_neg=[geom.lx_cc_neg/4, geom.lx_cc_neg/4],
    hx_cc_pos=[geom.lx_cc_pos/4, geom.lx_cc_pos/4],
    hx_el_neg=[geom.lx_el_neg/8, geom.lx_el_neg/4, geom.lx_el_neg/32],
    hx_el_pos=[geom.lx_el_pos/32, geom.lx_el_pos/4, geom.lx_el_pos/8],
    hy_cell=[geom.ly_cell/128, geom.ly_cell/128, geom.ly_cell/128]) where T<:AbstractFloat

    Mesh2D(hx_cc_neg=hx_cc_neg, hx_cc_pos=hx_cc_pos,
           hx_el_neg=hx_el_neg, hx_el_pos=hx_el_pos, hy_cell=hy_cell)
end

function cell_thickness(geom::FlowCellGeometry{T}) where T<:AbstractFloat
    geom.lx_cc_neg + geom.lx_el_neg + geom.lx_sep + geom.lx_el_pos + geom.lx_cc_pos
end

# Constructor of FlowCellGeometry
function flow_cell_geometry(; lx_cc_neg::T, lx_el_neg::T, lx_sep::T,
                              lx_el_pos::T, lx_cc_pos::T, ly_cell::T) where {T}
    FlowCellGeometry{T}(lx_cc_neg=lx_cc_neg, lx_el_neg=lx_el_neg, lx_sep=lx_sep,
        lx_el_pos=lx_el_pos, lx_cc_pos=lx_cc_pos, ly_cell=ly_cell)
end
