using StaticArrays
using ComponentArrays

@kwdef struct VariableIndices{VECTOR_TYPE, MATRIX_TYPE}
    p::VECTOR_TYPE
    c::MATRIX_TYPE
end


function run_test()
    bnd2dom = Dict(:cc_neg_left => :cc_neg, :cc_pos_right => :cc_pos,
    :el_neg_outflow => :el_neg, :el_pos_outflow => :el_pos)


    dom = (el_neg=1,el_pos=2)

    p=@SVector [1,2]
    c=@SMatrix [1 2;
                3 4]

    pn = (el_neg=1,el_pos=2)
    cn = (el_neg=[1,2],el_pos=[3,4])

    #vi = VariableIndices(p=p,c=c)
    #vi = VariableIndices(p=pn,c=cn)
    vi = (el_neg=VariableIndices(p=1,c=[1,2]),
          el_pos=VariableIndices(p=2,c=[3,4]))

    alloc = @allocated begin
        vi[dom[bnd2dom[:el_neg_outflow]]].p * vi[2].c[2]
        # p, c = @unpack vi[...]
        # vi.p[dom_ids[:el_neg]] * vi.c[dom_ids[:el_neg]][1]
        #vi.p[dom_ids[:el_neg]]+vi.c[1,2]-vi.p[2]*vi.c[2,2]
    end
    alloc > 0 && @show alloc

end

run_test()


# data.idx[data.dom[:el_pos]].p

# data.idx.p[:el_pos]
# data.idx[:el_pos].p
