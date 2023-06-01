# Bernoulli function
function fvc_b(z)
    if abs(z)<sqrt(eps(typeof(z)))
        return one(z) - z/2 + z*z/12
    end
    return z/(exp(z)-one(z))
end

function fvc_flux_Δx(ul, ur, v, d, Δx)
    p = v*Δx/d # Cell Peclet number

    # Evaluate Bernoulli function
    bm = fvc_b(-p)
    bp = fvc_b(p)

    return d * (bm * ul - bp * ur)
end

function fvc_flux_Δx(ul, ur, vh, d)
    p = vh/d # Cell Peclet number

    # Evaluate Bernoulli function
    bm = fvc_b(-p)
    bp = fvc_b(p)

    return d * (bm * ul - bp * ur)
end
