function ERK_model_1(v, dv, comvar)
    fgfr  = comvar
    gatafgf = zero(dv[2])
    nanoginh = zero(dv[2])
    @. dv[1] = (an) / (1 + (v[4]/(Ken))^coefm) - (dn)*v[1]
    @. dv[2] = (((ag)*v[4]^coefp)/((Keg)^coefp + v[4]^coefp)) - (dg)*v[2]
    @. dv[3]   = (((af)*v[1]^coefq)/((Knf)^coefq + v[1]^coefq)) - (df)*v[3]
    # @. dv[4]   = ((ae+aern[i])*v[2]*fgf) / (fgfn*(1 + (v[1]/(Kne+Knern[i]))^coefn)) - (de+dern[i])*v[4] + bE+bErn[i]
    @. gatafgf   = (v[2]*fgfr^coefo)/(((Kfe)^coefo + fgfr^coefo))
    @. nanoginh  = ((1 + (v[1]/(Kne))^coefn))
    @. dv[4]   = (ae) * gatafgf / nanoginh - (de)*v[4]
    # println(dv[4][1:5])
    return dv
end