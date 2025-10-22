function saiz_hadjantonakis_2020(v, dv, comvar)
    fgfr  = comvar
    #NANOG = v[1]
    #GATA6 = v[2]
    #FGF   = v[3]
    @. dv[1] = an/(1 + ((fgfr*v[2])/Kg)^coefm) - dn*v[1]
    @. dv[2] = (ag*fgfr)/(1 + (v[1]/Kn)^coefn) - dg*v[2]
    @. dv[3] = af*v[1] - df*v[3]
    return dv
end