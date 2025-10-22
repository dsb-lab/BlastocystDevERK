function raina_schroter_2021(v, dv, comvar)
    fgfr  = comvar
    #NANOG = v[1]
    #GATA6 = v[2]
    #FGF   = v[3]
    @. dv[1] = (an*(1/(1 + v[2]^coefb)) + anf*(1/(1+fgfr^coefe)) - v[1]) * l
    @. dv[2] = (ag*(1/(1 + v[1]^coefg)) - v[2]) * l
    @. dv[3] = (af*(1/(1 + v[2]^coefd)) - v[3]) * l
    return dv
end