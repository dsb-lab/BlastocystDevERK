function ERK_model_2(v, dv, comvar)
    fgfr  = comvar
    @. dv[1] = (an*(Kn1^coefb/(Kn1^coefb + v[2]^coefb)) + ane*(1/(1+v[4]^coefe)) - dn*v[1])
    @. dv[2] = (ag*(Kg^coefg/(Kg^coefg + v[1]^coefg)) - dg*v[2])
    @. dv[3] = (af*(Kf^coefd/(Kf^coefd + v[2]^coefd)) - df*v[3])
    @. dv[4] = (ae*(fgfr^coefo)/(fgfr^coefo + Ke^coefo) - de*v[4])
    return dv
end
