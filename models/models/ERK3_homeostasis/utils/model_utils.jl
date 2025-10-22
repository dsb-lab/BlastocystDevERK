function test(v, dv, comvar)
    fgfr  = comvar
    #NANOG = v[1]
    #GATA6 = v[2]
    #FGF   = v[3]
    @. dv[1] = ((an*(1/(1 + v[2]^coefb)))*(ane*(1/(1+v[4]^coefe))) - v[1]) * l
    @. dv[2] = (ag*(1/(1 + v[1]^coefg)) - v[2]) * l
    @. dv[3] = (af*(1/(1 + v[2]^coefd)) - v[3]) * l
    @. dv[4] = (ae*(fgfr^coefo)/(fgfr^coefo + 1) -v[4]) *l
    return dv
end
function ERK_model_3(v, dv, comvar)
    fgfr  = comvar
    @. dv[1] = (an*(Kn1^coefb/(Kn1^coefb + v[2]^coefb)) * ane*(1/(1+v[4]^coefe)) - dn*v[1])
    @. dv[2] = (ag*(Kg^coefg/(Kg^coefg + v[1]^coefg)) - dg*v[2])
    @. dv[3] = (af*(Kf^coefd/(Kf^coefd + v[2]^coefd)) - df*v[3])
    ae_combined = v[5].+1.0
    @. dv[4] = ((ae*ae_combined)*(fgfr^coefo)/(fgfr^coefo + Ke^coefo) - de*v[4])
    # FGFR2
    @. dv[5] = (afr*(v[2]^coeffrg)/(kffrg^coeffrg+v[2]^coeffrg) - v[5]*dfgfr2)

    return dv
end