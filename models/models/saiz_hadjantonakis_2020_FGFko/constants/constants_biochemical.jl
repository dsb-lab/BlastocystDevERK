const alpha = 10.0
const K     = 0.9

s  = 1.0

relative_sd_params = 1.0

const dn = 1.0
const Kn = 1.0
const an = alpha*dn*Kn

const ag = 10.0
const dg = 2.0
const Kg = 5.0

const df = dg
const af = (df / (K*Kn)) * sqrt((dg*Kg)/ag)

const coefn = 2.0
const coefm = 2.0

const states0   = [5.0, 8.5, 0.01]
const states0sd = [0.0, 0.0, 0.0]
const sigma_var=0.3
#const states_divsd = [0.3, 0.3*16/10, 0.1]

