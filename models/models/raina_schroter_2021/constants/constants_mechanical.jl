const b   = exp10(-6) # we use the negative sign directly here.
const F0  = exp10(-4)
const F0red = 1.0/1.5
const Nmax = 50
const N_start  =20
const Nstartsd = 0.000001
const mu = 2
const cfr = 1.2 #coupling factor range
const minit = exp10(-6)
const rinit = 5
const rdiv  = 1.0/(2.0^(1.0/3.0))
const tdiv = 30.0
const sdiv = 0.5

const nth = 0.7
const gth = 0.4
const nth2 = 0.95
const gth2 = 0.05

const fatevarmax = 2.0
const upperth = nth*fatevarmax
const lowerth = gth*fatevarmax