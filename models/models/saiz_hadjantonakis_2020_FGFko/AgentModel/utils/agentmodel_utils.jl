using Statistics
using LinearAlgebra
using Distributions

function agentsimEM(; h=0.001  #integration time-step
           , mh=0.5   #data measuring time-step
           , simtime=300.0
           , keep=false)
    times = range(0.0, step=h, stop=simtime)
    steps = length(times)

    uni = Uniform()
    nor = Normal()
    # Biochemical circuit variables
    #nanog = zeros(Nmax)
    NANOG = zeros(Nmax)
    GATA6 = zeros(Nmax)
    FGF   = zeros(Nmax)
    # derivative of Biochemical variables
    #dnanog = zeros(Nmax)
    dNANOG = zeros(Nmax)
    dGATA6 = zeros(Nmax)
    dFGF   = zeros(Nmax)

    # intermideate Biochemical variables
    #nanogi = zeros(Nmax)
    NANOGi = zeros(Nmax)
    GATA6i = zeros(Nmax)
    FGFi   = zeros(Nmax)

    # derivative of intermideate Biochemical variables
    #dnanogi = zeros(Nmax)
    dNANOGi = zeros(Nmax)
    dGATA6i = zeros(Nmax)
    dFGFi   = zeros(Nmax)

    FGFrng= zeros(Nmax*Nmax) # FGF range of effect for each cell, should be 1 or 0
    fgf  = 0
    fgfn = 0

    cfate = zeros(Nmax) #set according to the biochem vars

    #max nanog value for the model
    nanogmax = alpha/(1+1/(2*K)^(2*coefm))
    NANOGmax = nanogmax*Kn

    # Set initial values
    NANOG[1] = NANOG0*(1.0+NANOG0sd*(2.0*rand()-1.0))
    GATA6[1] = GATA60*(1.0+GATA60sd*(2.0*rand()-1.0))
    FGF[1]   = FGF0*(1.0+FGF0sd*(2.0*rand()-1.0))
    if NANOG[1]>nth*NANOGmax
        cfate[1] = 1
    elseif (NANOG[1]<gth*NANOGmax)
        cfate[1] = 2
    else
        cfate[1] = 0
    end
    # Spatial variables
    x = zeros(Nmax)
    y = zeros(Nmax)
    z = zeros(Nmax)
    vx = zeros(Nmax)
    vy = zeros(Nmax)
    vz = zeros(Nmax)

    Fx = zeros(Nmax*Nmax)
    Fy = zeros(Nmax*Nmax)
    Fz = zeros(Nmax*Nmax)

    # derivative of Spatial variables
    dx = zeros(Nmax)
    dy = zeros(Nmax)
    dz = zeros(Nmax)
    dvx = zeros(Nmax)
    dvy = zeros(Nmax)
    dvz = zeros(Nmax)

    # intermidiate Spatial variables
    xi = zeros(Nmax)
    yi = zeros(Nmax)
    zi = zeros(Nmax)
    vxi = zeros(Nmax)
    vyi = zeros(Nmax)
    vzi = zeros(Nmax)

    Fxi = zeros(Nmax*Nmax)
    Fyi = zeros(Nmax*Nmax)
    Fzi = zeros(Nmax*Nmax)

    # derivative of intermideate Spatial variables
    dxi = zeros(Nmax)
    dyi = zeros(Nmax)
    dzi = zeros(Nmax)
    dvxi = zeros(Nmax)
    dvyi = zeros(Nmax)
    dvzi = zeros(Nmax)

    #Celular vaiables
    m   = zeros(Nmax)
    im  = zeros(Nmax)
    r   = zeros(Nmax)
    bm  = zeros(Nmax)
    F0m = zeros(Nmax)
    nextdiv = zeros(Nmax)
    ndiv  = zeros(Nmax)

    r[1] = rinit
    m[1] = minit
    ndiv[1] = 1
    rng = rand(nor)
    rnu1 = rand(uni)
    rnu2 = rand(uni)

    Ncells   = 1
    Ncurrent = Ncells
    Nstart   = round(Int, N_start + Nstartsd*rand(nor))

    ct   = h
    nextdiv[1] = (ndiv[1] - sdiv + rnu1*2.0*sdiv)*tdiv
    d    = 0.0
    rij  = 0.0
    id   = 0.0
    ij   = 0
    ji   = 0
    rijd = 0.0
    h2   = h*0.5

    # MEASURED VARIABLES
    minterval = round(Int,mh/h)
    stepsm = round(Int, steps/minterval)
    timesm = zeros(Float64, stepsm)

    #nanogm = zeros(Nmax, stepsm)
    fBP  = zeros(stepsm)
    fEPI = zeros(stepsm)
    fPRE = zeros(stepsm)

    if keep
        NANOGm = zeros(Nmax, stepsm)
        GATA6m = zeros(Nmax, stepsm)
        FGFm   = zeros(Nmax, stepsm)
        xm = zeros(Nmax,stepsm)
        ym = zeros(Nmax,stepsm)
        zm = zeros(Nmax,stepsm)
        R  = zeros(Nmax, stepsm)
        CFATES  = zeros(Nmax, stepsm)
    end
    stepm = 1
    step  = 1
    while Ncells<Nstart

        ## Integration of cell movement. Heuns method.

        # First step of the integration
        # Compute F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((x[i] - x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                ij   = (i-1)*Ncells + j
                Fpre = (rijd - 1)*(mu*rijd-1)*id
                if (d < mu*rij)
                    Fx[ij] = Fpre*(x[i]-x[j])
                    Fy[ij] = Fpre*(y[i]-y[j])
                    Fz[ij] = Fpre*(z[i]-z[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Ncells+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]
            end
        end

        # Compute differentials
        dx .= vx
        dy .= vy
        dz .= vz
        im .= -1.0./m
        bm .= b.*im
        dvx .= bm.*vx
        dvy .= bm.*vy
        dvz .= bm.*vz

        F0m .= F0 ./ m
        for i=1:Ncells
            for j=1:Ncells
                if i!=j
                    ij=(i-1)*Ncells+j;
                    dvx[i] += F0m[i]*Fx[ij]
                    dvy[i] += F0m[i]*Fy[ij]
                    dvz[i] += F0m[i]*Fz[ij]
                end
            end
        end

        # Compute intermediate states
        @. xi = x + h*dx
        @. yi = y + h*dy
        @. zi = z + h*dz
        @. vxi = vx + h*dvx
        @. vyi = vy + h*dvy
        @. vzi = vz + h*dvz

        # Compute intermediate F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((xi[i] - xi[j])^2 + (yi[i]-yi[j])^2 + (zi[i]-zi[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                ij   = (i-1)*Ncells + j
                Fpre = (rijd -1)*(mu*rijd-1)*id
                if (d < mu*rij)
                    Fx[ij] = Fpre*(xi[i]-xi[j])
                    Fy[ij] = Fpre*(yi[i]-yi[j])
                    Fz[ij] = Fpre*(zi[i]-zi[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Ncells+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]
            end
        end

        # Compute differentials of the intermidiates
        dxi .= vxi
        dyi .= vyi
        dzi .= vzi

        #@. bm = b*-1.0/m #mass has not changed
        dvxi .= bm.*vxi
        dvyi .= bm.*vyi
        dvzi .= bm.*vzi

        #F0m .= F0 ./ m #Mass has not changed
        for i=1:Ncells
            for j=1:Ncells
                if j!=i
                    ij=(i-1)*Ncells+j
                    dvxi[i] += F0m[i]*Fx[ij]
                    dvyi[i] += F0m[i]*Fy[ij]
                    dvzi[i] += F0m[i]*Fz[ij]
                end
            end
        end

        # Compute final states
        @. x = x + h2*(dx+dxi)
        @. y = y + h2*(dy+dyi)
        @. z = z + h2*(dz+dzi)
        @. vx = vx + h2*(dvx+dvxi)
        @. vy = vy + h2*(dvy+dvyi)
        @. vz = vz + h2*(dvz+dvzi)

        ## Check for cell divisions
        ct = h*step
        Ncurrent = Ncells
        for i=1:Ncurrent
            #check if its division time
            if ct >= nextdiv[i]
                Ncells +=1 # We have now one more cell
                ndiv[i]+=1 #current cell has divided one more time
                ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
                rnu1 = rand(uni)
                rnu2 = rand(uni)

                #nextdiv[Ncells] = nextdiv[i] + tdiv + (rnu2*2*sdiv - sdiv)
                nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv

                #nextdiv[i] =  nextdiv[i] + tdiv + (rnu1*2*sdiv - sdiv)
                nextdiv[i] = (ndiv[i] - sdiv + rnu1*2.0*sdiv)*tdiv

                # Now we compute the displacement of both daughter cells
                rnu1 = rand(uni)*2.0*pi # Uniform distributed random #bw 0 to 2pi
                rnu2 = rand(uni)*2.0*pi
                x[Ncells] = x[i] + r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[Ncells] = y[i] + r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[Ncells] = z[i] + r[i]*0.5*cos(rnu1)
                # displacement is opposite
                x[i] = x[i] - r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[i] = y[i] - r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[i] = z[i] - r[i]*0.5*cos(rnu1)
                #velocity is inherited
                vx[Ncells] = vx[i]
                vy[Ncells] = vy[i]
                vz[Ncells] = vz[i]

                r[i] = rdiv*r[i]
                r[Ncells] = r[i]
                m[i] = m[i]*0.5
                m[Ncells] = m[i]

                NANOG[Ncells] = NANOG[i]*(1.0+NANOG_divsd*rand())
                GATA6[Ncells] = GATA6[i]*(1.0+GATA6_divsd*rand())
                FGF[Ncells]   = FGF[i]*(1.0+FGF_divsd*rand())

                cfate[Ncells] = cfate[i]
                NANOG[i] = NANOG[i]*(1.0+NANOG_divsd*rand())
                GATA6[i] = GATA6[i]*(1.0+GATA6_divsd*rand())
                FGF[i]   = FGF[i]*(1.0+FGF_divsd*rand())
            end
        end

        #Check if its measuring time
        if step%minterval==0
            nbp  = 0
            nepi = 0
            npre = 0
            for i=1:Ncells
                if cfate[i]==0
                    nbp += 1
                elseif cfate[i]==1
                    nepi += 1
                elseif cfate[i]==2
                    npre += 1
                end
            end
            timesm[stepm] = ct
            fBP[stepm]  = nbp
            fEPI[stepm] = nepi
            fPRE[stepm] = npre
            if keep
                NANOGm[:,stepm] .= NANOG
                GATA6m[:,stepm] .= GATA6
                FGFm[:,stepm]   .= FGF
                xm[:,stepm] .= x
                ym[:,stepm] .= y
                zm[:,stepm] .= z
                R[:,stepm]  .= r
                CFATES[:,stepm] .= cfate
            end
            stepm+=1
        end
        step += 1
    end
    seguimos = true

    while seguimos

        ## Integration of cell movement. Heuns method.

        # First step of the integration

        # Compute F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((x[i] - x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                ij   = (i-1)*Ncells + j
                Fpre = (rijd - 1)*(mu*rijd-1)*id
                if d < mu*rij
                    Fx[ij] = Fpre*(x[i]-x[j])
                    Fy[ij] = Fpre*(y[i]-y[j])
                    Fz[ij] = Fpre*(z[i]-z[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Ncells+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]
                #Change which cells are in range of FGF signalling
                if (d<(cfr*rij))
                    FGFrng[ij] = 1
                else
                    FGFrng[ij] = 0
                end
                FGFrng[ji] = FGFrng[ij]
            end
        end

        # Compute differentials
        dx .= vx
        dy .= vy
        dz .= vz

        im .= -1.0./m
        bm .= b.*im

        dvx .= bm.*vx
        dvy .= bm.*vy
        dvz .= bm.*vz

        F0m .= F0 ./ m
        for i=1:Ncells
            fgf  = FGF[i]
            fgfn = 1
            for j=1:Ncells
                if j!=1
                    ij=(i-1)*Ncells+j;
                    if cfate[i] == cfate[j]
                        dvx[i] += F0m[i]*Fx[ij]
                        dvy[i] += F0m[i]*Fy[ij]
                        dvz[i] += F0m[i]*Fz[ij]
                    else
                        dvx[i] += F0red*F0m[i]*Fx[ij]
                        dvy[i] += F0red*F0m[i]*Fy[ij]
                        dvz[i] += F0red*F0m[i]*Fz[ij]
                    end
                    if FGFrng[ij]==1
                        fgf  += FGF[j]
                        #fgf  += nanog[j]
                        fgfn += 1
                    end
                end
            end

            #Integrate biochemical model
            dNANOG[i] = an / (1 + ((fgf*GATA6[i])/(Kg*fgfn))^coefm) - dn*NANOG[i]
            dGATA6[i] = ag*fgf / ((1 + (NANOG[i]/Kn)^coefn)*fgfn) - dg*GATA6[i]
            dFGF[i]   = af*NANOG[i] - df*FGF[i]

            NANOGi[i] = NANOG[i] + h*dNANOG[i]
            GATA6i[i] = GATA6[i] + h*dGATA6[i]
            FGFi[i]   = FGF[i] + h*dFGF[i]
        end

        # Compute intermediate states
        @. xi = x + h*dx
        @. yi = y + h*dy
        @. zi = z + h*dz
        @. vxi = vx + h*dvx
        @. vyi = vy + h*dvy
        @. vzi = vz + h*dvz

        # Compute intermediate F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((xi[i] - xi[j])^2 + (yi[i]-yi[j])^2 + (zi[i]-zi[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                ij   = (i-1)*Ncells + j
                Fpre = (rijd -1)*(mu*rijd-1)*id
                if d < mu*rij
                    Fx[ij] = Fpre*(xi[i]-xi[j])
                    Fy[ij] = Fpre*(yi[i]-yi[j])
                    Fz[ij] = Fpre*(zi[i]-zi[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Ncells+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]

                #Change which cells are in range of FGF signalling
                if (d<cfr*rij)
                    FGFrng[ij] = 1
                else
                    FGFrng[ij] = 0
                end
                FGFrng[ji] = FGFrng[ij]
            end
        end

        # Compute differentials of the intermidiates
        dxi .= vxi
        dyi .= vyi
        dzi .= vzi

        #@. bm = b/m #mass has not changed
        dvxi .= bm.*vxi
        dvyi .= bm.*vyi
        dvzi .= bm.*vzi

        #F0m .= F0 ./ m #Mass has not changed
        for i=1:Ncells
            #fgf  = nanogi[i]
            fgf  = FGFi[i]
            fgfn = 1
            for j=1:Ncells
                if j!=i
                    ij=(i-1)*Ncells+j
                    if cfate[i] == cfate[j]
                        dvx[i] += F0m[i]*Fx[ij]
                        dvy[i] += F0m[i]*Fy[ij]
                        dvz[i] += F0m[i]*Fz[ij]
                    else
                        dvx[i] += F0red*F0m[i]*Fx[ij]
                        dvy[i] += F0red*F0m[i]*Fy[ij]
                        dvz[i] += F0red*F0m[i]*Fz[ij]
                    end
                    if FGFrng[ij]==1
                        fgf  += FGFi[j]
                        fgfn += 1
                    end
                end
            end
            #Integrate biochemical model
            dNANOGi[i] = an / (1 + ((fgf*GATA6i[i])/(Kg*fgfn))^coefm) - dn*NANOGi[i]
            dGATA6i[i] = ag*fgf / ((1 + (NANOGi[i]/Kn)^coefn)*fgfn) - dg*GATA6i[i]
            dFGFi[i]   = af*NANOGi[i] - df*FGFi[i]

            #if (NANOG[i]>gth2*NANOGmax)&(NANOG[i]<nth2*NANOGmax)
                NANOG[i] = NANOG[i] + h2*(dNANOG[i] + dNANOGi[i])
                GATA6[i] = GATA6[i] + h2*(dGATA6[i] + dGATA6i[i])
                FGF[i]   = FGF[i]   + h2*(dFGF[i]   + dFGFi[i])
            #end

            if (NANOG[i]>nth*NANOGmax) #& (cfate[i]==0)
                cfate[i] = 1
            elseif (NANOG[i]<gth*NANOGmax) #& (cfate[i]==0)
                cfate[i] = 2
            end
        end

        # Compute final states
        @. x = x + h2*(dx+dxi)
        @. y = y + h2*(dy+dyi)
        @. z = z + h2*(dz+dzi)
        @. vx = vx + h2*(dvx+dvxi)
        @. vy = vy + h2*(dvy+dvyi)
        @. vz = vz + h2*(dvz+dvzi)

        ## Check for cell divisions
        ct = h*step
        Ncurrent = Ncells

        for i=1:Ncurrent
            #check if its division time
            if ct >= nextdiv[i]
                Ncells +=1 # We have now one more cell
                ndiv[i]+=1 #current cell has divided one more time
                ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
                rnu1 = rand(uni)
                rnu2 = rand(uni)

                #nextdiv[Ncells] = nextdiv[i] + tdiv + (rnu2*2*sdiv - sdiv)
                nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv

                #nextdiv[i] =  nextdiv[i] + tdiv + (rnu1*2*sdiv - sdiv)
                nextdiv[i] = (ndiv[i] - sdiv + rnu1*2.0*sdiv)*tdiv

                # Now we compute the displacement of both daughter cells
                rnu1 = rand(uni)*2.0*pi # Uniform distributed random #bw 0 to 2pi
                rnu2 = rand(uni)*2.0*pi
                x[Ncells] = x[i] + r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[Ncells] = y[i] + r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[Ncells] = z[i] + r[i]*0.5*cos(rnu1)
                # displacement is opposite
                x[i] = x[i] - r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[i] = y[i] - r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[i] = z[i] - r[i]*0.5*cos(rnu1)
                #velocity is inherited
                vx[Ncells] = vx[i]
                vy[Ncells] = vy[i]
                vz[Ncells] = vz[i]

                r[i] = rdiv*r[i]
                r[Ncells] = r[i]
                m[i] = m[i]*0.5
                m[Ncells] = m[i]

                NANOG[Ncells] = NANOG[i]*(1.0+NANOG_divsd*rand(uni))
                GATA6[Ncells] = GATA6[i]*(1.0+GATA6_divsd*rand(uni))
                FGF[Ncells]   = FGF[i]*(1.0+FGF_divsd*rand(uni))

                cfate[Ncells] = cfate[i]
                NANOG[i] = NANOG[i]*(1.0+NANOG_divsd*rand(uni))
                GATA6[i] = GATA6[i]*(1.0+GATA6_divsd*rand(uni))
                FGF[i]   = FGF[i]*(1.0+FGF_divsd*rand(uni))

                if (Ncells >= Nmax) | (step == steps)
                    seguimos=false
                    nbp  = 0
                    nepi = 0
                    npre = 0
                    for i=1:Ncells
                        if cfate[i]==0
                            nbp += 1
                        elseif cfate[i]==1
                            nepi += 1
                        elseif cfate[i]==2
                            npre += 1
                        end
                    end
                    timesm[stepm] = ct

                    fBP[stepm]  = nbp
                    fEPI[stepm] = nepi
                    fPRE[stepm] = npre
                    if keep
                        NANOGm[:,stepm] .= NANOG
                        GATA6m[:,stepm] .= GATA6
                        FGFm[:,stepm]   .= FGF
                        xm[:,stepm] .= x
                        ym[:,stepm] .= y
                        zm[:,stepm] .= z
                        R[:,stepm]  .= r
                        CFATES[:,stepm] .= cfate
                    end
                    break
                end
            end
        end

        #Check if its measuring time
        if step%minterval==0
            nbp  = 0
            nepi = 0
            npre = 0
            for i=1:Ncells
                if cfate[i]==0
                    nbp += 1
                elseif cfate[i]==1
                    nepi += 1
                elseif cfate[i]==2
                    npre += 1
                end
            end
            timesm[stepm] = ct
            fBP[stepm]  = nbp
            fEPI[stepm] = nepi
            fPRE[stepm] = npre
            if keep
                NANOGm[:,stepm] .= NANOG
                GATA6m[:,stepm] .= GATA6
                FGFm[:,stepm]   .= FGF
                xm[:,stepm] .= x
                ym[:,stepm] .= y
                zm[:,stepm] .= z
                R[:,stepm]  .= r
                CFATES[:,stepm] .= cfate
            end
            stepm+=1
        end
        step += 1
    end
    if keep
        return timesm, fBP, fEPI, fPRE, NANOGm, GATA6m, FGFm, xm, ym, zm, R, CFATES
    else
        return timesm, fBP, fEPI, fPRE
    end
end

function agentsimRM(; h=0.001  #integration time-step
           , mh=0.5   #data measuring time-step
           , simtime=300.0
           , keep=false)

    times = range(0.0, step=h, stop=simtime)
    steps = length(times)

    uni = Uniform()
    nor = Normal()
    # Biochemical circuit variables
    nanog = zeros(Nmax)

    # derivative of Biochemical variables
    dnanog = zeros(Nmax)

    # intermideate Biochemical variables
    nanogi = zeros(Nmax)

    # derivative of intermideate Biochemical variables
    dnanogi = zeros(Nmax)

    FGFrng= zeros(Nmax*Nmax) # FGF range of effect for each cell, should be 1 or 0
    fgf  = 0
    fgfn = 0

    cfate = zeros(Nmax) #set according to the biochem vars

    #max nanog value for the model
    nanogmax = alpha/(1+1/(2*K)^(2*coefm))

    # Set initial values
    nanog[1] = nanog0 + nanog0sd*(2.0*rand(uni)-1.0)
    if nanog[1]>nth*nanogmax
        cfate[1] = 1
    elseif nanog[1]<gth*nanogmax
        cfate[1] = 2
    else
        cfate[1] = 0
    end

    # Spatial variables
    x = zeros(Nmax)
    y = zeros(Nmax)
    z = zeros(Nmax)
    vx = zeros(Nmax)
    vy = zeros(Nmax)
    vz = zeros(Nmax)

    Fx = zeros(Nmax*Nmax)
    Fy = zeros(Nmax*Nmax)
    Fz = zeros(Nmax*Nmax)

    # derivative of Spatial variables
    dx = zeros(Nmax)
    dy = zeros(Nmax)
    dz = zeros(Nmax)
    dvx = zeros(Nmax)
    dvy = zeros(Nmax)
    dvz = zeros(Nmax)

    # intermidiate Spatial variables
    xi = zeros(Nmax)
    yi = zeros(Nmax)
    zi = zeros(Nmax)
    vxi = zeros(Nmax)
    vyi = zeros(Nmax)
    vzi = zeros(Nmax)

    Fxi = zeros(Nmax*Nmax)
    Fyi = zeros(Nmax*Nmax)
    Fzi = zeros(Nmax*Nmax)

    # derivative of intermideate Spatial variables
    dxi = zeros(Nmax)
    dyi = zeros(Nmax)
    dzi = zeros(Nmax)
    dvxi = zeros(Nmax)
    dvyi = zeros(Nmax)
    dvzi = zeros(Nmax)

    #Celular vaiables
    m   = zeros(Nmax)
    im  = zeros(Nmax)
    r   = zeros(Nmax)
    bm  = zeros(Nmax)
    F0m = zeros(Nmax)
    nextdiv = zeros(Nmax)
    ndiv  = zeros(Nmax)

    r[1] = rinit
    m[1] = minit
    ndiv[1] = 1
    rng = rand(nor)
    rnu1 = rand(uni)
    rnu2 = rand(uni)

    Ncells   = 1
    Ncurrent = Ncells
    Nstart   = round(Int, N_start + Nstartsd*rand(nor))

    ct   = h
    nextdiv[1] = (ndiv[1] - sdiv + rnu1*2.0*sdiv)*tdiv
    d    = 0.0
    rij  = 0.0
    id   = 0.0
    ij   = 0
    ji   = 0
    rijd = 0.0
    h2   = h*0.5

    # MEASURED VARIABLES
    minterval = round(Int,mh/h)
    stepsm = round(Int, steps/minterval)
    timesm = zeros(Float64, stepsm)

    if keep
        nanogm = zeros(Nmax, stepsm)
        xm = zeros(Nmax,stepsm)
        ym = zeros(Nmax,stepsm)
        zm = zeros(Nmax,stepsm)
        R  = zeros(Nmax, stepsm)
        CFATES  = zeros(Nmax, stepsm)
    end

    fBP  = zeros(stepsm)
    fEPI = zeros(stepsm)
    fPRE = zeros(stepsm)

    stepm = 1
    step  = 1
    while Ncells<Nstart

        ## Integration of cell movement. Heuns method.

        # First step of the integration
        # Compute F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((x[i] - x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                ij   = (i-1)*Nmax + j
                Fpre = (rijd - 1)*(mu*rijd-1)*id
                if (d < mu*rij)
                    Fx[ij] = Fpre*(x[i]-x[j])
                    Fy[ij] = Fpre*(y[i]-y[j])
                    Fz[ij] = Fpre*(z[i]-z[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Nmax+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]
            end
        end

        # Compute differentials
        dx .= vx
        dy .= vy
        dz .= vz
        im .= -1.0./m
        bm .= b.*im
        dvx .= bm.*vx
        dvy .= bm.*vy
        dvz .= bm.*vz

        F0m .= F0 ./ m
        for i=1:Ncells
            for j=1:Ncells
                if i!=j
                    ij=(i-1)*Nmax+j;
                    dvx[i] += F0m[i]*Fx[ij]
                    dvy[i] += F0m[i]*Fy[ij]
                    dvz[i] += F0m[i]*Fz[ij]
                end
            end
        end

        # Compute intermediate states
        @. xi = x + h*dx
        @. yi = y + h*dy
        @. zi = z + h*dz
        @. vxi = vx + h*dvx
        @. vyi = vy + h*dvy
        @. vzi = vz + h*dvz

        # Compute intermediate F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((xi[i] - xi[j])^2 + (yi[i]-yi[j])^2 + (zi[i]-zi[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                ij   = (i-1)*Nmax + j
                Fpre = (rijd -1)*(mu*rijd-1)*id
                if (d < mu*rij)
                    Fx[ij] = Fpre*(xi[i]-xi[j])
                    Fy[ij] = Fpre*(yi[i]-yi[j])
                    Fz[ij] = Fpre*(zi[i]-zi[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Nmax+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]
            end
        end

        # Compute differentials of the intermidiates
        dxi .= vxi
        dyi .= vyi
        dzi .= vzi

        #@. bm = b*-1.0/m #mass has not changed
        dvxi .= bm.*vxi
        dvyi .= bm.*vyi
        dvzi .= bm.*vzi

        #F0m .= F0 ./ m #Mass has not changed
        for i=1:Ncells
            for j=1:Ncells
                if j!=i
                    ij=(i-1)*Nmax+j
                    dvxi[i] += F0m[i]*Fx[ij]
                    dvyi[i] += F0m[i]*Fy[ij]
                    dvzi[i] += F0m[i]*Fz[ij]
                end
            end
        end

        # Compute final states
        @. x = x + h2*(dx+dxi)
        @. y = y + h2*(dy+dyi)
        @. z = z + h2*(dz+dzi)
        @. vx = vx + h2*(dvx+dvxi)
        @. vy = vy + h2*(dvy+dvyi)
        @. vz = vz + h2*(dvz+dvzi)

        ## Check for cell divisions
        ct = h*step
        Ncurrent = Ncells
        for i=1:Ncurrent
            #check if its division time
            if ct >= nextdiv[i]
                Ncells +=1 # We have now one more cell
                ndiv[i]+=1 #current cell has divided one more time
                ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
                rnu1 = rand(uni)
                rnu2 = rand(uni)

                nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv
                nextdiv[i] = (ndiv[i] - sdiv + rnu2*2.0*sdiv)*tdiv

                # Now we compute the displacement of both daughter cells
                rnu1 = rand(uni)*2.0*pi # Uniform distributed random #bw 0 to 2pi
                rnu2 = rand(uni)*2.0*pi
                x[Ncells] = x[i] + r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[Ncells] = y[i] + r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[Ncells] = z[i] + r[i]*0.5*cos(rnu1)
                # displacement is opposite
                x[i] = x[i] - r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[i] = y[i] - r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[i] = z[i] - r[i]*0.5*cos(rnu1)
                #velocity is inherited
                vx[Ncells] = vx[i]
                vy[Ncells] = vy[i]
                vz[Ncells] = vz[i]

                r[i] = rdiv*r[i]
                r[Ncells] = r[i]
                m[i] = m[i]*0.5
                m[Ncells] = m[i]

                nanog[Ncells] = nanog[i]*(1.0+nanog_divsd*rand(uni))

                cfate[Ncells] = cfate[i]
                nanog[i] = nanog[i]*(1.0+nanog_divsd*rand(uni))
            end
        end

        #Check if its measuring time
        if step%minterval==0
            nbp  = 0
            nepi = 0
            npre = 0
            for i=1:Ncells
                if cfate[i]==0
                    nbp += 1
                elseif cfate[i]==1
                    nepi += 1
                elseif cfate[i]==2
                    npre += 1
                end
            end
            timesm[stepm] = ct
            fBP[stepm]  = nbp
            fEPI[stepm] = nepi
            fPRE[stepm] = npre
            if keep
                nanogm[:,stepm] .= nanog
                xm[:,stepm] .= x
                ym[:,stepm] .= y
                zm[:,stepm] .= z
                R[:,stepm]  .= r
                CFATES[:,stepm] .= cfate
            end
            stepm+=1
        end
        step += 1
    end
    seguimos = true
    while seguimos

        ## Integration of cell movement. Heuns method.

        # First step of the integration

        # Compute F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((x[i] - x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                #if rijd > 2.0
                #    rijd = 2.0
                #    id = rij*0.5
                #end
                ij   = (i-1)*Nmax + j
                Fpre = (rijd - 1)*(mu*rijd-1)*id
                if d < mu*rij
                    Fx[ij] = Fpre*(x[i]-x[j])
                    Fy[ij] = Fpre*(y[i]-y[j])
                    Fz[ij] = Fpre*(z[i]-z[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Nmax+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]

                #Change which cells are in range of FGF signalling
                if (d<(cfr*rij))
                    FGFrng[ij] = 1
                else
                    FGFrng[ij] = 0
                end
                FGFrng[ji] = FGFrng[ij]
            end
        end

        # Compute differentials
        dx .= vx
        dy .= vy
        dz .= vz

        im .= -1.0./m
        bm .= b.*im

        dvx .= bm.*vx
        dvy .= bm.*vy
        dvz .= bm.*vz

        F0m .= F0 ./ m
        for i=1:Ncells
            fgf  = nanog[i]
            #fgf  = FGF[i]
            fgfn = 1
            for j=1:Ncells
                if j!=1
                    ij=(i-1)*Nmax+j;
                    if cfate[i] == cfate[j]
                        dvx[i] += F0m[i]*Fx[ij]
                        dvy[i] += F0m[i]*Fy[ij]
                        dvz[i] += F0m[i]*Fz[ij]
                    else
                        dvx[i] += F0red*F0m[i]*Fx[ij]
                        dvy[i] += F0red*F0m[i]*Fy[ij]
                        dvz[i] += F0red*F0m[i]*Fz[ij]
                    end
                    if FGFrng[ij]==1
                        fgf  += nanog[j]
                        fgfn += 1
                    end
                end
            end

            #Integrate biochemical model
            dnanog[i] = alpha * ((1.0+nanog[i]^coefn)^coefm) / ((1.0 + nanog[i]^coefn)^coefm + (fgf/(K*fgfn))^(2*coefm)) - nanog[i]
            nanogi[i] = nanog[i] + h*dnanog[i]
        end

        # Compute intermediate states
        @. xi = x + h*dx
        @. yi = y + h*dy
        @. zi = z + h*dz
        @. vxi = vx + h*dvx
        @. vyi = vy + h*dvy
        @. vzi = vz + h*dvz

        # Compute intermediate F
        for i=2:Ncells # Lower triangular excluding diagonal
            for j=1:(i-1)
                d    = sqrt((x[i] - x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                rijd = rij*id
                #if rijd > 2.0
                #    rijd = 2.0
                #    id = rij*0.5
                #end
                ij   = (i-1)*Nmax + j
                Fpre = (rijd - 1)*(mu*rijd-1)*id
                if d < mu*rij
                    Fx[ij] = Fpre*(xi[i]-xi[j])
                    Fy[ij] = Fpre*(yi[i]-yi[j])
                    Fz[ij] = Fpre*(zi[i]-zi[j])
                else
                    Fx[ij] = 0.0
                    Fy[ij] = 0.0
                    Fz[ij] = 0.0
                end
                ji = (j-1)*Nmax+i
                Fx[ji] = -Fx[ij]
                Fy[ji] = -Fy[ij]
                Fz[ji] = -Fz[ij]

                #Change which cells are in range of FGF signalling
                if (d<cfr*rij)
                    FGFrng[ij] = 1
                else
                    FGFrng[ij] = 0
                end
                FGFrng[ji] = FGFrng[ij]
            end
        end

        # Compute differentials of the intermidiates
        dxi .= vxi
        dyi .= vyi
        dzi .= vzi

        #@. bm = b/m #mass has not changed
        dvxi .= bm.*vxi
        dvyi .= bm.*vyi
        dvzi .= bm.*vzi

        #F0m .= F0 ./ m #Mass has not changed
        for i=1:Ncells
            fgf  = nanogi[i]
            #fgf  = FGFi[i]
            fgfn = 1
            for j=1:Ncells
                if j!=i
                    ij=(i-1)*Nmax+j
                    if cfate[i] == cfate[j]
                        dvx[i] += F0m[i]*Fx[ij]
                        dvy[i] += F0m[i]*Fy[ij]
                        dvz[i] += F0m[i]*Fz[ij]
                    else
                        dvx[i] += F0red*F0m[i]*Fx[ij]
                        dvy[i] += F0red*F0m[i]*Fy[ij]
                        dvz[i] += F0red*F0m[i]*Fz[ij]
                    end
                    if FGFrng[ij]==1
                        fgf  += nanog[j]
                        fgfn += 1
                    end
                end
            end
            #Integrate biochemical model

            dnanogi[i] = alpha*((1+nanogi[i]^coefn)^coefm) / ((1 + nanogi[i]^coefn)^coefm + (fgf/(K*fgfn))^(2*coefm)) - nanogi[i]

            if (nanog[i]>gth2*nanogmax) & (nanog[i]<nth2*nanogmax)
                nanog[i] = nanog[i] + h2*(dnanog[i]+dnanogi[i])
            end

            if (nanog[i]>nth*nanogmax) & (cfate[i]==0)
                cfate[i] = 1
            elseif (nanog[i]<gth*nanogmax) & (cfate[i]==0)
                cfate[i] = 2
            #else
            #    cfate[i] = 0
            end
        end

        # Compute final states
        @. x = x + h2*(dx+dxi)
        @. y = y + h2*(dy+dyi)
        @. z = z + h2*(dz+dzi)
        @. vx = vx + h2*(dvx+dvxi)
        @. vy = vy + h2*(dvy+dvyi)
        @. vz = vz + h2*(dvz+dvzi)

        ## Check for cell divisions
        ct = h*step
        Ncurrent = Ncells

        for i=1:Ncurrent
            #check if its division time
            if ct >= nextdiv[i]
                Ncells +=1 # We have now one more cell
                ndiv[i]+=1 #current cell has divided one more time
                ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
                rnu1 = rand(uni)
                rnu2 = rand(uni)

                nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv

                nextdiv[i] = (ndiv[i] - sdiv + rnu2*2.0*sdiv)*tdiv

                # Now we compute the displacement of both daughter cells
                rnu1 = rand(uni)*2.0*pi # Uniform distributed random #bw 0 to 2pi
                rnu2 = rand(uni)*2.0*pi
                x[Ncells] = x[i] + r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[Ncells] = y[i] + r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[Ncells] = z[i] + r[i]*0.5*cos(rnu1)
                # displacement is opposite
                x[i] = x[i] - r[i]*0.5*sin(rnu1)*cos(rnu2)
                y[i] = y[i] - r[i]*0.5*sin(rnu1)*sin(rnu2)
                z[i] = z[i] - r[i]*0.5*cos(rnu1)
                #velocity is inherited
                vx[Ncells] = vx[i]
                vy[Ncells] = vy[i]
                vz[Ncells] = vz[i]

                r[i] = rdiv*r[i]
                r[Ncells] = r[i]
                m[i] = m[i]*0.5
                m[Ncells] = m[i]

                nanog[Ncells] = nanog[i]*(1.0+nanog_divsd*rand(uni))
                cfate[Ncells] = cfate[i]
                nanog[i] = nanog[i]*(1.0+nanog_divsd*rand(uni))

                if (Ncells >= Nmax) | (step == steps)
                    seguimos=false
                    nbp  = 0
                    nepi = 0
                    npre = 0
                    for i=1:Ncells
                        if cfate[i]==0
                            nbp += 1
                        elseif cfate[i]==1
                            nepi += 1
                        elseif cfate[i]==2
                            npre += 1
                        end
                    end
                    timesm[stepm] = ct

                    fBP[stepm]  = nbp
                    fEPI[stepm] = nepi
                    fPRE[stepm] = npre
                    if keep
                        nanogm[:,stepm] .= nanog
                        xm[:,stepm] .= x
                        ym[:,stepm] .= y
                        zm[:,stepm] .= z
                        R[:,stepm]  .= r
                        CFATES[:,stepm] .= cfate
                    end
                    break
                end
            end
        end

        #Check if its measuring time
        if step%minterval==0
            nbp  = 0
            nepi = 0
            npre = 0
            for i=1:Ncells
                if cfate[i]==0
                    nbp += 1
                elseif cfate[i]==1
                    nepi += 1
                elseif cfate[i]==2
                    npre += 1
                end
            end
            timesm[stepm] = ct
            if keep
                nanogm[:,stepm] .= nanog
                xm[:,stepm] .= x
                ym[:,stepm] .= y
                zm[:,stepm] .= z
                R[:,stepm]  .= r
                CFATES[:,stepm] .= cfate
            end
            fBP[stepm]  = nbp
            fEPI[stepm] = nepi
            fPRE[stepm] = npre

            stepm+=1
        end
        step += 1
    end
    if keep
        return timesm, fBP, fEPI, fPRE, nanogm, xm, ym, zm, R, CFATES
    end
    return timesm, fBP, fEPI, fPRE
end

function process_results(times, fBP, fEPI, fPRE)
    idx = 0
    for i=2:length(times)
        if times[i] == 0.0
            idx=i-1
            break
        end
    end
    totals = fBP .+ fEPI .+ fPRE
    return times[1:idx], fBP[1:idx], fEPI[1:idx], fPRE[1:idx]
end

function runs(modelused::Function, n::Int)
    _times, _fBP, _fEPI, _fPRE = modelused(h=0.001,mh = 0.5)
    times = deepcopy(_times)
    Ncells = _fBP .+ _fEPI .+ _fPRE
    _FBP   = zeros(length(_times), n)
    _FEPI  = zeros(length(_times), n)
    _FPRE  = zeros(length(_times), n)
    NCELLS = zeros(length(_times), n)
    IDX = zeros(Int, n)
    idx = -Inf
    for i=1:n
        _times, _FBP[:,i], _FEPI[:,i], _FPRE[:,i] = modelused(h=0.001,mh = 0.5)
        for j in eachindex(_times)
            if _times[j] == 0.0
                IDX[i] = j
                if j > idx+1
                    idx = j-1
                    times .= _times
                end
                break
            end
        end
    end
    for i=1:n
        for j=IDX[i]:length(_times)
            _FBP[j,i]  = _FBP[j-1,i]
            _FEPI[j,i] = _FEPI[j-1,i]
            _FPRE[j,i] = _FPRE[j-1,i]
        end
        NCELLS[:,i] = _FBP[:,i] .+ _FEPI[:,i] .+ _FPRE[:,i]
    end
    _FBP ./= NCELLS
    _FEPI ./= NCELLS
    _FPRE ./= NCELLS

    FBP  = mean(_FBP, dims=2)
    FEPI = mean(_FEPI, dims=2)
    FPRE = mean(_FPRE, dims=2)
    FBPsd  = std(_FBP, dims=2)
    FEPIsd = std(_FEPI, dims=2)
    FPREsd = std(_FPRE, dims=2)

    return times[1:idx], FBP[1:idx], FBPsd[1:idx], FEPI[1:idx], FEPIsd[1:idx], FPRE[1:idx], FPREsd[1:idx], NCELLS[1:idx,:]
end

function select_idxs(totals, n)
    stp  = totals[end]/n
    idxs = Vector{Int64}()
    j=1
    for i=1:length(totals)
        if totals[i] > j*stp
            push!(idxs, round(Int, i))
            j+=1
        end
    end
    push!(idxs, length(totals))
    return idxs
end

function select_values(NANOG, GATA6, FGF, CFATES)
    cfates = Vector{Vector{Int64}}()
    cols = Vector{Vector{Any}}()
    nanog = Vector{Vector{Float64}}()
    gata6 = Vector{Vector{Float64}}()
    fgf = Vector{Vector{Float64}}()
    for j=1:length(NANOG[1,:])
        ng = Vector{Float64}()
        gt = Vector{Float64}()
        fg = Vector{Float64}()
        cf = Vector{Int64}()
        co = Vector{Any}()
        for i=1:length(NANOG[:,1])
            if NANOG[i,j] != 0
                push!(ng, NANOG[i,j])
                push!(gt, GATA6[i,j])
                push!(fg, FGF[i,j])
                push!(cf, CFATES[i,j])
                if CFATES[i,j] == 0
                    push!(co, "purple")
                elseif CFATES[i,j] == 1
                    push!(co, "red")
                elseif CFATES[i,j] ==2
                    push!(co, "blue")
                end
            else
                break
            end
        end
        push!(cfates, cf)
        push!(nanog, ng)
        push!(gata6, gt)
        push!(fgf, fg)
        push!(cols, co)
    end
    return nanog, gata6, fgf, cfates, cols
end

function mycumsum!(x, X::Matrix, idxs)
    x[1] += sum(X[:,1:idxs[1]])
    for i=1:length(idxs)-1
        x[i+1] += sum(X[:,idxs[i]+1:idxs[i+1]])
    end
end

function pdf_extraction!(x, X::AbstractVector, idxs)
    for i=1:length(idxs)
        for j=1:idxs[i]
            for k=1:length(X[j])
                push!(x[i],X[j][k])
            end
        end
    end
end

function compute_pdfs(nbins, nn)
    NANOG = Vector{Vector{Float64}}()
    GATA6 = Vector{Vector{Float64}}()
    FGF = Vector{Vector{Float64}}()
    for bin=1:nbins
        push!(NANOG, Vector{Float64}())
        push!(GATA6, Vector{Float64}())
        push!(FGF, Vector{Float64}())
    end

    _times, _fBP, _fEPI, _fPRE, _NANOG, _GATA6, _FGF, _X, _Y, _Z, R, CFATES = agentsimEM(h=0.001,mh = 0.5, keep=true,simtime=300.0)

    times_, fBP, fEPI, fPRE = process_results(_times, _fBP, _fEPI, _fPRE)
    totals = fBP .+ fEPI .+ fPRE
    totals = totals[1:length(times_)]
    NANOG_ = _NANOG[:,1:length(times_)]
    GATA6_ = _GATA6[:,1:length(times_)]
    FGF_   = _FGF[:,1:length(times_)]
    CFATES_= CFATES[:,1:length(times_)]
    nanog, gata6, fgf, cfates, cols = select_values(NANOG_, GATA6_, FGF_, CFATES_)

    idxs = select_idxs(totals, nbins)
    tots  = totals[idxs]
    pdf_extraction!(NANOG,nanog,idxs)
    pdf_extraction!(GATA6,gata6,idxs)
    pdf_extraction!(FGF,fgf,idxs)
    for i=1:nn-1
        println(i)
        _times, _fBP, _fEPI, _fPRE, _NANOG, _GATA6, _FGF, _X, _Y, _Z, R, CFATES = agentsimEM(h=0.001,mh = 0.5, keep=true,simtime=300.0)

        times_, fBP, fEPI, fPRE = process_results(_times, _fBP, _fEPI, _fPRE)
        totals = fBP .+ fEPI .+ fPRE
        totals = totals[1:length(times_)]
        NANOG_ = _NANOG[:,1:length(times_)]
        GATA6_ = _GATA6[:,1:length(times_)]
        FGF_   = _FGF[:,1:length(times_)]
        CFATES_= CFATES[:,1:length(times_)]
        nanog, gata6, fgf, cfates, cols = select_values(NANOG_, GATA6_, FGF_, CFATES_)

        idxs = select_idxs(totals, nbins)
        pdf_extraction!(NANOG,nanog,idxs)
        pdf_extraction!(GATA6,gata6,idxs)
        pdf_extraction!(FGF,fgf,idxs)
    end
    return NANOG, GATA6, FGF, tots
end

function select_nanog(nanog; nanogth=5.0)
    ii = 0
    ik = true
    jj = 0
    jk = true
    for i=1:length(nanog[:,1])
        if ik
            if nanog[i,end] > nanogth
                ii = i
                ik = false
            end
        end
        if jk
            if nanog[i,end] < nanogth
                jj = i
                jk = false
            end
        end
    end
    return ii, jj
end

function find_start(totals, Nstart)
    for i=1:length(totals)
        if totals[i] >= Nstart
            return i
        end
    end
end
