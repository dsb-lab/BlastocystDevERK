using Distributions
using ElasticArrays

function agentsimICM(model, varsnames, comvarname; h=0.001  #integration time-step
           , mh=0.5   #data measuring time-step
           , Fth=Inf
           , comext=0.0
           , comKO = false
           , test_cfr=nothing)

    ndivs  = ceil(Int64, log2(Nmax))
    _steps = ceil(Int64, (ndivs*tdiv)/h)
    times  = Vector{Float64}()
    sizehint!(times, _steps)

    uni = Uniform()
    nor = Normal()

    # Number of variables
    nvars  = length(varsnames)

    # Biochemical circuit variables
    states = Vector{Vector{Float64}}()
    comstate = zeros(Nmax)

    # derivative of Biochemical variables
    dstates = Vector{Vector{Float64}}()

    # intermideate Biochemical variables
    statesi = Vector{Vector{Float64}}()

    # derivative of intermideate Biochemical variables 
    dstatesi = Vector{Vector{Float64}}()

    if test_cfr===nothing
        used_cfr=cfr
    else
        used_cfr=test_cfr
    end
    print(used_cfr)
    for var=1:nvars
        push!(states, zeros(Nmax))
        push!(dstates, zeros(Nmax))
        push!(statesi, zeros(Nmax))
        push!(dstatesi, zeros(Nmax))
    end

    # Index of the communication variable (presumably FGF)
    comvaridx = findfirst(comvarname .== varsnames)
    FGFrng= zeros(Nmax*Nmax) # FGF range of effect for each cell, should be 1 or 0
    fgf  = 0
    fgfn = 0

    cfate = zeros(Nmax) #set according to the biochem vars

    for var in eachindex(states0)
        states[var][1] = states0[var]*(1.0+states0sd[var]*(2.0*rand()-1.0))
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
    Fxabs = zeros(Nmax*Nmax)
    Fyabs = zeros(Nmax*Nmax)
    Fzabs = zeros(Nmax*Nmax)
    Fxidx = BitVector(zeros(Nmax*Nmax))
    Fyidx = BitVector(zeros(Nmax*Nmax))
    Fzidx = BitVector(zeros(Nmax*Nmax))

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

    # derivative of intermideate Spatial variables
    dxi = zeros(Nmax)
    dyi = zeros(Nmax)
    dzi = zeros(Nmax)
    dvxi = zeros(Nmax)
    dvyi = zeros(Nmax)
    dvzi = zeros(Nmax)

    # Celular vaiables
    m   = zeros(Nmax)
    im  = zeros(Nmax)
    r   = zeros(Nmax)
    bm  = zeros(Nmax)
    F0m = zeros(Nmax)
    nextdiv = zeros(Nmax)
    ndiv  = zeros(Nmax)
    tdivs = Vector{Vector{Float64}}()
    tdivs_idx = Vector{Vector{Int64}}()
    for i=1:Nmax
        push!(tdivs, ElasticArray{Float64}(undef, 0))
        sizehint!(tdivs[i], ndivs)
    end
    for i=1:Nmax
        push!(tdivs_idx, ElasticArray{Int64}(undef, 0))
        sizehint!(tdivs_idx[i], ndivs)
    end
    r[1] = rinit
    m[1] = minit
    ndiv[1] = 1
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
    _stepsm   = round(Int, _steps/minterval)
    stepsm    = 0
    timesm    = ElasticArray{Float64}(undef, 0)
    sizehint!(timesm, _stepsm)
        
    fBP  = ElasticArray{Int64}(undef, 0)
    fEPI = ElasticArray{Int64}(undef, 0)
    fPRE = ElasticArray{Int64}(undef, 0)

    vars = []
    for var=1:nvars
        currentvar = ElasticArray{Float64}(undef, Nmax, 0)
        sizehint!(currentvar, (Nmax, _stepsm))
        push!(vars, currentvar)
    end
    comvar = ElasticArray{Float64}(undef, Nmax, 0)
    xm = ElasticArray{Float64}(undef, Nmax, 0)
    ym = ElasticArray{Float64}(undef, Nmax, 0)
    zm = ElasticArray{Float64}(undef, Nmax, 0)
    R  = ElasticArray{Float64}(undef, Nmax, 0)
    CFATES  = ElasticArray{Float64}(undef, Nmax, 0)
    sizehint!(comvar, (Nmax, _stepsm))
    sizehint!(xm, (Nmax, _stepsm))
    sizehint!(ym, (Nmax, _stepsm))
    sizehint!(zm, (Nmax, _stepsm))
    sizehint!(R, (Nmax, _stepsm))
    sizehint!(CFATES, (Nmax, _stepsm))

    stepm = 1
    step  = 1

    if Nstart > Nmax
        Nstart = Nmax
    end
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
        Fxabs .= abs.(Fx)
        Fyabs .= abs.(Fy)
        Fzabs .= abs.(Fz)
        for idx in eachindex(Fxidx)
            Fxidx[idx] = Fxabs[idx]>Fth
            Fyidx[idx] = Fxabs[idx]>Fth
            Fzidx[idx] = Fxabs[idx]>Fth
        end
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
        Fxabs .= abs.(Fx)
        Fyabs .= abs.(Fy)
        Fzabs .= abs.(Fz)
        for idx in eachindex(Fxidx)
            Fxidx[idx] = Fxabs[idx]>Fth
            Fyidx[idx] = Fxabs[idx]>Fth
            Fzidx[idx] = Fxabs[idx]>Fth
        end
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
            #cfate[i] = check_fate(states[1][i])
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
                push!(tdivs[i], ct)
                push!(tdivs_idx[i], stepm)
                Ncells +=1 # We have now one more cell
                ndiv[i]+=1 #current cell has divided one more time
                ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
                rnu1 = rand(uni)
                rnu2 = rand(uni)
                nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv
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
                m[i] = m[i]*0.5
                r[Ncells] = r[i]
                m[Ncells] = m[i]
                # for var=1:nvars
                #     states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     while states[var][Ncells] < 0.0
                #         states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     end
                #     states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     while states[var][i] < 0.0
                #         states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     end
                # end
                for var=1:nvars
                    noise = states[var][i]*sigma_var*(rand()-0.5)
                    var1 = states[var][i] + noise
                    var2 = states[var][i] - noise     
                    while var1<0.0 || var2<0.0
                        println()
                        println(var1)
                        println(var2)
                        noise = states[var][i]*sigma_var[var]*(rand()-0.5)
                        var1 = states[var][i] + noise
                        var2 = states[var][i] - noise
                        println(var1)
                        println(var2)
                    end
                    states[var][Ncells] = var1
                    states[var][i] = var2
                end
                cfate[Ncells] = cfate[i]
            end
        end
        step += 1
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
            stepsm += 1
            push!(timesm, ct)
            push!(fBP, nbp)
            push!(fEPI, nepi)
            push!(fPRE, npre)
            for var=1:nvars
                append!(vars[var], states[var])
            end
            append!(comvar, comstate)
            append!(xm, x)
            append!(ym, y)
            append!(zm, z)
            append!(R, r)
            append!(CFATES, cfate)
        end
    end
    seguimos = true
    if Ncells >= Nmax
        seguimos = false
    end
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
                if (d<(used_cfr*rij))
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
        Fxabs .= abs.(Fx)
        Fyabs .= abs.(Fy)
        Fzabs .= abs.(Fz)
        for idx in eachindex(Fxidx)
            Fxidx[idx] = Fxabs[idx]>Fth
            Fyidx[idx] = Fxabs[idx]>Fth
            Fzidx[idx] = Fxabs[idx]>Fth
        end
        for i=1:Ncells
            fgf  = states[comvaridx][i]
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
                        fgf  += states[comvaridx][j]
                        fgfn += 1
                    end
                end
            end
            comstate[i] = !comKO*fgf/fgfn + comext
        end
        #Integrate biochemical model
        dstates .= model(states, dstates, comstate)
        for var in eachindex(states)
            statesi[var][1:Ncells] .= states[var][1:Ncells] .+ h.*dstates[var][1:Ncells]
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
                if (d<(used_cfr*rij))
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

        Fxabs .= abs.(Fx)
        Fyabs .= abs.(Fy)
        Fzabs .= abs.(Fz)
        for idx in eachindex(Fxidx)
            Fxidx[idx] = Fxabs[idx]>Fth
            Fyidx[idx] = Fxabs[idx]>Fth
            Fzidx[idx] = Fxabs[idx]>Fth
        end
        #F0m .= F0 ./ m #Mass has not changed
        for i=1:Ncells
            fgf  = states[comvaridx][i]
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
                        fgf  += states[comvaridx][j]
                        fgfn += 1
                    end
                end
            end
            comstate[i] = !comKO*fgf/fgfn + comext
        end
        #Integrate biochemical model
        dstates .= model(statesi, dstatesi, comstate)
        for var in eachindex(states)
            states[var][1:Ncells]  .= states[var][1:Ncells] .+ h2.*(dstates[var][1:Ncells] .+ dstatesi[var][1:Ncells])
        end 
        for i=1:Ncells
            if length(states)==1
                cfate[i] = check_fate(states[1][i])
            else
                cfate[i] = check_fate2(states[1][i], states[2][i])
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
                push!(tdivs[i], ct)
                push!(tdivs_idx[i], stepm)
                Ncells +=1 # We have now one more cell
                ndiv[i]+=1 #current cell has divided one more time
                ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
                rnu1 = rand(uni)
                rnu2 = rand(uni)

                nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv

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
                m[i] = m[i]*0.5
                r[Ncells] = r[i]
                m[Ncells] = m[i]

                # for var=1:nvars
                #     states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     while states[var][Ncells] < 0.0
                #         states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     end
                #     states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     while states[var][i] < 0.0
                #         states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
                #     end
                # end
                for var=1:nvars
                    noise = states[var][i]*sigma_var*(rand()-0.5)
                    var1 = states[var][i] + noise
                    var2 = states[var][i] - noise     
                    while var1<0.0 || var2<0.0
                        noise = states[var][i]*sigma_var*(rand()-0.5)
                        var1 = states[var][i] + noise
                        var2 = states[var][i] - noise
                    end
                    states[var][Ncells] = var1
                    states[var][i] = var2
                end
                cfate[Ncells] = cfate[i]
                cfate[Ncells] = cfate[i]

                if (Ncells >= Nmax)
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
                    stepsm += 1
                    push!(timesm, ct)
                    push!(fBP, nbp)
                    push!(fEPI, nepi)
                    push!(fPRE, npre)
                    for var=1:nvars
                        append!(vars[var], states[var])
                    end
                    append!(comvar, comstate)
                    append!(xm, x)
                    append!(ym, y)
                    append!(zm, z)
                    append!(R, r)
                    append!(CFATES, cfate)
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
            stepsm += 1
            push!(timesm, ct)
            push!(fBP, nbp)
            push!(fEPI, nepi)
            push!(fPRE, npre)
            for var=1:nvars
                append!(vars[var], states[var])
            end
            append!(comvar, comstate)
            append!(xm, x)
            append!(ym, y)
            append!(zm, z)
            append!(R, r)
            append!(CFATES, cfate)
        end
        step += 1
    end
    totals = fBP .+ fEPI .+ fPRE
    return SimResults(timesm, fBP, fEPI, fPRE, totals, vars, comvar, varsnames,  comvarname,Positions(xm, ym, zm), R, CFATES, tdivs, tdivs_idx)
end

function check_fate(var)
    if var>upperth #& (cfate[i]==0)
        return 1
    elseif var<lowerth #& (cfate[i]==0)
        return 2
    else
        return 0
    end
end

function check_fate2(var1, var2)
    if var1>upperth && var2<upperth 
        return 1
    elseif var2>upperth && var1<upperth 
        return 2
    else
        return 0
    end
end


function agentsimICM_ESCs(model, varsnames, comvarname; h=0.001  #integration time-step
    , mh=0.5   #data measuring time-step
    , Fth=Inf
    , comext=0.0
    , comKO = false
    , test_cfr=nothing
    , NESCs = 0.0
    , fixedFGF=0.0)

    NmaxESC = Nmax + NESCs
    ndivs  = ceil(Int64, log2(NmaxESC))
    _steps = ceil(Int64, (ndivs*tdiv)/h)
    times  = Vector{Float64}()
    sizehint!(times, _steps)

    uni = Uniform()
    nor = Normal()

    # Number of variables
    nvars  = length(varsnames)

    # Biochemical circuit variables
    states = Vector{Vector{Float64}}()
    comstate = zeros(NmaxESC)

    # derivative of Biochemical variables
    dstates = Vector{Vector{Float64}}()

    # intermideate Biochemical variables
    statesi = Vector{Vector{Float64}}()

    # derivative of intermideate Biochemical variables 
    dstatesi = Vector{Vector{Float64}}()

    if test_cfr===nothing
    used_cfr=cfr
    else
    used_cfr=test_cfr
    end
    print(used_cfr)
    for var=1:nvars
    push!(states, zeros(NmaxESC))
    push!(dstates, zeros(NmaxESC))
    push!(statesi, zeros(NmaxESC))
    push!(dstatesi, zeros(NmaxESC))
    end

    # Index of the communication variable (presumably FGF)
    comvaridx = findfirst(comvarname .== varsnames)
    FGFrng= zeros(NmaxESC*NmaxESC) # FGF range of effect for each cell, should be 1 or 0
    fgf  = 0
    fgfn = 0

    cfate = zeros(NmaxESC) #set according to the biochem vars

    for var in eachindex(states0)
    states[var][1] = states0[var]*(1.0+states0sd[var]*(2.0*rand()-1.0))
    end

    # Spatial variables
    x = zeros(NmaxESC)
    y = zeros(NmaxESC)
    z = zeros(NmaxESC)
    vx = zeros(NmaxESC)
    vy = zeros(NmaxESC)
    vz = zeros(NmaxESC)

    Fx = zeros(NmaxESC*NmaxESC)
    Fy = zeros(NmaxESC*NmaxESC)
    Fz = zeros(NmaxESC*NmaxESC)
    Fxabs = zeros(NmaxESC*NmaxESC)
    Fyabs = zeros(NmaxESC*NmaxESC)
    Fzabs = zeros(NmaxESC*NmaxESC)
    Fxidx = BitVector(zeros(NmaxESC*NmaxESC))
    Fyidx = BitVector(zeros(NmaxESC*NmaxESC))
    Fzidx = BitVector(zeros(NmaxESC*NmaxESC))

    # derivative of Spatial variables
    dx = zeros(NmaxESC)
    dy = zeros(NmaxESC)
    dz = zeros(NmaxESC)
    dvx = zeros(NmaxESC)
    dvy = zeros(NmaxESC)
    dvz = zeros(NmaxESC)

    # intermidiate Spatial variables
    xi = zeros(NmaxESC)
    yi = zeros(NmaxESC)
    zi = zeros(NmaxESC)
    vxi = zeros(NmaxESC)
    vyi = zeros(NmaxESC)
    vzi = zeros(NmaxESC)

    # derivative of intermideate Spatial variables
    dxi = zeros(NmaxESC)
    dyi = zeros(NmaxESC)
    dzi = zeros(NmaxESC)
    dvxi = zeros(NmaxESC)
    dvyi = zeros(NmaxESC)
    dvzi = zeros(NmaxESC)

    # Celular vaiables
    m   = zeros(NmaxESC)
    im  = zeros(NmaxESC)
    r   = zeros(NmaxESC)
    bm  = zeros(NmaxESC)
    F0m = zeros(NmaxESC)
    nextdiv = zeros(NmaxESC)
    ndiv  = zeros(NmaxESC)
    tdivs = Vector{Vector{Float64}}()
    tdivs_idx = Vector{Vector{Int64}}()
    for i=1:NmaxESC
    push!(tdivs, ElasticArray{Float64}(undef, 0))
    sizehint!(tdivs[i], ndivs)
    end
    for i=1:NmaxESC
    push!(tdivs_idx, ElasticArray{Int64}(undef, 0))
    sizehint!(tdivs_idx[i], ndivs)
    end
    r[1] = rinit
    m[1] = minit
    ndiv[1] = 1
    rnu1 = rand(uni)
    rnu2 = rand(uni)

    Ncells   = 1
    Ncurrent = Ncells
    Nstart   = round(Int, (N_start+NESCs) + Nstartsd*rand(nor))
    ESC_cell_ids = Vector{Int64}(range(Nstart-NESCs-1, Nstart))
    
    println()
    println("NMAX is $NmaxESC")
    println("Nstart is $Nstart")

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
    _stepsm   = round(Int, _steps/minterval)
    stepsm    = 0
    timesm    = ElasticArray{Float64}(undef, 0)
    sizehint!(timesm, _stepsm)
    
    fBP  = ElasticArray{Int64}(undef, 0)
    fEPI = ElasticArray{Int64}(undef, 0)
    fPRE = ElasticArray{Int64}(undef, 0)

    vars = []
    for var=1:nvars
    currentvar = ElasticArray{Float64}(undef, NmaxESC, 0)
    sizehint!(currentvar, (NmaxESC, _stepsm))
    push!(vars, currentvar)
    end
    comvar = ElasticArray{Float64}(undef, NmaxESC, 0)
    xm = ElasticArray{Float64}(undef, NmaxESC, 0)
    ym = ElasticArray{Float64}(undef, NmaxESC, 0)
    zm = ElasticArray{Float64}(undef, NmaxESC, 0)
    R  = ElasticArray{Float64}(undef, NmaxESC, 0)
    CFATES  = ElasticArray{Float64}(undef, NmaxESC, 0)
    sizehint!(comvar, (NmaxESC, _stepsm))
    sizehint!(xm, (NmaxESC, _stepsm))
    sizehint!(ym, (NmaxESC, _stepsm))
    sizehint!(zm, (NmaxESC, _stepsm))
    sizehint!(R, (NmaxESC, _stepsm))
    sizehint!(CFATES, (NmaxESC, _stepsm))

    stepm = 1
    step  = 1

    if Nstart > NmaxESC
    Nstart = NmaxESC
    end
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
            ij   = (i-1)*NmaxESC + j
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
            ji = (j-1)*NmaxESC+i
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
    Fxabs .= abs.(Fx)
    Fyabs .= abs.(Fy)
    Fzabs .= abs.(Fz)
    for idx in eachindex(Fxidx)
        Fxidx[idx] = Fxabs[idx]>Fth
        Fyidx[idx] = Fxabs[idx]>Fth
        Fzidx[idx] = Fxabs[idx]>Fth
    end
    for i=1:Ncells
        for j=1:Ncells
            if i!=j
                ij=(i-1)*NmaxESC+j;
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
            ij   = (i-1)*NmaxESC + j
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
            ji = (j-1)*NmaxESC+i
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
    Fxabs .= abs.(Fx)
    Fyabs .= abs.(Fy)
    Fzabs .= abs.(Fz)
    for idx in eachindex(Fxidx)
        Fxidx[idx] = Fxabs[idx]>Fth
        Fyidx[idx] = Fxabs[idx]>Fth
        Fzidx[idx] = Fxabs[idx]>Fth
    end
    #F0m .= F0 ./ m #Mass has not changed
    for i=1:Ncells
        for j=1:Ncells
            if j!=i
                ij=(i-1)*NmaxESC+j
                dvxi[i] += F0m[i]*Fx[ij]
                dvyi[i] += F0m[i]*Fy[ij]
                dvzi[i] += F0m[i]*Fz[ij]
            end
        end
        #cfate[i] = check_fate(states[1][i])
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
            push!(tdivs[i], ct)
            push!(tdivs_idx[i], stepm)
            Ncells +=1 # We have now one more cell
            ndiv[i]+=1 #current cell has divided one more time
            ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
            rnu1 = rand(uni)
            rnu2 = rand(uni)
            nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv
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
            m[i] = m[i]*0.5
            r[Ncells] = r[i]
            m[Ncells] = m[i]
            # for var=1:nvars
            #     states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     while states[var][Ncells] < 0.0
            #         states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     end
            #     states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     while states[var][i] < 0.0
            #         states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     end
            # end
            for var=1:nvars
                noise = states[var][i]*sigma_var*(rand()-0.5)
                var1 = states[var][i] + noise
                var2 = states[var][i] - noise     
                while var1<0.0 || var2<0.0
                    println()
                    println(var1)
                    println(var2)
                    noise = states[var][i]*sigma_var[var]*(rand()-0.5)
                    var1 = states[var][i] + noise
                    var2 = states[var][i] - noise
                    println(var1)
                    println(var2)
                end
                states[var][Ncells] = var1
                states[var][i] = var2
            end
            cfate[Ncells] = cfate[i]
            if i in ESC_cell_ids
                cfate[Ncells] = 1
                push!(ESC_cell_ids, Ncells)
            end
        end
    end
    step += 1
    if step%minterval==0
        nbp  = 0
        nepi = 0
        npre = 0
        for i=1:Ncells
            if i in ESC_cell_ids
                continue
            end
            if cfate[i]==0
                nbp += 1
            elseif cfate[i]==1
                nepi += 1
            elseif cfate[i]==2
                npre += 1
            end
        end
        stepsm += 1
        push!(timesm, ct)
        push!(fBP, nbp)
        push!(fEPI, nepi)
        push!(fPRE, npre)
        for var=1:nvars
            append!(vars[var], states[var])
        end
        append!(comvar, comstate)
        append!(xm, x)
        append!(ym, y)
        append!(zm, z)
        append!(R, r)
        append!(CFATES, cfate)
    end
    end
    seguimos = true
    if Ncells >= NmaxESC
    seguimos = false
    end
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
            ij   = (i-1)*NmaxESC + j
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
            ji = (j-1)*NmaxESC+i
            Fx[ji] = -Fx[ij]
            Fy[ji] = -Fy[ij]
            Fz[ji] = -Fz[ij]

            #Change which cells are in range of FGF signalling
            if (d<(used_cfr*rij))
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
    Fxabs .= abs.(Fx)
    Fyabs .= abs.(Fy)
    Fzabs .= abs.(Fz)
    for idx in eachindex(Fxidx)
        Fxidx[idx] = Fxabs[idx]>Fth
        Fyidx[idx] = Fxabs[idx]>Fth
        Fzidx[idx] = Fxabs[idx]>Fth
    end
    for i=1:Ncells
        fgf  = states[comvaridx][i]
        fgfn = 1
        for j=1:Ncells
            if j!=1
                ij=(i-1)*NmaxESC+j;
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
                    if i <= NESCs
                        fgf += fixedFGF
                    else
                        fgf += states[comvaridx][j]
                    end
                    fgfn += 1
                end
            end
        end
        comstate[i] = !comKO*fgf/fgfn + comext
    end
    #Integrate biochemical model
    dstates .= model(states, dstates, comstate)
    for var in eachindex(states)
        statesi[var][1:Ncells] .= states[var][1:Ncells] .+ h.*dstates[var][1:Ncells]
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
            ij   = (i-1)*NmaxESC + j
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
            ji = (j-1)*NmaxESC+i
            Fx[ji] = -Fx[ij]
            Fy[ji] = -Fy[ij]
            Fz[ji] = -Fz[ij]

            #Change which cells are in range of FGF signalling
            if (d<(used_cfr*rij))
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

    Fxabs .= abs.(Fx)
    Fyabs .= abs.(Fy)
    Fzabs .= abs.(Fz)
    for idx in eachindex(Fxidx)
        Fxidx[idx] = Fxabs[idx]>Fth
        Fyidx[idx] = Fxabs[idx]>Fth
        Fzidx[idx] = Fxabs[idx]>Fth
    end
    #F0m .= F0 ./ m #Mass has not changed
    for i=1:Ncells
        fgf  = states[comvaridx][i]
        fgfn = 1
        for j=1:Ncells
            if j!=1
                ij=(i-1)*NmaxESC+j;
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
                    if i <= NESCs
                        fgf += fixedFGF
                    else
                        fgf += states[comvaridx][j]
                    end
                    fgfn += 1
                end
            end
        end
        comstate[i] = !comKO*fgf/fgfn + comext
    end
    #Integrate biochemical model
    dstates .= model(statesi, dstatesi, comstate)
    for var in eachindex(states)
        states[var][1:Ncells]  .= states[var][1:Ncells] .+ h2.*(dstates[var][1:Ncells] .+ dstatesi[var][1:Ncells])
    end 
    for i=1:Ncells
        if length(states)==1
            cfate[i] = check_fate(states[1][i])
        else
            cfate[i] = check_fate2(states[1][i], states[2][i])
        end
        if i in ESC_cell_ids
            cfate[i] = 1
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
            push!(tdivs[i], ct)
            push!(tdivs_idx[i], stepm)
            Ncells +=1 # We have now one more cell
            ndiv[i]+=1 #current cell has divided one more time
            ndiv[Ncells] = ndiv[i]; #daughter cell is equal as mother
            rnu1 = rand(uni)
            rnu2 = rand(uni)

            nextdiv[Ncells] = (ndiv[Ncells] - sdiv + rnu1*2.0*sdiv)*tdiv

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
            m[i] = m[i]*0.5
            r[Ncells] = r[i]
            m[Ncells] = m[i]

            # for var=1:nvars
            #     states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     while states[var][Ncells] < 0.0
            #         states[var][Ncells] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     end
            #     states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     while states[var][i] < 0.0
            #         states[var][i] = states[var][i] + states_divsd[var]*(rand()-0.5)
            #     end
            # end
            for var=1:nvars
                noise = states[var][i]*sigma_var*(rand()-0.5)
                var1 = states[var][i] + noise
                var2 = states[var][i] - noise     
                while var1<0.0 || var2<0.0
                    noise = states[var][i]*sigma_var*(rand()-0.5)
                    var1 = states[var][i] + noise
                    var2 = states[var][i] - noise
                end
                states[var][Ncells] = var1
                states[var][i] = var2
            end
            cfate[Ncells] = cfate[i]
            if i in ESC_cell_ids
                push!(ESC_cell_ids, Ncells)
            end
            if (Ncells >= NmaxESC)
                seguimos=false
                nbp  = 0
                nepi = 0
                npre = 0
                for i=1:Ncells
                    if i in ESC_cell_ids
                        continue
                    end
                    if cfate[i]==0
                        nbp += 1
                    elseif cfate[i]==1
                        nepi += 1
                    elseif cfate[i]==2
                        npre += 1
                    end
                end
                stepsm += 1
                push!(timesm, ct)
                push!(fBP, nbp)
                push!(fEPI, nepi)
                push!(fPRE, npre)
                for var=1:nvars
                    append!(vars[var], states[var])
                end
                append!(comvar, comstate)
                append!(xm, x)
                append!(ym, y)
                append!(zm, z)
                append!(R, r)
                append!(CFATES, cfate)
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
            if i in ESC_cell_ids
                continue
            end
            if cfate[i]==0
                nbp += 1
            elseif cfate[i]==1
                nepi += 1
            elseif cfate[i]==2
                npre += 1
            end
        end
        stepsm += 1
        push!(timesm, ct)
        push!(fBP, nbp)
        push!(fEPI, nepi)
        push!(fPRE, npre)
        for var=1:nvars
            append!(vars[var], states[var])
        end
        append!(comvar, comstate)
        append!(xm, x)
        append!(ym, y)
        append!(zm, z)
        append!(R, r)
        append!(CFATES, cfate)
    end
    step += 1
    end
    totals = fBP .+ fEPI .+ fPRE
    return SimResults(timesm, fBP, fEPI, fPRE, totals, vars, comvar, varsnames,  comvarname,Positions(xm, ym, zm), R, CFATES, tdivs, tdivs_idx), ESC_cell_ids
end
