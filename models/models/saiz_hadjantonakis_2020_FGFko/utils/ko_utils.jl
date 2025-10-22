function process_results(times, fBP, fEPI, fPRE)
    idx = 0
    for i=2:length(times)
        if times[i] == 0.0
            idx=i-1
            break
        end
    end
    if idx==0
        idx=length(times)
    end
    totals = fBP .+ fEPI .+ fPRE
    return times[1:idx], fBP[1:idx], fEPI[1:idx], fPRE[1:idx]
end

function find_start(totals, Nstart)
    for i=1:length(totals)
        if totals[i] >= Nstart
            return i
        end
    end
end
function NGF_nullclines(FGFmedium, Nanogmax, Gata6max)
    Nrange = collect(range(0,step=0.1,stop=Nanogmax))
    Grange = collect(range(0,step=0.1,stop=Gata6max))
    N_Nnull = an./(dn.*(1.0 .+ (FGFmedium.*Grange./Kg).^coefm))
    G_Gnull = ag.*FGFmedium./(dg.*(1.0 .+ (Nrange./Kn).^coefn))
    return Nrange, Grange, N_Nnull, G_Gnull
end

function NGF_vfield(FGFmedium, Nanogmax, Gata6max;xylen=20)
    x = range(0, length=xylen, stop=Gata6max)
    y = range(0, length=xylen, stop=Nanogmax)
    X = x' .* ones(length(y))
    Y = ones(length(x))' .* y
    dx = ag.*FGFmedium./(1.0 .+ (Y./Kn).^coefn) .- dg.*X
    dy = an./(1.0 .+ (FGFmedium.*X./Kg).^coefm) .- dn.*Y
    return X, Y, dx, dy
end

function NGF_postprocessing(_times, _fBP, _fEPI, _fPRE, _NANOG, _GATA6, _FGF, _X, _Y, _Z, R, _CFATES)
    times_, fBP, fEPI, fPRE = process_results(_times, _fBP, _fEPI, _fPRE)
    _totals = fBP .+ fEPI .+ fPRE
    println(_totals)
    _totals = _totals[1:length(times_)]

    println(_totals)
    idx_start = find_start(_totals, N_start)-100
    NANOG = _NANOG[:,idx_start:length(times_)]
    GATA6 = _GATA6[:,idx_start:length(times_)]
    FGF = _FGF[:,idx_start:length(times_)]
    CFATES = _CFATES[:,idx_start:length(times_)]
    X = _X[:,idx_start:length(times_)]
    Y = _Y[:,idx_start:length(times_)]
    Z = _Z[:,idx_start:length(times_)]
    times = times_[idx_start:end]
    totals = _totals[idx_start:length(times_)]
    return fBP, fEPI, fPRE, _totals, NANOG, GATA6, FGF, CFATES, X, Y, Z, times, totals
end

function NG_FGFmed(state, params)
    dNANOG = an./(1.0 .+ (params[1].*state[2]./Kg).^coefm) .- dn.*state[1]
    dGATA6 = ag.*params[1]./(1.0 .+ (state[1]./Kn).^coefn) .- dg.*state[2]
	return dNANOG, dGATA6
end

function fate_phasespace(FGFmedium, Nanogmax, Gata6max;resolution=100, thrds=false)
    params=[FGFmedium]
    h = 0.001
    times = collect(range(0,step=h,stop=500))
    N0s = range(0.1, stop=Nanogmax+0.2*Nanogmax, length=resolution)
    G0s = range(0.1, stop=Gata6max+0.2*Gata6max, length=resolution)
    total_iterations = length(N0s) * length(G0s)

    if thrds
        N = zeros(Threads.nthreads())
        G = zeros(Threads.nthreads())
        state=Vector{Vector{Float64}}()
        for _=1:Threads.nthreads()
            push!(state, [0.0,0.0])
        end
    else
        N = 0.0
        G = 0.0
        state=[0.0, 0.0]
    end

    init_pointsEPIN = zeros(total_iterations)
    init_pointsEPIG = zeros(total_iterations)
    init_pointsPrEN = zeros(total_iterations)
    init_pointsPrEG = zeros(total_iterations)
    
    end_pointsEPIN = zeros(total_iterations)
    end_pointsEPIG = zeros(total_iterations)
    end_pointsPrEN = zeros(total_iterations)
    end_pointsPrEG = zeros(total_iterations)

    Nthr = nth*Nanogmax
    Gthr = gth*Nanogmax
    current_it = [1]
    cit = 1
    if thrds
        Threads.@threads for N0 in N0s
            id = Threads.threadid()
            for G0 in G0s
                cit = current_it[1]
                overprint("sim = $cit out of $total_iterations", i=cit)
                current_it[1] +=1
                state[id].=[N0, G0]
                for t=2:length(times) 
                    state[id] .= Heuns(NG_FGFmed, state[id], params, h)
                end
                if state[id][1]>Nthr
                    init_pointsEPIN[cit]=N0
                    init_pointsEPIG[cit]=G0
                    end_pointsEPIN[cit]=state[id][1]
                    end_pointsEPIG[cit]=state[id][2]
                elseif state[id][1]<Gthr
                    init_pointsPrEN[cit]=N0
                    init_pointsPrEG[cit]=G0
                    end_pointsPrEN[cit]=state[id][1]
                    end_pointsPrEG[cit]=state[id][2]
                end
            end
        end
    else
        for N0 in N0s            
            for G0 in G0s
                cit = current_it[1]
                overprint("sim = $cit out of $total_iterations", i=cit)
                current_it[1] +=1
                state.=[N0, G0]
                for t=2:length(times) 
                    state .= Heuns(NG_FGFmed, state, params, h)
                end
                if state[1]>Nthr
                    init_pointsEPIN[cit]=N0
                    init_pointsEPIG[cit]=G0
                    end_pointsEPIN[cit]=state[1]
                    end_pointsEPIG[cit]=state[2]
                elseif state[1]<Gthr
                    init_pointsPrEN[cit]=N0
                    init_pointsPrEG[cit]=G0
                    end_pointsPrEN[cit]=state[1]
                    end_pointsPrEG[cit]=state[2]
                end
            end
        end
    end
    return init_pointsEPIN, init_pointsEPIG, init_pointsPrEN, init_pointsPrEG, end_pointsEPIN, end_pointsEPIG, end_pointsPrEN, end_pointsPrEG
end