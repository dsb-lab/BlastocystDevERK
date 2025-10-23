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