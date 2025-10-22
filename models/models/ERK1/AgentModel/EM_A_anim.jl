if pwd() !== "/home/pablo/Desktop/PhD/projects/EmbryonicDev/models/AgentModel"
    println("changing wd to home/pablo/Desktop/PhD/projects/EmbryonicDev/models/AgentModel")
    cd()
    cd("/home/pablo/Desktop/PhD/projects/EmbryonicDev/models/AgentModel")
end

include("constants/constants_mechanical.jl")
include("constants/constants_biochemical_EM.jl")
include("utils/general_utils.jl")
include("utils/agentmodel_utils.jl")

@time _times, _fBP, _fEPI, _fPRE, _NANOG, _GATA6, _FGF, _X, _Y, _Z, R, _CFATES = agentsimEM(h=0.0001,mh = 0.01, keep=true,simtime=300.0)

times_, fBP, fEPI, fPRE = process_results(_times, _fBP, _fEPI, _fPRE)
totals = fBP .+ fEPI .+ fPRE
totals = totals[1:length(times_)]
R = R[:,1:length(times_)]
X = _X[:,1:length(times_)]
Y = _Y[:,1:length(times_)]
Z = _Z[:,1:length(times_)]

function nonan(A)
    for i=1:size(A)[1]
        for j=1:size(A)[2]
            if isnan(A[i,j])
                A[i,j] = 100.0
            end
        end
    end
    return A
end

X .= nonan(X)
Y .= nonan(Y)
Z .= nonan(Z)
times = times_ .*90.0./56.0

using DelimitedFiles
writedlm("results/times.dat",times',',')
writedlm("results/X.dat",X',',')
writedlm("results/Y.dat",Y',',')
writedlm("results/Z.dat",Z',',')
writedlm("results/R.dat",R',',')
