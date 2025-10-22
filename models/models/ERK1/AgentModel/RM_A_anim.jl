using PyPlot
using PyCall
@pyimport matplotlib.patches as patch
if pwd() !== "/home/pablo/Desktop/PhD/projects/ICMModels/Nestor/AgentModel/"
    println("setting wd to /home/pablo/Desktop/PhD/projects/ICMModels/Nestor/code/AgentModel")
    cd()
    cd("/home/pablo/Desktop/PhD/projects/ICMModels/Nestor/AgentModel/")
end

include("constants/constants_mechanical.jl")
include("constants/constants_biochemical_RM.jl")
include("utils/general_utils.jl")
include("utils/agentmodel_utils.jl")

@time _times, _fBP, _fEPI, _fPRE, _nanog, _X, _Y, _Z, _R, _CFATES = agentsimRM(h=0.001,mh = 0.05, simtime=300.0, keep=true)

times_, fBP, fEPI, fPRE = process_results(_times, _fBP, _fEPI, _fPRE)
totals = fBP .+ fEPI .+ fPRE
totals = totals[1:length(times_)]
R = _R[:,1:length(times_)]
X = _X[:,1:length(times_)]
Y = _Y[:,1:length(times_)]
Z = _Z[:,1:length(times_)]
CFATES = _CFATES[:,1:length(times_)]

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
writedlm("results/fates.dat",CFATES',',')
writedlm("results/times.dat",times',',')
writedlm("results/X.dat",X',',')
writedlm("results/Y.dat",Y',',')
writedlm("results/Z.dat",Z',',')
writedlm("results/R.dat",R',',')
