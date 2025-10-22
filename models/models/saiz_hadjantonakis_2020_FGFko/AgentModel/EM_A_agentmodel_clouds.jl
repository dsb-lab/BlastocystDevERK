using PyPlot
using PyCall
@pyimport matplotlib.patches as patch
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
_totals = fBP .+ fEPI .+ fPRE
_totals = _totals[1:length(times_)]

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

n = 4
idxs = select_idxs(totals, n)

NANOG = NANOG[:,idxs]
GATA6 = GATA6[:,idxs]
FGF   = FGF[:,idxs]
tots  = totals[idxs]
CFATES = CFATES[:,idxs]


nanog, gata6, fgf, cfates, cols = select_values(NANOG, GATA6, FGF, CFATES)

fig = figure("cloud plots",figsize=(20,10)) # Create a new blank figure
subplot(241)
scatter(nanog[1], gata6[1], color=cols[1])
ylabel("GATA6")
N = tots[1]
title("N = $N")

subplot(242)
scatter(nanog[2], gata6[2], color=cols[2])
N = tots[2]
title("N = $N")

subplot(243)
scatter(nanog[3], gata6[3], color=cols[3])
N = tots[3]
title("N = $N")

subplot(244)
scatter(nanog[4], gata6[4], color=cols[4])
N = tots[4]
title("N = $N")
p_patch = patch.Patch(color="purple", label="BP")
r_patch = patch.Patch(color="red", label="EPI")
b_patch = patch.Patch(color="blue", label="PrE")
PyPlot.legend(handles=[p_patch, r_patch, b_patch])

subplot(245)
scatter(nanog[1], fgf[1], color=cols[1])
xlabel("NANOG")
ylabel("FGF")

subplot(246)
scatter(nanog[2], fgf[2], color=cols[2])
xlabel("NANOG")

subplot(247)
scatter(nanog[3], fgf[3], color=cols[3])
xlabel("NANOG")

subplot(248)
scatter(nanog[4], fgf[4], color=cols[4])
xlabel("NANOG")
plt.tight_layout()
gcf()