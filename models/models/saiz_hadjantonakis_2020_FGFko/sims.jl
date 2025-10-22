if pwd() !== "/home/pablo/Desktop/PhD/projects/ICMModels/models/saiz_hadjantonakis_2020"
    println("changing wd to /home/pablo/Desktop/PhD/projects/ICMModels/models/saiz_hadjantonakis_2020")
    cd()
    cd("/home/pablo/Desktop/PhD/projects/ICMModels/models/saiz_hadjantonakis_2020")
end
PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",size=35)
@pyimport matplotlib.patches as patch
include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("../../utils/plot_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")

@time results = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
                          , h=0.0005  #integration time-step
                          , mh=0.001   #data measuring time-step
                          , Fth=Inf
                          , comext=0.0
                          , comKO=false);


rr = findfirst.(isequal.(unique(results.totals)), [results.totals])
N=10

tots = results.totals[rr]
FDPmat  = zeros(N,length(tots))
FEPImat = zeros(N,length(tots))
FPREmat = zeros(N,length(tots))
for n=1:N
@time _results = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
                        , h=0.0005  #integration time-step
                        , mh=0.001   #data measuring time-step
                        , Fth=Inf
                        , comext=0.0
                        , comKO=false);
    rr = findfirst.(isequal.(unique(_results.totals)), [_results.totals])
    tots = _results.totals[rr]
    FDPmat[n,:]  .= _results.fDP[rr]./tots
    FEPImat[n,:] .= _results.fEPI[rr]./tots
    FPREmat[n,:] .= _results.fPRE[rr]./tots
end
FDP  = mean(FDPmat, dims=1)[:].*100
FEPI = mean(FEPImat, dims=1)[:].*100
FPRE = mean(FPREmat, dims=1)[:].*100
FDPstd = std(FDPmat, dims=1)[:].*100
FEPIstd = std(FEPImat, dims=1)[:].*100
FPREstd = std(FPREmat, dims=1)[:].*100

plt.close("all")
plt.clf()
fig, ax1 = plt.subplots(figsize=(22,15)) # Create a new blank figure
ax1.plot(tots, FDP, color="purple", label="DP", linewidth=5)
ax1.fill_between(tots, FDP.-FDPstd, FDP.+FDPstd, color="purple", alpha=0.3)
ax1.plot(tots, FEPI, color="red", label="EPI", linewidth=5)
ax1.fill_between(tots, FEPI.-FEPIstd, FEPI.+FEPIstd, color="red", alpha=0.3)
ax1.plot(tots, FPRE, color="blue", label="PRE", linewidth=5)
ax1.fill_between(tots, FPRE.-FPREstd, FPRE.+FPREstd, color="blue", alpha=0.3)

ax1.spines["right"].set_visible(false)
ax1.spines["top"].set_visible(false)
ax1.set_ylabel(join([L"\%", " of ICM"]))
ax1.set_xlabel("Cell number")
ax1.set_ylim(-1,101)
#PyPlot.legend(bbox_to_anchor=[1.1,1.1])
#plt.title(join([L" $F_m$", "=$FGFmedium"]))
p_patch = patch.Patch(color="purple", label="BP")
r_patch = patch.Patch(color="red", label="EPI")
b_patch = patch.Patch(color="blue", label="PrE")
plt.legend(handles=[p_patch, r_patch, b_patch], loc=1)
plt.tight_layout()
gcf()

plt.savefig("sims.png")