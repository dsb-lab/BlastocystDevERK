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

@time _times, FBP, FBPsd, FEPI, FEPIsd, FPRE, FPREsd, NCELLS = runs(agentsimRM, 10)

times = _times.*90.0./56.0
PyPlot.close("all")
PyPlot.clf()
fig = figure(figsize=(10,7)) # Create a new blank figure
errorbar(times, FBP.*100, yerr=FBPsd.*100, color="purple", label="BP")
errorbar(times, FEPI.*100, yerr=FEPIsd.*100, color="red", label="EPI")
errorbar(times, FPRE.*100, yerr=FPREsd.*100, color="blue", label="PRE")
errorbar(times, ones(length(FBP)).*50.0, color="black", linestyle="--", label="50 %")
ylabel("% of ICM")
xlabel("scaled time")
#xlim(200,250)
ylim(-1,101)
PyPlot.legend(loc=0)
PyPlot.tight_layout()
gcf()
#PyPlot.savefig("/home/pablo/Desktop/RM_fates_times.png")

PyPlot.clf()
fig = figure(figsize=(10,7)) # Create a new blank figure

plot(mean(NCELLS,dims=2)[:], FBP.*100, color="purple", label="DP", linewidth=3)
plt.fill_between(mean(NCELLS,dims=2)[:], FBP.*100.0.-FBPsd.*100, FBP.*100.0.+FBPsd.*100, alpha=0.5, color="purple")

plot(mean(NCELLS,dims=2)[:], FEPI.*100, color="red", label="EPI", linewidth=3)
plt.fill_between(mean(NCELLS,dims=2)[:], FEPI.*100.0.-FEPIsd.*100, FEPI.*100.0.+FEPIsd.*100, alpha=0.5, color="red")

plot(mean(NCELLS,dims=2)[:], FPRE.*100, color="blue", label="PrE", linewidth=3)
plt.fill_between(mean(NCELLS,dims=2)[:], FPRE.*100.0.-FPREsd.*100, FPRE.*100.0.+FPREsd.*100, alpha=0.5, color="blue")

plot(mean(NCELLS,dims=2)[:], ones(length(FBP)).*50.0, color="black", linestyle="--", label="50 %")

ylabel("% of ICM")
xlabel("ICM cell count")
ylim(-1,101)

PyPlot.legend(loc=0)
PyPlot.tight_layout()
gcf()
#PyPlot.savefig("/home/pablo/Desktop/RM_fates_cellN.png")
