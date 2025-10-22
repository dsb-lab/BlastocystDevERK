using Plots
pyplot()
if pwd() !== "/home/pablo/Desktop/PhD/projects/EmbryonicDev/models/AgentModel"
    println("changing wd to home/pablo/Desktop/PhD/projects/EmbryonicDev/models/AgentModel")
    cd()
    cd("/home/pablo/Desktop/PhD/projects/EmbryonicDev/models/AgentModel")
end

include("constants/constants_mechanical.jl")
include("constants/constants_biochemical_EM.jl")
include("utils/general_utils.jl")
include("utils/agentmodel_utils.jl")

@time _times, FBP, FBPsd, FEPI, FEPIsd, FPRE, FPREsd, NCELLS = runs(agentsimEM, 100)

times = _times.*90.0./56.0
p = plot(times, FBP.*100, yerror=FBPsd.*100, color=:purple, label="BP")
plot!(times, FEPI.*100, yerror=FEPIsd.*100, color=:red, label="EPI")
plot!(times, FPRE.*100, yerror=FPREsd.*100, color=:blue, label="PRE")
plot!(times, ones(length(FBP)).*50.0, color=:black, linestyle=:dash, leg=true, label="50 %")
ylabel!("% of ICM")
xlabel!("scaled time")
xlims!(60,100)
ylims!(0,100)
#save_fig(p, "/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel/Figures/cellfates_vs_time")

p = plot(mean(NCELLS,dims=2), FBP.*100, yerror=FBPsd.*100, color=:purple, label="BP")
plot!(mean(NCELLS,dims=2), FEPI.*100, yerror=FEPIsd.*100, color=:red, label="EPI")
plot!(mean(NCELLS,dims=2), FPRE.*100, yerror=FPREsd.*100, color=:blue, label="PRE")
plot!(mean(NCELLS,dims=2), ones(length(FBP)).*50.0, color=:black, linestyle=:dash, leg=true, label="50 %")
ylabel!("% of ICM")
xlabel!("ICM cell count")
ylims!(0,100)
#save_fig(p, "/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel/Figures/cellfates_vs_cellN")
