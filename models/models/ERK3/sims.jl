if pwd() !== "/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/test"
    cd("/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/test")
end

include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("../../utils/plot_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")
PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",size=35)
@pyimport matplotlib.patches as patch
### SIMULATION ###

### SIMULATION ###
@time results = agentsimICM( test, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.0005  #integration time-step
                          , mh=0.001   #data measuring time-step
                          , Fth=Inf
                          , comext=0.0
                          , comKO=false);

### END SIMULATION ###
rr = findfirst.(isequal.(unique(results.totals)), [results.totals])

N=10

tots_pre = results.totals[rr]
FDPmat  = zeros(N,length(tots_pre))
FEPImat = zeros(N,length(tots_pre))
FPREmat = zeros(N,length(tots_pre))
ERK_EPI = Vector{Float64}()
ERK_PrE = Vector{Float64}()
for n=1:N
    @time _results = agentsimICM( test, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.0005  #integration time-step
                          , mh=0.001   #data measuring time-step
                          , Fth=Inf
                          , comext=0.0
                          , comKO=false);

    rr = findfirst.(isequal.(unique(_results.totals)), [_results.totals])
    tots = _results.totals[rr]
    if length(tots)!=length(tots_pre)
        continue
    end
    FDPmat[n,:]  .= _results.fDP[rr]./tots
    FEPImat[n,:] .= _results.fEPI[rr]./tots
    FPREmat[n,:] .= _results.fPRE[rr]./tots
    _ERK_EPI = var_distribution(_results, "ERK", 1, N_start)
    _ERK_PrE = var_distribution(_results, "ERK", 2, N_start)
    ERK_EPI = [ERK_EPI; _ERK_EPI]
    ERK_PrE = [ERK_PrE; _ERK_PrE]
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
ax1.plot(tots_pre, FDP, color="purple", label="DP", linewidth=5)
ax1.fill_between(tots_pre, FDP.-FDPstd, FDP.+FDPstd, color="purple", alpha=0.3)
ax1.plot(tots_pre, FEPI, color="red", label="EPI", linewidth=5)
ax1.fill_between(tots_pre, FEPI.-FEPIstd, FEPI.+FEPIstd, color="red", alpha=0.3)
ax1.plot(tots_pre, FPRE, color="blue", label="PRE", linewidth=5)
ax1.fill_between(tots_pre, FPRE.-FPREstd, FPRE.+FPREstd, color="blue", alpha=0.3)

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
#plot(FPRE_plot.*100, linestyle="--", color="blue")
#plot(FEPI_plot.*100, linestyle="--", color="red")

#name=join(["fates_cellN_FGFext_$FGFmedium", "wexamples", ".png"])
#PyPlot.savefig("../../figures/NGF/FGFmedium/fates/$name")

gcf()

fig,ax = plt.subplots()
ax.hist(ERK_PrE, color=[0.9,0,0.9,0.5], bins=30, density=true)
ax.hist(ERK_EPI, color=[0,0.8,0,0.5], bins=30, density=true)

ax.set_ylabel("pERK")
gcf()

fig,ax = plt.subplots()
ax.scatter(rand(length(ERK_EPI)), ERK_EPI, c="red",s=1)
ax.scatter(rand(length(ERK_PrE)).+2.0, ERK_PrE, c="blue",s=1)
ax.set_ylabel("ERK")
# ax.set_ylim(0.5,2.2)
ax.set_xticks([])
gcf()