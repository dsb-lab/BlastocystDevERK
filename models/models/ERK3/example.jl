if pwd() !== "/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/ERK3"
    cd("/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/ERK3")
end

using PyCall
using PyPlot
include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("../../utils/plot_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")
path_figures="/home/pablo/Desktop/PhD/projects/BlastocystDev/figures/ERK_model_4/"

### SIMULATION ###
@time results = agentsimICM( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.001  #integration time-step
                          , mh=0.5   #data measuring time-step
                          , Fth=Inf
                        #   , comext=1.15
                          , comext=0.0
                          , comKO=false);
### END SIMULATION ###

### PLOTING ###
PyPlot.close("all")

# Plot example of time traces on a 3x3 grid
plot_timeseries(results, all=false)
# path_figures*"$cfr.svg"
# plot 3d embryo 

for cn in unique(results.totals)
  plot_3D_embryo(results, cell_number=cn, path=path_figures*"progression/$cn.svg")
end
gcf()
# plot cell fates over cell number
plot_fates_cellN(results)

ERK_EPI = var_distribution(results, "ERK", 1, N_start)
ERK_PrE = var_distribution(results, "ERK", 2, N_start)

fig,ax = plt.subplots()
ax.scatter(rand(length(ERK_EPI)), ERK_EPI, c="red",s=1)
ax.scatter(rand(length(ERK_PrE)).+2.0, ERK_PrE, c="blue",s=1)
ax.set_ylabel("ERK")
ax.set_ylim(0.5,2.2)
ax.set_xticks([])
gcf()

fig,ax = plt.subplots()
plt.hist(ERK_PrE, color=[0.9,0,0.9,0.5], bins=30, density=true)
plt.hist(ERK_EPI, color=[0,0.8,0,0.5], bins=30, density=true)

plt.ylabel("pERK")
gcf()

ERK_EPI = Vector{Float64}()
ERK_PrE = Vector{Float64}()
N=10
for n in range(0, N)
    @time results = agentsimICM( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.0005  #integration time-step
                          , mh=0.5   #data measuring time-step
                          , Fth=Inf
                        #   , comext=1.15
                          , comext=0.0
                          , comKO=false);

    _ERK_EPI = var_distribution(results, "ERK", 1, N_start)
    _ERK_PrE = var_distribution(results, "ERK", 2, N_start)
    ERK_EPI = [ERK_EPI; _ERK_EPI]
    ERK_PrE = [ERK_PrE; _ERK_PrE]
end
path_figures="/home/pablo/Desktop/PhD/projects/BlastocystDev/figures/ERK_model_4/"

# Calculate the point density
using KernelDensity
using LaTeXStrings
using Random
using PyPlot
using PyCall
PyPlot.matplotlib[:rc]("mathtext",fontset="dejavusans")
PyPlot.matplotlib[:rc]("font",size=18)

fig,ax = plt.subplots(figsize=(5,4))

kde_epi = kde(ERK_EPI)
plt.plot(kde_epi.x, kde_epi.density, color=[0,0.7,0,1.0],lw=5)
plt.hist(ERK_EPI, color=[0,0.8,0,0.3], bins=30, density=true)

kde_pre = kde(ERK_PrE)
plt.plot(kde_pre.x, kde_pre.density, color=[0.8,0,0.8,1.0],lw=5)
plt.hist(ERK_PrE, color=[0.9,0,0.9,0.3], bins=30, density=true)

plt.xlabel("pERK")
plt.tight_layout()
plt.savefig(join([path_figures, "erkdists.svg"]))
gcf()
