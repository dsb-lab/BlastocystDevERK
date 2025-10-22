# Set path to the repository folder
include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("../../utils/plot_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")

### SIMULATION ###
@time _ = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
                    , h=0.1    #integration time-step
                    , mh=1000000000.0   #data measuring time-step
                    , Fth=10000.0
                    , comext=0.0
                    , comKO = false);

@time results = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
                          , h=0.001  #integration time-step
                          , mh=0.5   #data measuring time-step
                          , Fth=Inf
                          , comext=0.0
                          , comKO=false);
### END SIMULATION ###


cwd = pwd()
foldername = basename(cwd)
basepath = dirname(dirname(cwd))
save_dir = joinpath(basepath, "results", foldername)

### PLOTING ###
PyPlot.close("all")

# Plot example of time traces on a 3x3 grid
fig1, ax1 = plot_timeseries(results,all=false)
savepath = joinpath(save_dir, "timeseries.pdf")
fig1.savefig(savepath)

# plot 3d embryo 
fig_3d, ax_3d = plot_3D_embryo(results)
savepath = joinpath(save_dir, "3drender.pdf")
fig_3d.savefig(savepath)

# plot cell fates over cell number
fig2, ax2 = plot_fates_cellN(results)
savepath = joinpath(save_dir, "fate_progression_example.pdf")
fig2.savefig(savepath)

