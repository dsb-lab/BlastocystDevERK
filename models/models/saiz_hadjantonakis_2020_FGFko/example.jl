if pwd() !== "/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/saiz_hadjantonakis_2020_FGFko"
    cd("/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/saiz_hadjantonakis_2020_FGFko")
end
include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("../../utils/plot_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")
include("utils/ko_utils.jl")

### SIMULATION ###W
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
                          , comext=3.5
                          , comKO=true);
### END SIMULATION ###

### PLOTING ###
PyPlot.close("all")

# Plot example of time traces on a 3x3 grid
plot_timeseries(results,all=false)

# plot 3d embryo 
plot_3D_embryo(results)

# plot cell fates over cell number
plot_fates_cellN(results)

