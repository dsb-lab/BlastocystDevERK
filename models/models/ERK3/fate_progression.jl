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
                          , mh=0.01   #data measuring time-step
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

  # plot cell fates over cell number
plot_fates_cellN(results)
  
times = []
for cn in unique(results.totals)
  plot_3D_embryo(results, cell_number=cn, path=path_figures*"progression/$cn.svg")
  idx = findfirst(results.totals.>=cn)
  push!(times, idx)
end

N=10
ALL = zeros(N, length(times), 3, 3)
for nn=1:N
    @time results = agentsimICM( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.001  #integration time-step
                          , mh=0.01   #data measuring time-step
                          , Fth=Inf
                        #   , comext=1.15
                          , comext=0.0
                          , comKO=false);
                            
    times = []
    for cn in unique(results.totals)
        idx = findfirst(results.totals.>=cn)
        push!(times, idx)
    end
    for (tid, t) in enumerate(times)
        ### END SIMULATION ###
        cfates = results.CFATES[:,t]
        xi = results.positions.X[:, t]
        yi = results.positions.Y[:, t]
        zi = results.positions.Z[:, t]
        r = results.R[:, t]

        fates = zeros(3)
        for i=1:results.totals[t]
            fates[results.CFATES[i, t]+1]+=1
        end
        if fates[1]==0
            fates[1]=1
        end
        cell_neighs = zeros(3, results.totals[t])
        for i=1:results.totals[t] # Lower triangular excluding diagonal
            cell_neigh = zeros(3)
            for j=1:results.totals[t]
                if i==j
                    continue
                end
                d    = sqrt((xi[i] - xi[j])^2 + (yi[i]-yi[j])^2 + (zi[i]-zi[j])^2)
                id   = 1.0/d
                rij  = r[i]+r[j]
                if (d<(1.2*rij))
                    cell_neigh[results.CFATES[j, t]+1] += 1
                end
            end
            cell_neighs[:, i] .= cell_neigh./fates
            cell_neighs[:, i] ./= sum(cell_neighs[:, i])
        end

        for fate=1:3
            ids = findall(results.CFATES[1:results.totals[t], t].==fate-1)
            if isempty(ids)
                continue
            end
            cell_neighs_fate = cell_neighs[:, ids]
            ALL[nn, tid, fate, :] .+= mean(cell_neighs_fate, dims=2)
        end
    end
end

using NaNStatistics
data = nanmean(ALL,dims=1)[1,:,:,:]
data_std = nanstd(ALL,dims=1)[1,:,:,:]

fig, ax = plt.subplots(1,2,figsize=(14,5))
for i=2:3
ax[i-1].plot(unique(results.totals), data[:,i,1], color=[0.5,0.5,0.5])
ax[i-1].fill_between(unique(results.totals), data[:,i,1]-data_std[:,i,1],data[:,i,1]+data_std[:,i,1], color=[0.5,0.5,0.5], alpha=0.3)
ax[i-1].plot(unique(results.totals), data[:,i,2], color=[0.0,0.8,0.0])
ax[i-1].fill_between(unique(results.totals), data[:,i,2]-data_std[:,i,2],data[:,i,2]+data_std[:,i,2], color=[0.0,0.8,0.0], alpha=0.3)
ax[i-1].plot(unique(results.totals), data[:,i,3], color=[0.8,0.0,0.8])
ax[i-1].fill_between(unique(results.totals), data[:,i,3]-data_std[:,i,3],data[:,i,3]+data_std[:,i,3], color=[0.8,0.0,0.8], alpha=0.3)
ax[i-1].set_xlim(N_start,Nmax)
end
gcf()
