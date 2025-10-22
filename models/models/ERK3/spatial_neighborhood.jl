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

### SIMULATION ###
@time results = agentsimICM( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.0005  #integration time-step
                          , mh=0.5   #data measuring time-step
                          , Fth=Inf
                        #   , comext=1.15
                          , comext=0.0
                          , comKO=false);
### END SIMULATION ###

results.totals
DATA_MEANS = []
DATA_STDS = []

N = 20

test_cfr =1.2
ALL = zeros(N, 3, 3)
Threads.@threads for nn=1:N
    @time results = agentsimICM( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                            , h=0.001  #integration time-step
                            , mh=0.5   #data measuring time-step
                            , Fth=Inf
                            #   , comext=1.15
                            , comext=0.0
                            , comKO=false
                            , test_cfr=test_cfr);

    t = length(results.CFATES[1,:])
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
    println(fates)
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
        ALL[nn, fate, :] .+= mean(cell_neighs_fate, dims=2)
    end
end


ALL_mean = mean(ALL, dims=1)[1,:,:];
ALL_std = std(ALL, dims=1)[1,:,:];
push!(DATA_MEANS, ALL_mean)
push!(DATA_STDS, ALL_std)

test_cfr = 2.2
ALL = zeros(N, 3, 3)
Threads.@threads for nn=1:N

    @time results = agentsimICM( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                            , h=0.001  #integration time-step
                            , mh=0.5   #data measuring time-step
                            , Fth=Inf
                            #   , comext=1.15
                            , comext=0.0
                            , comKO=false
                            , test_cfr=test_cfr);

    t = length(results.CFATES[1,:])
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
        ALL[nn, fate, :] .+= mean(cell_neighs_fate, dims=2)
    end
end

ALL_mean = mean(ALL, dims=1)[1,:,:];
ALL_std = std(ALL, dims=1)[1,:,:];
push!(DATA_MEANS, ALL_mean)
push!(DATA_STDS, ALL_std)

DATA_MEANS[1]
DATA_MEANS[2]

using LaTeXStrings
using Random
using PyPlot
using PyCall
PyPlot.matplotlib[:rc]("mathtext",fontset="dejavusans") 
PyPlot.matplotlib[:rc]("font",size=20)

fig, ax = plt.subplots(figsize=(9,4.5))
barwidth = 0.5
colors = [[0.5, 0.5, 0.5], [0.0, 0.7, 0.0], [0.7, 0.0, 0.7]]

labs = ["Epi", "PrE", "Epi", "PrE"]

ALL_mean = DATA_MEANS[1]'
ALL_std = DATA_STDS[1]'

Y1 = vec(ALL_mean[1, 2:end])
Y2 = vec(ALL_mean[2, 2:end])
Y3 = vec(ALL_mean[3, 2:end])

Y1_std = vec(ALL_std[1, 2:end])
Y2_std = vec(ALL_std[2, 2:end])
Y3_std = vec(ALL_std[3, 2:end])

xx = collect(range(1, 4).* 1.5)
x = xx[1:2]
for i in 2:3
    if i == 2
        offset = -barwidth/2
        data=Y2
        data_std = Y2_std
        ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    elseif i == 3
        offset = +barwidth/2
        data=Y3
        data_std = Y3_std
        ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    end
end

ALL_mean = DATA_MEANS[2]'
ALL_std = DATA_STDS[2]'

Y1 = vec(ALL_mean[1, 2:end])
Y2 = vec(ALL_mean[2, 2:end])
Y3 = vec(ALL_mean[3, 2:end])

Y1_std = vec(ALL_std[1, 2:end])
Y2_std = vec(ALL_std[2, 2:end])
Y3_std = vec(ALL_std[3, 2:end])

x = xx[3:4]
for i in 2:3
    if i == 2
        offset = -barwidth/2
        data=Y2
        data_std = Y2_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    elseif i == 3
        offset = +barwidth/2
        data=Y3
        data_std = Y3_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    end
end

ax.set_xticks(xx)
ax.set_xticklabels(labs)
# ax.set_ylabel(L"neighborhood composition")
ax.set_xlabel("cell population")
path_figures="/home/pablo/Desktop/PhD/projects/BlastocystDev/figures/ERK_model_4/"
savefig("$path_figures neighborhood.svg")
gcf()



fig, ax = plt.subplots(figsize=(9,4.5))
barwidth = 0.5
colors = [[0.5, 0.5, 0.5], [0.0, 0.7, 0.0], [0.7, 0.0, 0.7]]

labs = ["Epi", "PrE", "Epi", "PrE"]

ALL_mean = DATA_MEANS[1]
ALL_std = DATA_STDS[1]

Y1 = vec(ALL_mean[1, 2:end])
Y2 = vec(ALL_mean[2, 2:end])
Y3 = vec(ALL_mean[3, 2:end])

Y1_std = vec(ALL_std[1, 2:end])
Y2_std = vec(ALL_std[2, 2:end])
Y3_std = vec(ALL_std[3, 2:end])

xx = collect(range(1, 4).* 2.0)
x = xx[1:2]
for i in 1:3
    if i == 1
        offset = -barwidth
        data=Y1
        data_std = Y1_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    elseif i == 2
        offset = 0
        data=Y2
        data_std = Y2_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    elseif i == 3
        offset = +barwidth
        data=Y3
        data_std = Y3_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    end
end

ALL_mean = DATA_MEANS[2]
ALL_std = DATA_STDS[2]

Y1 = vec(ALL_mean[1, 2:end])
Y2 = vec(ALL_mean[2, 2:end])
Y3 = vec(ALL_mean[3, 2:end])

Y1_std = vec(ALL_std[1, 2:end])
Y2_std = vec(ALL_std[2, 2:end])
Y3_std = vec(ALL_std[3, 2:end])

x = xx[3:4]
for i in 1:3
    if i == 1
        offset = -barwidth
        data=Y1
        data_std = Y1_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    elseif i == 2
        offset = 0
        data=Y2
        data_std = Y2_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    elseif i == 3
        offset = +barwidth
        data=Y3
        data_std = Y3_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    end
end

ax.set_xticks(xx)
ax.set_xticklabels(labs)
ax.set_ylabel(L"neighborhood composition")
ax.set_xlabel("cell population")
gcf()



# fig, ax = plt.subplots(1, 2, figsize=(9,4.5))

# ax1 = ax[1].twinx()
# ax1.plot(results.totals, color="black")
# ax[1].plot(ALL_mean[2, 2,:], color="green", lw=3)
# ax[1].fill_between(range(1, length(ALL_mean[2, 2,:])), ALL_mean[2, 2,:]-ALL_std[2, 2,:], ALL_mean[2, 2,:]+ALL_std[2, 2,:], color="grey", alpha=0.3)

# ax[1].plot( ALL_mean[2, 3,:], color="magenta", lw=3)
# ax[1].fill_between(range(1, length(ALL_mean[2, 3,:])), ALL_mean[2, 3,:]-ALL_std[2, 3,:], ALL_mean[2, 3,:]+ALL_std[2, 3,:], color="grey", alpha=0.3)

# ax[1].set_title("Epi neighborhood")
# ax[1].set_ylim(-0.1,1)
# ax2 = ax[2].twinx()
# ax2.plot(results.totals, color="black")
# ax[2].plot(ALL_mean[3, 2,:], color="green", lw=3)
# ax[1].fill_between(range(1, length(ALL_mean[3, 2,:])), ALL_mean[3, 2,:]-ALL_std[3, 2,:], ALL_mean[3, 2,:]+ALL_std[3, 2,:], color="grey", alpha=0.3)

# ax[2].plot(ALL_mean[3, 3,:], color="magenta", lw=3)
# ax[1].fill_between(range(1, length(ALL_mean[3, 3,:])), ALL_mean[3, 3,:]-ALL_std[3, 3,:], ALL_mean[3, 3,:]+ALL_std[3, 3,:], color="grey", alpha=0.3)

# ax[2].set_title("PrE neighborhood")
# ax[2].set_ylim(-0.1,1)

# id = findfirst(results.totals.>=N_start)
# ax[1].set_xlim(id, length(results.totals))
# ax[2].set_xlim(id, length(results.totals))
# tight_layout()
# gcf()