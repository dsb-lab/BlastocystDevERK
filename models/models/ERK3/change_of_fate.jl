if pwd() !== "/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/ERK3"
    cd("/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/ERK3")
end

include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
# include("../../utils/plot_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")

N = 20

FATE_CHANGE_PER_CELL = zeros(N, 3)
### SIMULATION ###
@time results = agentsimICM( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.001  #integration time-step
                          , mh=0.1   #data measuring time-step
                          , Fth=Inf
                        #   , comext=1.15
                          , comext=0.0
                          , comKO=false);


cfrs = [1.2, 2.2, 3.2]
for ii in eachindex(cfrs)
    test_cfr =  cfrs[ii]
    for n=1:N
        @time _results = agentsimICM(ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                            , h=0.001  #integration time-step
                            , mh=0.1   #data measuring time-step
                            , Fth=Inf
                            #   , comext=1.15
                            , comext=0.0
                            , comKO=false
                            , test_cfr=test_cfr);
        fc_per_cell = zeros(N_start)
        for c=1:N_start
            fc_per_cell[c]+=sum(diff(_results.CFATES[c,:]).>0)
        end
        FATE_CHANGE_PER_CELL[n, ii] = mean(fc_per_cell)
    end
end

path_figures="/home/pablo/Desktop/PhD/projects/BlastocystDev/figures/ERK_model_4/"
using DelimitedFiles
writedlm(path_figures*"fgfr2_KO_all", FATE_CHANGE_PER_CELL)

using LaTeXStrings
using Random
using PyPlot
using PyCall
PyPlot.matplotlib[:rc]("mathtext",fontset="dejavusans") 
PyPlot.matplotlib[:rc]("font",size=20)


f_change = readdlm(path_figures*"fgfr2_all")
f_change_KO = readdlm(path_figures*"fgfr2_KO_all")

f_change_mean = mean(f_change, dims=1)
f_change_mean_KO = mean(f_change_KO, dims=1)
f_change_std = std(f_change, dims=1)
f_change_std_KO = std(f_change_KO, dims=1)

fig, ax = plt.subplots(figsize=(5.5,4.5))
barwidth = 0.5
colors = [[0.6, 0.6, 0.6], [255/255.0, 100.0/255.0, 0.0]]

cfrs = [1.2, 2.2, 3.2]

Y = vec(f_change_mean)
Z =  vec(f_change_mean_KO)

Y_std = vec(f_change_std)
Z_std =  vec(f_change_std_KO)

x = collect(range(1, length(f_change_mean)).* 1.5)
for i in 1:2
    println(i)
    if i==1
        offset=-barwidth/2
        data=Y
        data_std = Y_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-WT")
        ax.errorbar(x .+ offset, data, yerr=data_std, capsize=5, fmt="o", color="k")

    elseif i == 2
        offset = +barwidth/2
        data=Z
        data_std = Z_std
        rects = ax.bar(x .+ offset, data, barwidth, color=colors[i], edgecolor="k", label="FGFR2-KO")
        ax.errorbar(x .+ offset,  data, yerr=data_std, capsize=5, fmt="o", color="k")
    end
end

ax.set_xticks(x)
ax.set_xticklabels(cfrs)
ax.set_ylabel(L"$F_{\mathrm{change}}$ / cell")
ax.set_xlabel(L"$f_{range}$")
savefig("$path_figures coupling_range.svg")

using HypothesisTests
f_change = readdlm(path_figures*"fgfr2_all")[:, 1]
f_change_KO = readdlm(path_figures*"fgfr2_KO_all")[:, 1]
res = EqualVarianceTTest(f_change, f_change_KO)
pvalue(res)