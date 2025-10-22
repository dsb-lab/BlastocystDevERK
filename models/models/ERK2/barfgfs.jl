if pwd() !== "/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/ERK2"
    cd("/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models/ERK2")
end
using LaTeXStrings
using Random
using PyPlot
using PyCall
PyPlot.matplotlib[:rc]("mathtext",fontset="dejavusans")        #computer modern font 
PyPlot.matplotlib[:rc]("font",size=18)
@pyimport matplotlib.patches as patch

include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("../../utils/plot_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")
include("utils/ko_utils.jl")

path_figures = "/home/pablo/Desktop/PhD/projects/BlastocystDev/figures/ERK_model_1/"

nsimulations = 20

labs_step=0.05
labs_end = 1.3
labs_start = 0.8
labels = collect(range(labs_start, step=labs_step, stop=labs_end))
labels = [labels; 0; 0]
fatesDPs  = zeros(length(labels), nsimulations)
fatesEPIs = zeros(length(labels), nsimulations)
fatesPrEs = zeros(length(labels), nsimulations)

### SIMULATION ###
@time simresults = agentsimICM( ERK_model_1, ["NANOG", "GATA6", "FGF", "ERK"], "FGF"
                    , h=0.001  #integration time-step
                    , mh=0.5   #data measuring time-step
                    , Fth=Inf
                    , comext=0.0
                    , comKO=false);


for (i,FGFmedium) in enumerate(labels)
    if i==length(labels)-1
        continue
    end
    println("\n######\n $FGFmedium \n######\n")
    Threads.@threads for n=1:nsimulations
        @time results = agentsimICM( ERK_model_1, ["NANOG", "GATA6", "FGF", "ERK"], "FGF"
        , h=0.001  #integration time-step
        , mh=0.01   #data measuring time-step
        , Fth=Inf
        , comext=FGFmedium
        , comKO=true);
        println(n)

        _NANOG = results.vars[1]
        _GATA6 = results.vars[2]
        _FGF = results.vars[3]

        _fBP = results.fDP
        _fEPI = results.fEPI
        _fPRE = results.fPRE
        _times = results.times

        times_, fBP, fEPI, fPRE = process_results(_times, _fBP, _fEPI, _fPRE)
        _totals = fBP .+ fEPI .+ fPRE
        _totals = _totals[1:length(times_)]
        
        _CFATES = results.CFATES
        _X = results.positions.X
        _Y = results.positions.Y
        _Z = results.positions.Z

        idx_start = find_start(_totals, N_start)-100
        NANOGsim = _NANOG[:,idx_start:length(times_)]
        GATA6sim = _GATA6[:,idx_start:length(times_)]
        FGF = _FGF[:,idx_start:length(times_)]
        CFATES = _CFATES[:,idx_start:length(times_)]
        X = _X[:,idx_start:length(times_)]
        Y = _Y[:,idx_start:length(times_)]
        Z = _Z[:,idx_start:length(times_)]
        times  = times_[idx_start:end]
        totals = _totals[idx_start:length(times_)]
        fatesDPs[i,n]  = fBP[end]./_totals[end]
        fatesEPIs[i,n] = fEPI[end]./_totals[end]
        fatesPrEs[i,n] = fPRE[end]./_totals[end]
    end
end

# Means
fatesDPms   = reshape(mean(fatesDPs, dims=2), length(labels))
fatesEPIms  = reshape(mean(fatesEPIs, dims=2), length(labels))
fatesPrEms  = reshape(mean(fatesPrEs, dims=2), length(labels))

# Stds
fatesDPstds  = reshape(std(fatesDPs, dims=2), length(labels))
fatesEPIstds = reshape(std(fatesEPIs, dims=2), length(labels))
fatesPrEstds = reshape(std(fatesPrEs, dims=2), length(labels))

width = 0.04 # the width of the bars: can also be len(x) sequence
PyPlot.close("all")
fig, ax = plt.subplots(figsize=(7,6)) 
labels[end] = labs_end + 2* labs_step
ax.bar(labels, fatesEPIms.+fatesPrEms.+fatesDPms, width, label="EPI", color=[0.0, 0.8, 0.0], edgecolor="black")
ax.bar(labels, fatesPrEms.+fatesDPms, width, label="PrE", color=[0.8, 0.0, 0.8], edgecolor="black")
ax.bar(labels, fatesDPms, width, label="DP", color=[0.5, 0.5, 0.5], edgecolor="black")
ax.set_xlim(labels[1]-labs_step, labels[end]+labs_step)
ax.set_ylabel("Fraction of ICM")
ax.set_xlabel(L" $F_m$")
ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)
#ax.legend()
# PyPlot.tight_layout()
ax.set_xticks([0.8, 1.0, 1.2, 1.4])
ax.set_xticklabels([0.8, 1.0, 1.2, "MEKi"])
name=join([path_figures, "fates_vs_FGF", ".svg"])
PyPlot.tight_layout()
savefig(name)
gcf()

