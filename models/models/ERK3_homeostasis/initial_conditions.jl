
cwd = pwd()
foldername = basename(cwd)
basepath = dirname(dirname(cwd))
save_dir = joinpath(basepath, "results", foldername)
mkpath(save_dir)

include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")

import PyPlot
using PyCall

PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",size=40)
@pyimport matplotlib.patches as patch
### SIMULATION ###
@time results, ESC_cell_ids = agentsimICM_ESCs( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                          , h=0.001  #integration time-step
                          , mh=0.01   #data measuring time-step
                          , Fth=Inf
                        #   , comext=1.15
                          , comext=0.0
                          , comKO=false
                          , NESCs=10
                          , fixedFGF=2.0);
### END SIMULATION ###
results.totals

### PLOTING ###
PyPlot.close("all")

for ESC_Number in [0, 2, 5, 7, 10]
    nanogs_init = Vector{Float64}()
    gata6s_init = Vector{Float64}()
    colors = []
    N=100
    for n in range(1, N)
    @time _results, _ESC_cell_ids = agentsimICM_ESCs( ERK_model_3, ["NANOG", "GATA6", "FGF", "ERK", "FGFR2"], "FGF"
                            , h=0.001  #integration time-step
                            , mh=0.01   #data measuring time-step
                            , Fth=Inf
                            #   , comext=1.15
                            , comext=0.0
                            , comKO=false
                            , NESCs=ESC_Number
                            , fixedFGF=2.0);

    for c in range(2,N_start)
        if c in _ESC_cell_ids
        continue
        end
        push!(nanogs_init, _results.vars[1][c, :][findfirst(_results.vars[1][c, :].>0)])
        push!(gata6s_init, _results.vars[2][c, :][findfirst(_results.vars[2][c, :].>0)])
        if _results.CFATES[c,end] == 1
        push!(colors, "red")
        elseif _results.CFATES[c,end] == 2
        push!(colors, "blue")
        else
        push!(colors, "purple")
        end
    end
    end

    # Start with a square Figure.
    fig = PyPlot.figure(figsize=(10, 10))
    # fig.suptitle("% ESCs = 25", fontsize=35)

    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal Axes and the main Axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(2, 2,  width_ratios=(3, 1), height_ratios=(1, 3),
                        left=0.1, right=0.9, bottom=0.1, top=0.9,
                        wspace=0.05, hspace=0.05)
    # Create the Axes.
    ax = fig.add_subplot(gs[2, 1])
    ax_histx = fig.add_subplot(gs[1, 1], sharex=ax)
    ax_histy = fig.add_subplot(gs[2, 2], sharey=ax)
    # Draw the scatter plot and marginals.
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=false)
    ax_histy.tick_params(axis="y", labelleft=false)

    # Calculate the point density

    sc = ax.scatter(gata6s_init, nanogs_init, c=colors, s=50, alpha=0.3)

    # PyPlot.colorbar(sc, ax=ax, label="Number of points per pixel")
    ax.set_xlabel("Gata6")
    ax.set_ylabel("Nanog")

    ax.set_xlim(0.6, 2.0)
    ax.set_ylim(0.6, 2.0)
    ax.set_xticks(range(0.8, 1.9, step=0.5))
    ax.set_yticks(range(0.8, 1.9, step=0.5))
    PyPlot.tight_layout()
    # now determine nice limits by hand:

    gata6_pre = [gata6s_init[i] for i in eachindex(gata6s_init) if colors[i]=="blue"]
    _n, _bins, _ = ax_histx.hist(gata6_pre, bins=30, color="blue", alpha=0.3, density=true, label="PrE")
    nanog_pre = [nanogs_init[i] for i in eachindex(nanogs_init) if colors[i]=="blue"]
    _n, _bins, _ = ax_histy.hist(nanog_pre, bins=30, orientation="horizontal", color="blue", alpha=0.3, density=true)

    gata6_epi = [gata6s_init[i] for i in eachindex(gata6s_init) if colors[i]=="red"]
    _n, _bins, _ = ax_histx.hist(gata6_epi, bins=30, color="red", alpha=0.3, density=true, label="Epi")
    nanog_epi = [nanogs_init[i] for i in eachindex(nanogs_init) if colors[i]=="red"]
    _n, _bins, _ = ax_histy.hist(nanog_epi, bins=30, orientation="horizontal", color="red", alpha=0.3, density=true)
    ax_histy.set_xticks([])
    ax_histx.set_yticks([])
    name=join([save_dir, "/$ESC_Number", ".pdf"])
    fig.savefig(name)
end