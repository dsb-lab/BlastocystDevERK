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

path_figures = "/home/pablo/Desktop/PhD/projects/BlastocystDev/figures/saiz_FGFko/"
fgfs = [2.0,4.0,6.0]
ms = [1.5,4.0]

const coefm=ms[2]
for FGFmedium in fgfs
    Nanogmax=10.0

    Gata6max = ag*FGFmedium/dg

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
                            , comext=FGFmedium
                            , comKO=true);
    ### END SIMULATION ###
    @time Nrange, Grange, N_Nnull, G_Gnull = NGF_nullclines(FGFmedium, Nanogmax, Gata6max);
    @time XX, YY, dx, dy = NGF_vfield(FGFmedium, Nanogmax, Gata6max);
    _NANOG = results.vars[1]
    _GATA6 = results.vars[2]
    _FGF = results.vars[3]
    positions=results.positions
    @time fBP, fEPI, fPRE, _totals, NANOG, GATA6, FGF, CFATES, X, Y, Z, times, totals = NGF_postprocessing(results.times, results.fDP, results.fEPI, results.fPRE, _NANOG, _GATA6, _FGF, positions.X, positions.Y, positions.Z, results.R, results.CFATES);

    if true
    @time init_pointsEPIN, init_pointsEPIG, init_pointsPrEN, init_pointsPrEG, end_pointsEPIN, end_pointsEPIG, end_pointsPrEN, end_pointsPrEG = fate_phasespace(FGFmedium, Nanogmax, Gata6max, resolution=100);
    end

    using LaTeXStrings
    using Random
    using PyPlot
    using PyCall
    PyPlot.matplotlib[:rc]("mathtext",fontset="dejavusans")        #computer modern font 
    PyPlot.matplotlib[:rc]("font",size=19)
    @pyimport matplotlib.patches as patch

    PyPlot.close("all")
    PyPlot.clf()
    fig, ax = plt.subplots(figsize=(7,6)) # Create a new blank figure
    ax.plot(G_Gnull, Nrange, label="G null", linewidth=4, zorder=20)
    ax.plot(Grange, N_Nnull, label="N null", linewidth=4, zorder=20)
    ax.set_xlabel("GATA6")
    ax.set_ylabel("NANOG")
    xtran = 0.1*Nanogmax
    xtrag = 0.1*Gata6max
    ax.set_xlim(-xtrag,Gata6max+xtrag)
    ax.set_ylim(-xtran,Nanogmax+xtran)

    ax.streamplot(XX, YY, dx, dy,
                arrowsize=2,arrowstyle="->", linewidth=1, color=(0.2,0.2,0.2))

    ax.scatter([], [], color=[0.0, 0.8, 0.0,0.2], s=120, label="EPI basin of attraction")
    ax.scatter([], [], color=[0.8, 0.0, 0.8,0.2], s=120, label="PrE basin of attraction")
    ax.scatter([], [], color="grey", edgecolor="black", s=80, label="Initial conditions")

    idx = findall(init_pointsEPIG.!=0.0)
    ax.scatter([init_pointsEPIG[idx]], [init_pointsEPIN[idx]], color=(0.0, 0.8, 0.0,0.2), s=20.0)
    idx = findall(init_pointsPrEG.!=0.0)
    ax.scatter([init_pointsPrEG[idx]], [init_pointsPrEN[idx]], color=(0.8, 0.0, 0.8,0.2), s=20.0)

    idxs = findall(NANOG.<0.0001)
    NANOG[idxs] .= ones(length(idxs))*.-10
    idxs = findall(GATA6.<0.0001)
    GATA6[idxs] .= ones(length(idxs))*.-10

    ax.scatter(GATA6[1:N_start-1,1], NANOG[1:N_start-1,1], color="grey", edgecolor="black", s=80.0, zorder=5)
    if sum(end_pointsEPIG) != 0
        idx = findfirst(0.0 .!= end_pointsEPIG)
        ax.scatter([end_pointsEPIG[idx]], [end_pointsEPIN[idx]],  color=[0.0, 0.8, 0.0,1.0], edgecolor="black", s=120,zorder=25, label="Epi FP")
    end
    if sum(end_pointsPrEG) != 0
        idx = findfirst(0.0 .!= end_pointsPrEG)
        ax.scatter([end_pointsPrEG[idx]], [end_pointsPrEN[idx]], color=[0.8, 0.0, 0.8,1.0], edgecolor="black", s=120, zorder=25, label="PrE FP")
    end

    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    # PyPlot.legend(bbox_to_anchor=[1.0,0.75])
    plt.title(join([L" $F_m$", "=$FGFmedium", L" ; $m=$","$coefm"]))
    # plt.title(join([L" $F_m$", "=$FGFmedium"]))
    PyPlot.tight_layout()

    name=join([path_figures,"initialconditions_FGFext$FGFmedium$coefm", ".svg"])
    # name=join([path_figures,"finalconditions_FGFext$FGFmedium", ".png"])
    savefig(name)
end
