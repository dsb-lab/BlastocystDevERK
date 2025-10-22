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


FGFmedium = 3.5
const Nanogmax = 10.0
if FGFmedium < 0.5
    const Gata6max = 10.0
else
    const Gata6max = ag*FGFmedium/dg
end

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
@time Nrange, Grange, N_Nnull, G_Gnull = NGF_nullclines(FGFmedium);
@time XX, YY, dx, dy = NGF_vfield(FGFmedium);
_NANOG = results.vars[1]
_GATA6 = results.vars[2]
_FGF = results.vars[3]
positions=results.positions
@time fBP, fEPI, fPRE, _totals, NANOGsim, GATA6sim, FGFsim, CFATES, X, Y, Z, times, totals = NGF_postprocessing(results.times, results.fDP, results.fEPI, results.fPRE, _NANOG, _GATA6, _FGF, positions.X, positions.Y, positions.Z, results.R, results.CFATES);

if true
@time init_pointsEPIN, init_pointsEPIG, init_pointsPrEN, init_pointsPrEG, end_pointsEPIN, end_pointsEPIG, end_pointsPrEN, end_pointsPrEG = fate_phasespace(FGFmedium, resolution=5);
end

using Random
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch
@pyimport matplotlib.animation as animation
plt.rcParams["animation.ffmpeg_path"] = "/usr/bin/ffmpeg"


PyPlot.close("all")
fig, ax = PyPlot.subplots(figsize=(18,12))
ax.plot(G_Gnull, Nrange, label="G null", linewidth=4, zorder=1)
ax.plot(Grange, N_Nnull, label="N null", linewidth=4, zorder=2)
ax.set_xlabel("GATA6")
ax.set_ylabel("NANOG")
xtran = 0.2*Nanogmax
xtrag = 0.2*Gata6max
ax.set_xlim(-xtrag/2,Gata6max+xtrag)
ax.set_ylim(-xtran/2,Nanogmax+xtran)

x = range(0, length=20, stop=Gata6max.+0.2*Gata6max)
y = range(0, length=20, stop=Nanogmax.+0.2*Nanogmax)
X = x' .* ones(length(y))
Y = ones(length(x))' .* y

dx = ag.*FGFmedium./(1.0 .+ (Y./Kn).^coefn) .- dg.*X
dy = an./(1.0 .+ (FGFmedium.*X./Kg).^coefm) .- dn.*Y
ax.streamplot(X, Y, dx, dy,
               arrowsize=1,arrowstyle="->", linewidth=2, color=(0.2,0.2,0.2))

ax.scatter([], [], color="red", s=3, label="EPI starting points")
#ax.plot([],[], color="red", linewidth=1,label="EPI trajectories")
ax.scatter([], [], color="blue", s=3, label="PrE starting points")
#ax.plot([],[], color="blue", linewidth=1,label="PrE trajectories")
ax.scatter([], [], color="black", s=3, label="FP")

ax.scatter([init_pointsEPIG], [init_pointsEPIN], color=(1.0,0.0,0.0,0.3), s=4.0)
ax.scatter([init_pointsPrEG], [init_pointsPrEN], color=(0.0,0.0,1.0,0.3), s=4.0)

sc = ax.scatter(GATA6sim[:,1], NANOGsim[:,1], color="green", s=60.0)
if sum(end_pointsEPIG) != 0
    idx = findfirst(0.0 .!= end_pointsEPIG)
    ax.scatter([end_pointsEPIG[idx]], [end_pointsEPIN[idx]], color="black", s=60,zorder=3)
end
if sum(end_pointsPrEG) != 0
    idx = findfirst(0.0 .!= end_pointsPrEG)
    ax.scatter([end_pointsPrEG[idx]], [end_pointsPrEN[idx]], color="black", s=60, zorder=4)
end
ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)
PyPlot.legend(bbox_to_anchor=[1.0,1.0])
PyPlot.title("FGF in medium = $FGFmedium")
PyPlot.tight_layout()
gcf()
function update(frame)
    println(frame)
    sc._offsets = hcat(GATA6sim[:,frame], NANOGsim[:,frame])
    PyPlot.draw()
end
anim = animation.FuncAnimation(fig, update, vcat(10:5:300 ,305:20:length(GATA6sim[1,:])), repeat=false, blit=false)
FFwriter=animation.FFMpegWriter(fps=1, bitrate=-1,extra_args=["-vcodec", "libx264"])
name=join(["anims/anim_FGF$FGFmedium", ".gif"])
anim.save("../../figures/NGF/FGFko/phaseplane/$name", fps=10)
