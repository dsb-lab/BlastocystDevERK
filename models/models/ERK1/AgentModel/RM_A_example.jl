using PyPlot
using PyCall
@pyimport matplotlib.patches as patch
if pwd() !== "/home/pablo/Desktop/PhD/projects/ICMModels/Nestor/AgentModel/"
    println("setting wd to /home/pablo/Desktop/PhD/projects/ICMModels/Nestor/code/AgentModel")
    cd()
    cd("/home/pablo/Desktop/PhD/projects/ICMModels/Nestor/AgentModel/")
end

include("constants/constants_mechanical.jl")
include("constants/constants_biochemical_RM.jl")
include("utils/general_utils.jl")
include("utils/agentmodel_utils.jl")

@time _times, _fBP, _fEPI, _fPRE, _nanog, _X, _Y, _Z, _R, _CFATES = agentsimRM(h=0.001,mh = 0.05, simtime=300.0, keep=true)

times_, fBP, fEPI, fPRE = process_results(_times, _fBP, _fEPI, _fPRE)
_totals = fBP .+ fEPI .+ fPRE
_totals = _totals[1:length(times_)]

idx_start = find_start(_totals, N_start)

nanog = _nanog[:,idx_start:length(times_)]
CFATES = _CFATES[:,idx_start:length(times_)]
X = _X[:,idx_start:length(times_)]
Y = _Y[:,idx_start:length(times_)]
Z = _Z[:,idx_start:length(times_)]
times = times_[idx_start:end]
totals = _totals[idx_start:length(times_)]
times = times .*90.0./56.0

PyPlot.close("all")
fig1, ax1 = PyPlot.subplots(3,3, figsize=(30,30), sharey=false, sharex=true)
rang = collect(Int64, range(1,step=1,stop=50))
k=[0.0]

for i=1:3
    for j=1:3
        jj = rand(rang)
        ax1[i,j].plot(times, nanog[jj,:], label="NANOG", color="red")
        ax1[i,j].set_ylim(0, maximum(nanog))

        k[1]+=1
        if k[1] in [7,8,9]
            ax1[i,j].set_xlabel("time")
        end
        if k[1] in [1, 4, 7]
            ax1[i,j].set_ylabel("NANOG")
        end
        if CFATES[jj,end] == 1
            ax1[i,j].set_facecolor((1.0, 0.8, 0.8))
        elseif CFATES[jj,end] == 2
            ax1[i,j].set_facecolor((0.8, 1.0, 1.0))
        end
        if k[1]==3
            nanog_patch = patch.Patch(color="red", label="NANOG")
            ax1[i,j].legend(handles=[nanog_patch], bbox_to_anchor=[1.15,1],loc=2,borderaxespad=0)
        end
        if k[1] == 9
            p_patch = patch.Patch(color=(1.0, 1.0, 1.0), label="BP")
            r_patch = patch.Patch(color=(1.0, 0.8, 0.8), label="EPI")
            b_patch = patch.Patch(color=(0.8, 1.0, 1.0), label="PrE")
            PyPlot.legend(handles=[p_patch, r_patch, b_patch], bbox_to_anchor=[1.15,1],loc=2,borderaxespad=0)
            #ax[:set_position]([0.06,0.06,0.71,0.91])
        end
    end
end

fig2, ax2 = plt.subplots()

ax2.plot(_totals, fBP./_totals, color="purple", label="BP")
ax2.plot(_totals, fEPI./_totals, color="red", label="EPI")
ax2.plot(_totals, fPRE./_totals, color="blue", label="PRE")
ax2.set_ylim(-0.1,1.1)
ax2.set_ylabel("% of ICM")
xlabel("ICM cell count")
p_patch = patch.Patch(color="purple", label="BP")
r_patch = patch.Patch(color="red", label="EPI")
b_patch = patch.Patch(color="blue", label="PrE")
PyPlot.legend(handles=[p_patch, r_patch, b_patch], loc=1)

idx = round(Int, length(X[1,:])/1)
colors = Vector{String}()
for i=1:length(CFATES[:,idx])
    if CFATES[i, idx] == 1 #EPI
        push!(colors, "red")
    end
    if CFATES[i, idx] == 2 #PrE
        push!(colors, "blue")
    end
    if CFATES[i, idx] == 0 #DP
        push!(colors, "purple")
    end
end

x = Vector{Float64}()
y = Vector{Float64}()
z = Vector{Float64}()
for i=1:length(X[:,idx])
    push!(x, X[i,idx])
    push!(z, Y[i,idx])
    push!(y, Z[i,idx])
end
fig_3d, ax_3d = PyPlot.subplots()
scatter3D(x, y, z, s=2000, c=colors)
gcf()