using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

function plot_3D_embryo(simresults::SimResults; cell_number=nothing, path=nothing)
    if cell_number===nothing
        cell_number=Nmax
    end
    print(cell_number)

    idx = findfirst(simresults.totals.>=cell_number)
    colrs = []
    for i=1:length(simresults.CFATES[:,idx])
        if simresults.CFATES[i, idx] == 1 #EPIa
            push!(colrs, [0.0, 0.8, 0.0])
        end
        if simresults.CFATES[i, idx] == 2 #PrE
            push!(colrs, [0.8,0.0,0.8])
        end
        if simresults.CFATES[i, idx] == 0 #DP
            push!(colrs,  [0.5,0.5,0.5])
        end
    end

    x = Vector{Float64}()
    y = Vector{Float64}()
    z = Vector{Float64}()
    for i=1:length(simresults.positions.X[:,idx])
        push!(x, simresults.positions.X[i,idx])
        push!(z, simresults.positions.Y[i,idx])
        push!(y, simresults.positions.Z[i,idx])
    end
    fig_3d, ax_3d = plt.subplots()
    scatter3D(x, y, z, s=2000, c=colrs)
    # Hide grid lines
    plt.grid(false)

    # Hide axes ticks
    plt.axis("off")
    plt.tight_layout()
    if path===nothing
        gcf()
    else
        gcf()
        savefig(path)
    end

    return fig_3d, ax_3d
end

function plot_fates_cellN(simresults::SimResults)
    fig2, ax1 = plt.subplots()
    ax1.plot(simresults.totals, simresults.fDP./simresults.totals, color="purple", label="BP")
    ax1.plot(simresults.totals, simresults.fEPI./simresults.totals, color="red", label="EPI")
    ax1.plot(simresults.totals, simresults.fPRE./simresults.totals, color="blue", label="PRE")
    ax1.set_ylim(-0.1,1.1)
    ax1.set_ylabel("% of ICM")
    xlabel("ICM cell count")
    p_patch = patch.Patch(color="purple", label="BP")
    r_patch = patch.Patch(color="red", label="EPI")
    b_patch = patch.Patch(color="blue", label="PrE")
    plt.legend(handles=[p_patch, r_patch, b_patch], loc=1)
    plt.tight_layout
    return fig2, ax1
end

function plot_timeseries(simresults::SimResults; all=false)
    fig1, ax = plt.subplots(3,3, figsize=(10,8), sharey=false, sharex=true)
    rang = collect(Int64, range(1,step=1,stop=Nmax))
    k=[0.0]
    totplt   = simresults.totals./maximum(simresults.totals)
    clrs = ["red", "blue", "green", "yellow", "grey", "orange"]
    comvaridx = findfirst(simresults.comvarname .== simresults.varsnames)
    comvarname = simresults.comvarname
    max_val = 0.0
    min_val = Inf
    for i=1:3
        for j=1:3
            jj = rand(rang)
            if i==1
                if j==1
                    jj = 1
                end
            end
            if all 
                start = 1 
            else 
                start = find_start(simresults.totals, N_start)
            end
            for vid in eachindex(simresults.vars)
                max_val = max(max_val, maximum(simresults.vars[vid][jj,start:end]))
                min_val = min(min_val, minimum(simresults.vars[vid][jj,start:end]))
                ax[i,j].plot(simresults.times[start:end], simresults.vars[vid][jj,start:end], label=simresults.varsnames[vid], color=clrs[vid])
            end
            ax[i,j].plot(simresults.times[start:end], simresults.comvar[jj,start:end], label="$comvarname received", color=clrs[comvaridx], linestyle="--")
            ax2 = ax[i,j].twinx()
            ax2.plot(simresults.times[start:end], totplt[start:end], color="black")
            k[1]+=1
            if k[1] in [7,8,9]
                ax[i,j].set_xlabel("time")
            else
                ax[i,j].set_xticks([])
            end
            if k[1] in [1, 4, 7]
                ax[i,j].set_ylabel("VARIABLES")
            else
                ax[i,j].set_yticks([])
            end
            if k[1] in [3,6,9]
                ax2.set_ylabel("# cells")
            else
                ax2.set_yticks([])
            end
            if simresults.CFATES[jj,end] == 1
                ax[i,j].set_facecolor((1.0, 0.8, 0.8))
            elseif simresults.CFATES[jj,end] == 2
                ax[i,j].set_facecolor((0.8, 1.0, 1.0))
            end
            if k[1]==3
                handls = []
                for vid in eachindex(simresults.vars)
                    push!(handls, patch.Patch(color=clrs[vid], label=simresults.varsnames[vid]))
                end
                ax[i,j].legend(handles=handls, bbox_to_anchor=[1.15,1],loc=2,borderaxespad=0)
            end
            if k[1] == 9
                p_patch = patch.Patch(color=(1.0, 1.0, 1.0), label="DP")
                r_patch = patch.Patch(color=(1.0, 0.8, 0.8), label="EPI")
                b_patch = patch.Patch(color=(0.8, 1.0, 1.0), label="PrE")
                ax[i,j].legend(handles=[p_patch, r_patch, b_patch], bbox_to_anchor=[1.15,1],loc=2,borderaxespad=0)
            end
        end
    end
    for i=1:3
        for j=1:3
            ax[i,j].set_ylim(min_val, max_val)
        end
    end
    plt.tight_layout()
    gcf()
    return fig1, ax 
end

