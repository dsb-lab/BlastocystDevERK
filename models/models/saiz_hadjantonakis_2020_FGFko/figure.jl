using Plots
gr()

using LaTeXStrings

plot_font = "DejaVu Sans"
default(fontfamily=plot_font,
        linewidth=3, label=nothing, grid=false,
        xtickfont=font(18), 
        ytickfont=font(18), 
        guidefont=font(18),
        titlefont=font(18),
        legendfontsize=16)

include("../../types/types.jl")
include("../../utils/agentsim_utils.jl")
include("../../utils/general_utils.jl")
include("../../utils/analysis_utils.jl")
include("constants/constants_mechanical.jl")
include("constants/constants_biochemical.jl")
include("utils/model_utils.jl")

### SIMULATION ###
@time _ = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
                    , h=0.1    #integration time-step
                    , mh=1000000000.0   #data measuring time-step
                    , Fth=10000.0
                    , comext=0.0
                    , comKO = false);

@time simresults = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
                          , h=0.001  #integration time-step
                          , mh=0.5   #data measuring time-step
                          , Fth=Inf
                          , comext=0.0
                          , comKO=false);


@time results = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
, h=0.0005  #integration time-step
, mh=0.001   #data measuring time-step
, Fth=Inf
, comext=0.0
, comKO=false);

rr = findfirst.(isequal.(unique(results.totals)), [results.totals])
N=20
tots = results.totals[rr]
FDPmat  = zeros(N,length(tots))
FEPImat = zeros(N,length(tots))
FPREmat = zeros(N,length(tots))
for n=1:N
@time _results = agentsimICM( saiz_hadjantonakis_2020, ["NANOG", "GATA6", "FGF"], "FGF"
                        , h=0.001  #integration time-step
                        , mh=0.001   #data measuring time-step
                        , Fth=Inf
                        , comext=0.0
                        , comKO=false);
    _rr = findfirst.(isequal.(unique(_results.totals)), [_results.totals])
    _tots = _results.totals[_rr]
    FDPmat[n,:]  .= _results.fDP[_rr]./_tots
    FEPImat[n,:] .= _results.fEPI[_rr]./_tots
    FPREmat[n,:] .= _results.fPRE[_rr]./_tots
end

FDP  = mean(FDPmat, dims=1)[:].*100
FEPI = mean(FEPImat, dims=1)[:].*100
FPRE = mean(FPREmat, dims=1)[:].*100
FDPstd = std(FDPmat, dims=1)[:].*100
FEPIstd = std(FEPImat, dims=1)[:].*100
FPREstd = std(FPREmat, dims=1)[:].*100

# ### END SIMULATION ###

### Plotting ###

l = grid(1, 3)
p2=plot()
p3=plot(right_margin=20Plots.mm)
p31 = twinx()
p4=plot(right_margin=20Plots.mm)
p41 = twinx()
ax = [p2,p3,p4]

rang = collect(Int64, range(1,step=1,stop=Nmax))
k=[0.0]

clrs = [
    RGB(0.0, 0.8, 0.0), 
    RGB(0.8, 0.0, 0.8), 
    RGB(101.0/230, 67/230, 33/230),
    RGB(0.8, 0.8, 0.0), 
    RGB(0.5, 0.5, 0.5)]

comvaridx = findfirst(simresults.comvarname .== simresults.varsnames)
comvarname = simresults.comvarname
global max_val = 0.0
global min_val = Inf

cells = [1,2,3]
all=false

jj1=3
jj2=2

if all 
    start = 1 
else 
    start = find_start(simresults.totals, N_start)
end
if simresults.CFATES[jj1,end] == 1
    FC = RGBA(0.87, 1.0, 0.87, 1.0)
    tit="time series Epi"
elseif simresults.CFATES[jj1,end] == 2
    FC = RGBA(1.0, 0.97, 1.0, 1.0)
    tit="time series PrE"
end
for vid in eachindex(simresults.vars)
    global max_val = max(max_val, maximum(simresults.vars[vid][jj1,start:end]))
    global min_val = min(min_val, minimum(simresults.vars[vid][jj1,start:end]))
    plot!(ax[2], simresults.times[start:end], simresults.vars[vid][jj1,start:end], label=simresults.varsnames[vid], color=clrs[vid], background_color_inside=FC, background_color_outside="white", legend=false)
end

plot(ax[2], simresults.times[start:end], simresults.comvar[jj1,start:end], label="$comvarname received", color=clrs[comvaridx],linestyle=:dash)
title!(ax[2], tit)
ylabel!(ax[2],"variables (a. u.)")
xlabel!(ax[2], "time")
xticks!(ax[2], [0], [""])

plot!(p31, simresults.times[start:end], simresults.totals[start:end], color="black", label="cell number")
ylims!(p31, N_start, Nmax)
ylabel!(p31, "cell number")

if simresults.CFATES[jj2,end] == 1
    FC = RGBA(0.87, 1.0, 0.87, 1.0)
    tit="time series Epi"
elseif simresults.CFATES[jj2,end] == 2
    FC = RGBA(1.0, 0.97, 1.0, 1.0)
    tit="time series PrE"
end
for vid in eachindex(simresults.vars)
    global max_val = max(max_val, maximum(simresults.vars[vid][jj2,start:end]))
    global min_val = min(min_val, minimum(simresults.vars[vid][jj2,start:end]))
    plot!(ax[3], simresults.times[start:end], simresults.vars[vid][jj2,start:end], label=simresults.varsnames[vid], color=clrs[vid], background_color_inside=FC, background_color_outside="white", legend=false)
end
plot(ax[3], simresults.times[start:end], simresults.comvar[jj2,start:end], label="$comvarname received", color=clrs[comvaridx],linestyle=:dash)
title!(ax[3], tit)
ylabel!(ax[3],"variables (a. u.)")
xlabel!(ax[3], "time")
xticks!(ax[3], [0], [""])

plot!(p41, simresults.times[start:end], simresults.totals[start:end], color="black", label="cell number")
ylims!(p41, N_start, Nmax)
ylabel!(p41, "cell number")


plot!(ax[1], tots, FDP, ribbon = FDPstd, color=RGBA(0.5, 0.5, 0.5, 0.3), label="DP", linewidth=5)
plot!(ax[1], tots, FEPI, ribbon = FEPIstd, color=RGBA(0.0, 0.8, 0.0, 0.3), label="Epi", linewidth=5)
plot!(ax[1], tots, FPRE, ribbon = FPREstd, color=RGBA(0.8, 0.0, 0.8), label="PrE", linewidth=5)

ylabel!(ax[1], join([L"\%", " of ICM"]))
xlabel!(ax[1], "cell number")
ylims!(ax[1], -1,101)

p = plot(p2,p3,p4,
 layout=l, 
 size=(1600,500), 
 bottom_margin=12Plots.mm, 
 top_margin=12Plots.mm,
#  right_margin=5Plots.mm,
 left_margin=5Plots.mm
 )

 # Save figure
cwd = pwd()
foldername = basename(cwd)
basepath = dirname(dirname(cwd))
save_dir = joinpath(basepath, "results", foldername)
mkpath(save_dir)

savepath = joinpath(save_dir, "fate_progression.pdf")
savefig(p, savepath)
