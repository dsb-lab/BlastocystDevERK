include("/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models_EmbryonicDev/NGF/utils/utils.jl")
include("/home/pablo/Desktop/PhD/projects/BlastocystDev/ICMModels/models_EmbryonicDev/NGF/utils/general_utils.jl")


using LaTeXStrings
using BifurcationKit
using Revise, ForwardDiff, Parameters, Plots, Setfield, LinearAlgebra
const BK = BifurcationKit
using Plots.PlotMeasures

# Save figure
cwd = pwd()
foldername = basename(cwd)
basepath = dirname(dirname(cwd))
save_dir = joinpath(basepath, "results", foldername)
mkpath(save_dir)


# sup norm
norminf(x) = norm(x, Inf)

# vector field
function model!(dz, z, p, t)
	@unpack dn, Kn, an, ag, dg, Kg, coefn, coefm, FGFmedium = p
    NANOG, GATA6 = z
    dz[2] = ag.*FGFmedium./(1.0 .+ (NANOG./Kn).^coefn) .- dg.*GATA6
    dz[1] = an./(1.0 .+ (FGFmedium.*GATA6./Kg).^coefm) .- dn.*NANOG
	dz
end

# out of place method
model(z, p) = model!(similar(z), z, p, 0)

# we group the differentials together
dmodel(z,p) = ForwardDiff.jacobian(x -> model(x,p), z)

# parameter values
alpha=10.0
K=0.9
par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 2.0, coefm = 2.0, FGFmedium=3.5)

# initial condition
z0 = [0.1, 16.0]

# continuation options
opts_FGF = ContinuationPar(pMin = 0.0, pMax = 15.0,
	# parameters to have a smooth result
	dsmin=0.0001, ds = 0.001, dsmax = 0.001,
	# this is to detect bifurcation points precisely with bisection
	detectBifurcation = 3, plotEveryStep=1,
	# Optional: bisection options for locating bifurcations
	nInversion = 8, maxBisectionSteps = 100, nev = 3)

# continuation of equilibria
br1NANOG, = continuation(model, dmodel, z0, par_tm, (@lens _.FGFmedium), opts_FGF;
    bothside=true,
	recordFromSolution = (x, p) -> (NANOG = x[1], GATA6 = x[2]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

br1GATA6, = continuation(model, dmodel, z0, par_tm, (@lens _.FGFmedium), opts_FGF;
    bothside=true,
	recordFromSolution = (x, p) -> (GATA6 = x[2], NANOG = x[1]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)
# continuation options
opts_n = ContinuationPar(pMin = 1.0, pMax = 4.0,
    # parameters to have a smooth result
    dsmin=0.000001, ds = 0.00001, dsmax = 0.0001,
    maxSteps=200000,
    # this is to detect bifurcation points precisely with bisection
    detectBifurcation = 2, plotEveryStep=1,
    # Optional: bisection options for locating bifurcations
    nInversion = 8, maxBisectionSteps = 10000, nev = 3)

par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 1.0, coefm = 2.0, FGFmedium=3.5)

br21NANOG, = continuation(model, dmodel, [0.1,16], par_tm, (@lens _.coefn), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (NANOG = x[1], GATA6 = x[2]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 1.5, coefm = 2.0, FGFmedium=3.5)
br22NANOG, = continuation(model, dmodel, [5.0, 5.0], par_tm, (@lens _.coefn), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (NANOG = x[1], GATA6 = x[2]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 1.0, coefm = 2.0, FGFmedium=3.5)
br21GATA6, = continuation(model, dmodel, z0, par_tm, (@lens _.coefn), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (GATA6 = x[2], NANOG = x[1]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 1.5, coefm = 2.0, FGFmedium=3.5)
br22GATA6, = continuation(model, dmodel, [5.0, 5.0], par_tm, (@lens _.coefn), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (GATA6 = x[2], NANOG = x[1]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

# Bifurcations for parameter m
par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 2.0, coefm = 1.0, FGFmedium=3.5)
brm1NANOG, = continuation(model, dmodel, [16, 0.1], par_tm, (@lens _.coefm), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (NANOG = x[1], GATA6 = x[2]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 2.0, coefm = 1.5, FGFmedium=3.5)
brm2NANOG, = continuation(model, dmodel, [2.0, 3.0], par_tm, (@lens _.coefm), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (NANOG = x[1], GATA6 = x[2]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 2.0, coefm = 1.0, FGFmedium=3.5)
brm1GATA6, = continuation(model, dmodel, [16, 0.1], par_tm, (@lens _.coefm), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (GATA6 = x[2], NANOG = x[1]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)

par_tm = (dn = 1.0, Kn = 1.0, an = alpha*1.0*1.0, ag = 10.0, dg = 2.0, Kg = 5.0, coefn = 2.0, coefm = 1.5, FGFmedium=3.5)
brm2GATA6, = continuation(model, dmodel, [5.0, 5.0], par_tm, (@lens _.coefm), opts_n;
    bothside=true,
	recordFromSolution = (x, p) -> (GATA6 = x[2], NANOG = x[1]),
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)


plot_font = "DejaVu Sans"
default(fontfamily=plot_font,
        linewidth=3, label=nothing, grid=false,
        xtickfont=font(18), 
        ytickfont=font(18), 
        guidefont=font(18),
        titlefont=font(18),
        legendfontsize=16)

p1NANOG = plot(br1NANOG, plotfold=true, markersize=5, legend=false, linecolor="red", markercolor="black")
p1GATA6 = plot(br1GATA6, plotfold=true, markersize=5, legend=:topright, linecolor="blue", markercolor="black")
name=join(["BD_FGFmedium", ".png"])
#savefig(plot(p1NANOG, p1GATA6, layout=2, size=(1000,500)), "../../figures/NGF/FGFko/bifurcation/$name")

p2NANOG = plot(brm1NANOG, plotfold=true, markersize=5, legend=false, linecolor="red", markercolor="black")
plot!(brm2NANOG, plotfold=true, markersize=5, legend=false, linecolor="red", markercolor="black")
p2GATA6 = plot(brm1GATA6, plotfold=true, markersize=5, legend=false, linecolor="blue", markercolor="black")
plot!(brm2GATA6, plotfold=true, markersize=5, legend=:topright, linecolor="blue", markercolor="black")
name=join(["BD_coefm", ".png"])
#savefig(plot(p2NANOG, p2GATA6, layout=2, size=(1000,500)), "../../figures/NGF/FGFko/bifurcation/$name")

p2NANOG = plot(br21NANOG, plotfold=true, markersize=5, legend=false, linecolor="red", markercolor="black")
plot!(br22NANOG, plotfold=true, markersize=5, legend=false, linecolor="red", markercolor="black")
p2GATA6 = plot(br21GATA6, plotfold=true, markersize=5, legend=false, linecolor="blue", markercolor="black")
plot!(br22GATA6, plotfold=true, markersize=5, legend=:topright, linecolor="blue", markercolor="black")
name=join(["BD_coefn", ".png"])
#savefig(plot(p2NANOG, p2GATA6, layout=2, size=(1000,500)), "../../figures/NGF/FGFko/bifurcation/$name")

# continuation options
opts = ContinuationPar(pMin = 0.1, pMax = 6.0,
    # parameters to have a smooth result
    dsmin=0.000001, ds = 0.00001, dsmax = 0.0001,
    maxSteps=200000,
    # this is to detect bifurcation points precisely with bisection
    detectBifurcation = 2, plotEveryStep=1,
    # Optional: bisection options for locating bifurcations
    nInversion = 8, maxBisectionSteps = 10000, nev = 3)
br31, = continuation(model, dmodel, br1NANOG, 1, (@lens _.coefm), opts;
	bothside=true,
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)
br32, = continuation(model, dmodel, br1NANOG, 2, (@lens _.coefm), opts;
	bothside=true,
	tangentAlgo = BorderedPred(),
	plot = false, normC = norminf)
P = plot(br31, plotfold=true, markersize=0, legend=:topleft, linewidth=5,linecolor="black", size=(1000,500), xtickfontsize=25, ytickfontsize=25,xguidefontsize=25, yguidefontsize=25)
plot!(br32, plotfold=true, markersize=0, legend=:topleft, linewitdh=5, linecolor="black", bottom_margin = 10mm, left_margin=10mm)
scatter!([1.5, 1.5, 1.5, 4.0, 4.0, 4.0], [2.0, 4.0, 6.0, 2.0, 4.0, 6.0], markersize=10, color="green", leg=false)
ylabel!(L"$F_m$")
xlabel!(L"$m$")

savefig(P, join([save_dir, "/codim2.svg"]))

#xlims!(0.0, 10)